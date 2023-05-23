module kernels

    implicit none
    
    private
    public :: weight_scatter, pathlength_scatter, test_kernel

contains
!###############################################################################
!                   KERNELS
    subroutine weight_scatter(input_file)

        !Shared data
        use iarray
        use constants, only : wp, CHANCE, THRESHOLD

        !subroutines
        use detector_mod,  only : dect_array, hit_t
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2, update_voxels
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vec4_class,    only : vec4
        use vector_class,  only : vector
        
        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in)   :: input_file
        
        integer                        :: numproc, id, j, i
        type(history_stack_t)          :: history
        type(pbar)                     :: bar
        type(photon)                   :: packet
        type(toml_table)               :: dict
        real(kind=wp),     allocatable :: distances(:), image(:,:,:)
        type(hit_t)                    :: hpoint
        type(vector)                   :: dir
        type(dect_array),  allocatable :: dects(:)
        type(sdf),         allocatable :: array(:)
        real(kind=wp)                  :: nscatt, start, weight_absorb
        type(tevipc)                   :: tev
        integer                        :: celli, cellj, cellk
        type(spectrum_t)               :: spectrum

        call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

#ifdef _OPENMP
        !is state%seed private, i dont think so...
        !$omp parallel default(none) shared(dict, array, numproc, start, state, bar, jmean, tev, dects, spectrum)&
        !$omp& private(id, distances, image, dir, hpoint, history, weight_absorb, cellk, cellj, celli) &
        !$omp& reduction(+:nscatt) firstprivate(packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        numproc = 1
        id = 0
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#endif
        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

        ! set seed for rnd generator. id to change seed for each process
        call init_rng(spread(state%iseed+id, 1, 8), fwd=.true.)

        bar = pbar(state%nphotons/ 10)
        !$OMP BARRIER
        !$OMP do
        !loop over photons 
        do j = 1, state%nphotons
            if(mod(j, 10) == 0)call bar%progress()

            ! Release photon from point source
            call packet%emit(spectrum, dict)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=minloc(abs(distances),dim=1)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ! Find scattering location
            call tauint2(state%grid, packet, array)

            do while(.not. packet%tflag)
                if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
                
                weight_absorb = packet%weight * (1._wp - array(packet%layer)%getAlbedo())

                packet%weight = packet%weight - weight_absorb
                call update_voxels(state%grid, &
                packet%pos + vector(state%grid%xmax, state%grid%ymax, state%grid%zmax), celli, cellj, cellk)

                if(celli < 1)then
                    packet%tflag = .true.
                    exit
                end if
                if(cellj < 1)then
                    packet%tflag = .true.
                    exit
                end if
                if(cellk < 1)then
                    packet%tflag = .true.
                    exit
                end if
                !$omp atomic
                jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + weight_absorb
                call packet%scatter(array(packet%layer)%gethgg(), array(packet%layer)%getg2(), dects)
                if(packet%weight < THRESHOLD)then
                    if(ran2() < CHANCE)then
                        packet%weight = packet%weight / CHANCE
                    else
                        packet%tflag = .true.
                        exit
                    end if
                end if

                ! !Find next scattering location
                call tauint2(state%grid, packet, array)
            end do

            dir = vector(packet%nxp, packet%nyp, packet%nzp)
            hpoint = hit_t(packet%pos, dir, packet%weight, packet%layer)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint, history)
            end do

            if(id == 0 .and. mod(j,1000) == 0)then
                if(state%tev)then
!$omp critical
                    image = reshape(jmean(:,100:100,:), [state%grid%nxg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["I"], 0, 0, .false., .false.)

                    image = reshape(jmean(100:100,:,:), [state%grid%nyg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["J"], 0, 0, .false., .false.)

                    image = reshape(jmean(:,:,100:100), [state%grid%nxg,state%grid%nyg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["K"], 0, 0, .false., .false.)
!$omp end critical
                end if
            end if
        end do

#ifdef _OPENMP
!$OMP end  do
!$OMP end parallel
#endif
    call finalise(dict, dects, nscatt, start, history)
    end subroutine weight_scatter

    subroutine pathlength_scatter(input_file)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detector_mod,  only : dect_array, hit_t
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vec4_class,    only : vec4
        use vector_class,  only : vector
        use piecewiseMod

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: numproc, id, j, i
        type(history_stack_t)         :: history
        type(pbar)                    :: bar
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(hit_t)                   :: hpoint
        type(vector)                  :: dir
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: ran, nscatt, start
        type(tevipc)                  :: tev
        type(seq)                     :: seqs(2)
        type(spectrum_t)              :: spectrum

        call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

#ifdef _OPENMP
        !is state%seed private, i dont think so...
        !$omp parallel default(none) shared(dict, array, numproc, start, state, bar, jmean, phasor, tev, dects, spectrum)&
        !$omp& private(ran, id, distances, image, dir, hpoint, history, seqs) reduction(+:nscatt) firstprivate(packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        numproc = 1
        id = 0
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#endif
        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

        ! set seed for rnd generator. id to change seed for each process
        call init_rng(spread(state%iseed+id, 1, 8), fwd=.true.)
        seqs = [seq((id+1)*(state%nphotons/numproc), 2),&
                seq((id+1)*(state%nphotons/numproc), 3)]

        bar = pbar(state%nphotons/ 10)
        !$OMP BARRIER
        !$OMP do
        !loop over photons 
        do j = 1, state%nphotons
            if(mod(j, 10) == 0)call bar%progress()

            ! Release photon from point source
            call packet%emit(spectrum, dict, seqs)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=minloc(abs(distances),dim=1)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ! Find scattering location
            call tauint2(state%grid, packet, array)

            do while(.not. packet%tflag)
                if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
                ran = ran2()

                if(ran < array(packet%layer)%getAlbedo())then!interacts with tissue
                    call packet%scatter(array(packet%layer)%gethgg(), &
                                        array(packet%layer)%getg2(), dects)
                    nscatt = nscatt + 1
                    packet%step = packet%step + 1
                else
                    packet%tflag = .true.
                    exit
                end if
                ! !Find next scattering location
                call tauint2(state%grid, packet, array)
            end do

            dir = vector(packet%nxp, packet%nyp, packet%nzp)
            hpoint = hit_t(packet%pos, dir, sqrt(packet%pos%x**2+packet%pos%y**2), packet%layer)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint, history)
            end do
            if(id == 0 .and. mod(j,1000) == 0)then
                if(state%tev)then
!$omp critical
                    image = reshape(jmean(:,100:100,:), [state%grid%nxg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["I"], 0, 0, .false., .false.)

                    image = reshape(phasor(100:100,:,:), [state%grid%nyg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["J"], 0, 0, .false., .false.)

                    image = reshape(phasor(:,:,100:100), [state%grid%nxg,state%grid%nyg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["K"], 0, 0, .false., .false.)
!$omp end critical
                end if
            end if
        end do

#ifdef _OPENMP
!$OMP end  do
!$OMP end parallel
#endif
    call finalise(dict, dects, nscatt, start, history)
    end subroutine pathlength_scatter

    subroutine test_kernel(input_file, end_early)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detector_mod,  only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vector_class,  only : vector
        
        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer :: numproc, id, j, i
        type(history_stack_t) :: history
        ! type(pbar)        :: bar
        type(photon)      :: packet
        type(toml_table)  :: dict
        real(kind=wp), allocatable :: distances(:), image(:,:,:)
        type(dect_array),  allocatable :: dects(:)

        type(sdf),   allocatable :: array(:)
        real(kind=wp) :: ran, nscatt, start
        type(tevipc)      :: tev
        type(vector)  :: pos(4), pos2(4)
        logical, intent(in) :: end_early
        type(spectrum_t) :: spectrum

        pos = vector(0.0_wp, 0.0_wp, 0.0_wp)
        pos2 = vector(0.0_wp, 0.0_wp, 0.0_wp)
        call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

        numproc = 1
        id = 0

        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

        ! set seed for rnd generator. id to change seed for each process
        call init_rng(spread(state%iseed+id, 1, 8), fwd=.true.)

        ! bar = pbar(state%nphotons/ 10)
        !loop over photons 
        do j = 1, state%nphotons
            ! if(mod(j, 10) == 0)call bar%progress()

            ! Release photon from point source
            call packet%emit(spectrum, dict)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=minloc(abs(distances),dim=1)

            ! Find scattering location
            call tauint2(state%grid, packet, array)


            do while(.not. packet%tflag)
                ran = ran2()
                if(ran < array(packet%layer)%getalbedo())then!interacts with tissue
                    call packet%scatter(array(packet%layer)%gethgg(), &
                                        array(packet%layer)%getg2())
                    nscatt = nscatt + 1
                    packet%step = packet%step + 1
                    if(packet%step == 1)then
                        pos(1) = pos(1) + packet%pos
                        pos2(1) = pos2(1) + packet%pos**2
                    elseif(packet%step == 2)then
                        pos(2) = pos(2) + packet%pos
                        pos2(2) = pos2(2) + packet%pos**2
                    elseif(packet%step == 3)then
                        pos(3) = pos(3) + packet%pos
                        pos2(3) = pos2(3) + packet%pos**2
                    elseif(packet%step == 4)then
                        pos(4) = pos(4) + packet%pos
                        pos2(4) = pos2(4) + packet%pos**2
                    else
                        if(end_early)packet%tflag = .true.
                    end if
                else
                    packet%tflag = .true.
                    exit
                end if
                ! !Find next scattering location
                call tauint2(state%grid, packet, array)
            end do
        end do

    open(newunit=j,file="positions.dat")
    do i = 1, 4
        write(j,*)10.*pos(i)%x/state%nphotons,10.*pos(i)%y/state%nphotons,10.*pos(i)%z/state%nphotons
    end do
    do i = 1,4
        write(j,*)100.*pos2(i)%x/state%nphotons,100.*pos2(i)%y/state%nphotons,100.*pos2(i)%z/state%nphotons
    end do
    close(j)
    call finalise(dict, dects, nscatt, start, history)
    end subroutine test_kernel


!####################################################################################################
!                           Setup and break down helper routines
    subroutine setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

        ! !shared data
        use iarray
        use constants, only : wp
        
        ! !subroutines
        use detector_mod,  only : dect_array
        use parse_mod,     only : parse_params
        use photonMod,     only : photon
        use random,        only : init_rng
        use piecewiseMod

        use sdfs,          only : sdf, render
        use sim_state_mod, only : state
        use subs,          only : setup_simulation, directory
        use utils,         only : get_time, print_time, str
        use vector_class,  only : vector
        ! !external deps
        use tev_mod, only : tevipc, tev_init
        use tomlf,   only : toml_table
        
        character(*),                  intent(in)  :: input_file
        type(sdf),        allocatable, intent(out) :: array(:)
        type(dect_array), allocatable, intent(out) :: dects(:)
        type(toml_table),              intent(out) :: dict
        type(tevipc),                  intent(out) :: tev
        type(photon),                  intent(out) :: packet
        real(kind=wp),    allocatable, intent(out) :: distances(:), image(:,:,:)
        real(kind=wp),                 intent(out) :: nscatt, start
        type(spectrum_t),              intent(out) :: spectrum
        
        ! mpi/mp variables
        integer       :: id
        real(kind=wp) :: chance, threshold
        
        chance = 1._wp/10._wp
        threshold = 1e-6_wp
        
        call directory()

        dict = toml_table()
        call parse_params("res/"//trim(input_file), packet, dects, spectrum, dict)
        allocate(image(state%grid%nxg,state%grid%nzg,1))
        
        call display_settings(state, input_file, packet, "Pathlength")

        if(state%tev)then
            !init TEV link
            tev = tevipc()
            call tev%close_image(state%experiment)
            call tev%create_image(state%experiment, state%grid%nxg, state%grid%nzg, ["I", "J", "K"], .true.)
        end if

        nscatt = 0._wp
        call init_rng(spread(state%iseed+0, 1, 8), fwd=.true.)
        
        call setup_simulation(array, dict)
        ! render geometry to voxel format for debugging
        if(state%render_geom)then
            print*,"Rendering geometry to file"
            call render(array, state)
        end if
        
        allocate(distances(size(array)))
        
        start = get_time()
        id = 0            
        
        if(id == 0)then
           print*,'# of photons to run',state%nphotons
        end if
end subroutine setup

subroutine finalise(dict, dects, nscatt, start, history)

    use constants,     only : wp, fileplace
    use detector_mod,  only : dect_array
    use historyStack,  only : history_stack_t
    use iarray,        only : phasor, phasorGLOBAL, jmean, jmeanGLOBAL, absorb, absorbGLOBAL
    use sim_state_mod, only : state
    use subs,          only : dealloc_array
    use writer_mod,    only : normalise_fluence, write_data, write_detected_photons
    
    use utils, only : get_time, print_time, str
    use tomlf, only : toml_table, set_value

    real(kind=wp),         intent(in)    :: nscatt, start
    type(dect_array),      intent(in)    :: dects(:)
    type(history_stack_t), intent(in)    :: history
    type(toml_table),      intent(inout) :: dict

    integer       :: id, numproc, i
    real(kind=wp) :: nscattGLOBAL, time_taken

    id = 0
    numproc = 1

#ifdef MPI
    ! collate fluence from all processes
    call mpi_reduce(jmean, jmeanGLOBAL, size(jmean),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(absorb, absorbGLOBAL, size(absorb),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(phasor, phasorGLOBAL, size(phasor),MPI_DOUBLE_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD)
#else
    jmeanGLOBAL = jmean
    absorbGLOBAL = absorb
    phasorGLOBAL = phasor
    nscattGLOBAL = nscatt
#endif

    if(id == 0)then
#ifdef _OPENMP
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons)
#else
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons*numproc)
        ! for testing purposes
        open(newunit=i,file="nscatt.dat")
        write(i,*)nscattGLOBAL/(state%nphotons)
        close(i)
#endif
        !write out files
        !create dict to store metadata and nrrd hdr info
        call set_value(dict, "grid_data", "fluence map")
        call set_value(dict, "real_size", str(state%grid%xmax,7)//" "//str(state%grid%ymax,7)//" "//str(state%grid%zmax,7))
        call set_value(dict, "nphotons", state%nphotons)
        call set_value(dict, "source", state%source)
        call set_value(dict, "experiment", state%experiment)

        call normalise_fluence(state%grid, jmeanGLOBAL, state%nphotons)
        call write_data(jmeanGLOBAL, trim(fileplace)//"jmean/"//state%outfile, state, dict)
        ! if(state%absorb)call write_data(absorbGLOBAL, trim(fileplace)//"deposit/"//state%outfile_absorb, state, dict)
        !INTENSITY
        ! call write_data(abs(phasorGLOBAL)**2, trim(fileplace)//"phasor/"//state%outfile, state, dict)    
    end if
    !write out detected photons
    if(size(dects) > 0)then
        call write_detected_photons(dects)
        block
            logical :: mask(size(dects))
            do i = 1, size(dects)
                mask(i) = dects(i)%p%trackHistory
            end do
            if(state%trackHistory)call history%finish()
        end block
    end if

    time_taken = get_time() - start
    call print_time(time_taken, 4)
#ifdef MPI
    call MPI_Finalize()
#endif
    call dealloc_array()
end subroutine finalise

subroutine display_settings(state, input_file, packet, kernel_type)

    use sim_state_mod, only : settings_t
    use photonMod,     only : photon
    use utils,         only : str

    type(settings_t), intent(IN) :: state
    character(*),     intent(IN) :: input_file, kernel_type
    type(photon),     intent(IN) :: packet

    print*,repeat("#", 20)//" Settings "//repeat("#", 20)
    print*,"# Config file: ",trim(input_file),repeat(" ", 50-16-len(trim(input_file))),"#"
    print*,"# Using: "//trim(kernel_type)//"kernel"//repeat(" ", 50-16-len(kernel_type)),"#"
    print*,"# Light source: "//trim(state%source)//repeat(" ", 50-17-len(trim(state%source))),"#"
    if(state%source == "point")then
        print*,"# Light Source Position: ["//str(packet%pos%x,4)//", "//str(packet%pos%y,4)//", "//str(packet%pos%z,4)// &
                                        "]"//repeat(" ", 6)//"#"
    else
        print*,"# Light direction: ["//str(packet%nxp,4)//", "//str(packet%nyp,4)//", "//str(packet%nzp,4)// &
                                  "]"//repeat(" ", 12)//"#"
    end if
    print*,"# Geometry: "//trim(state%experiment)//repeat(" ", 50-13-len(trim(state%experiment))),"#"
    print*,"# Seed: "//str(state%iseed,9)//repeat(" ", 32)//"#"
    if(state%tev)then
        print*,"# Tev enabled!"//repeat(" ", 35)//"#"
    end if
    if(state%render_geom)then
        print*,"# Render geometry to file enabled!"//repeat(" ", 15)//"#"
    end if
    if(state%overwrite)then
        print*,"# Overwrite Enabled!",repeat(" ", 29)//"#"
    end if
    if(state%absorb)then
        print*,"# Energy absorbed will be written to file."//repeat(" ", 7)//"#"
    end if
    print*,repeat("#", 50)
    print*,new_line("a")

end subroutine display_settings
end module kernels