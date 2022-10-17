module kernels

    implicit none
    
    private
    public :: weight_scatter, pathlength_scatter

contains
!###############################################################################
!                   KERNELS
    subroutine weight_scatter()

        error stop "Not yet implmented!"

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
        use random,        only : ran2, init_rng
        use sdfs,          only : container
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
        character(*), intent(in) :: input_file
        
        integer :: numproc, id, j, i
        ! type(history_stack_t) :: history
        type(pbar)        :: bar
        type(photon)      :: packet
        type(toml_table)  :: dict
        real(kind=wp), allocatable :: distances(:), image(:,:,:)
        type(hit_t)       :: hpoint
        type(vector)      :: dir
        type(dect_array),  allocatable :: dects(:)
        type(container),   allocatable :: array(:)
        real(kind=wp) :: ran, nscatt, start
        type(tevipc)      :: tev
    

        call setup(input_file, tev, dects, array, packet, dict, distances, image, nscatt, start)

#ifdef _OPENMP
        !is state%seed private, i dont think so...
        !$omp parallel default(none) shared(dict, array, numproc, start, state, bar, jmean, tev, dects)&
        !$omp& private(ran, id, distances, image, dir, hpoint) reduction(+:nscatt) firstprivate(packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        ! history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        numproc = 1
        id = 0
        ! history = history_stack_t(state%historyFilename, id)
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
            call packet%emit(dict)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%p%evaluate(packet%pos)
                if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=minloc(abs(distances),dim=1)

            ! Find scattering location
            call tauint2(state%grid, packet, array)

            dir = vector(packet%nxp, packet%nyp, packet%nzp)
            hpoint = hit_t(packet%pos, dir, packet%step, packet%layer)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint)!, history)
            end do

            do while(.not. packet%tflag)
                ran = ran2()
                if(ran < array(packet%layer)%p%albedo)then!interacts with tissue
                    call packet%scatter(array(packet%layer)%p%hgg, array(packet%layer)%p%g2)
                    ! call stokes(packet, array(packet%layer)%p%hgg, array(packet%layer)%p%g2)
                    nscatt = nscatt + 1
                    packet%step = packet%step + 1
                else
                    packet%tflag = .true.
                    exit
                end if
                ! !Find next scattering location
                call tauint2(state%grid, packet, array)
                ! call history%push(vec4(packet%pos, packet%step))
            end do

            dir = vector(packet%nxp, packet%nyp, packet%nzp)
            hpoint = hit_t(packet%pos, dir, packet%step, packet%layer)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint)!, history)
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

    call finalise(dict, dects, nscatt, start)
    end subroutine pathlength_scatter


!####################################################################################################
!                           Setup and break down helper routines
    subroutine setup(input_file, tev, dects, array, packet, dict, distances, image, nscatt, start)

        ! !shared data
        use iarray
        use constants, only : wp
        
        ! !subroutines
        use detector_mod,  only : dect_array
        use parse_mod,     only : parse_params
        use photonMod,     only : photon
        use random,        only : init_rng
        use sdfs,          only : container, render
        use sim_state_mod, only : state
        use subs,          only : setup_simulation
        use utils,         only : get_time, print_time
        use vector_class,  only : vector
        
        ! !external deps
        use tev_mod, only : tevipc, tev_init
        use tomlf,   only : toml_table
        
        character(*), intent(in) :: input_file

        type(container),  allocatable, intent(out) :: array(:)
        type(dect_array), allocatable, intent(out) :: dects(:)
        type(toml_table), intent(out)  :: dict
        type(tevipc),     intent(out)  :: tev
        type(photon),     intent(out)  :: packet
        real(kind=wp), allocatable, intent(out) :: distances(:), image(:,:,:)
        real(kind=wp), intent(out) :: nscatt
        real(kind=wp), intent(out) :: start
        
        ! mpi/mp variables
        integer       :: id
        real(kind=wp) :: chance, threshold
        
        chance = 1._wp/10._wp
        threshold = 1e-6_wp
        
        dict = toml_table()
        call parse_params("res/"//trim(input_file), packet, dects, dict)
        allocate(image(state%grid%nxg,state%grid%nzg,1))
        
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
            call render(array, vector(state%grid%xmax, state%grid%ymax, state%grid%zmax), state%render_size, fname=state%renderfile)
        end if
        
        allocate(distances(size(array)))
        
        start = get_time()
        id = 0            
        
        if(id == 0)then
           print*,'# of photons to run',state%nphotons
        end if
end subroutine setup

subroutine finalise(dict, dects, nscatt, start)

    use constants,     only : wp, fileplace
    use iarray,        only : jmean, jmeanGLOBAL
    use sim_state_mod, only : state
    
    ! use historyStack, only : history_stack_t
    use detector_mod, only : dect_array
    use writer_mod,   only : normalise_fluence, write_fluence, write_detected_photons
    
    use utils, only : get_time, print_time, str
    use tomlf, only : toml_table, set_value

    real(kind=wp),         intent(in) :: nscatt, start
    type(dect_array),      intent(in) :: dects(:)
    ! type(history_stack_t), intent(in) :: history
    type(toml_table),      intent(inout) :: dict

    integer :: id, numproc!, i
    real(kind=wp) :: nscattGLOBAL, time_taken

    id = 0
    numproc = 1

#ifdef MPI
    ! collate fluence from all processes
    call mpi_reduce(jmean, jmeanGLOBAL, size(jmean),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD)
#else
    jmeanGLOBAL = jmean
    nscattGLOBAL = nscatt
#endif

    if(id == 0)then
#ifdef _OPENMP
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons)
#else
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons*numproc)
#endif
        !write out files
        !create dict to store metadata and nrrd hdr info
        call set_value(dict, "grid_data", "fluence map")
        call set_value(dict, "real_size", str(state%grid%xmax,7)//" "//str(state%grid%ymax,7)//" "//str(state%grid%zmax,7))
        call set_value(dict, "nphotons", state%nphotons)
        call set_value(dict, "source", state%source)
        call set_value(dict, "experiment", state%experiment)

        jmeanGLOBAL = normalise_fluence(state%grid, jmeanGLOBAL, state%nphotons)
        call write_fluence(jmeanGLOBAL, trim(fileplace)//"jmean/"//state%outfile, dict)
    end if
    !write out detected photons
    if(size(dects) > 0)then
        call write_detected_photons(dects)
        ! block
        !     logical :: mask(size(dects))
        !     do i = 1, size(dects)
        !         mask(i) = dects(i)%p%trackHistory
        !     end do
        !     if(any(mask))call history%finish(numproc)
        ! end block
    end if

    time_taken = get_time() - start
    call print_time(time_taken, id)
#ifdef MPI
    call MPI_Finalize()
#endif

end subroutine finalise
end module kernels