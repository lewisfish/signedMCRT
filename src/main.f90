program mcpolar

#ifdef MPI
    use mpi_f08
#endif

#ifdef _OPENMP
    use omp_lib
#endif

!shared data
use iarray
use random,    only : ran2, init_rng, ranu
use constants, only : fileplace, wp

!subroutines
use sdfs
use subs
use inttau2
use hitStack
use parse_mod
use photonMod
use stokes_mod
use writer_mod
use detector_mod
use vector_class
use sim_state_mod, only : state
use string_utils,  only : str
use utils,         only : pbar, get_time, print_time

!external deps
use tev_mod
use tomlf
implicit none

type(container),   allocatable :: array(:)
type(dect_array),  allocatable :: dects(:)
character(len=64), allocatable :: args(:)
type(toml_table)  :: dict
type(tevipc)      :: tev
type(photon)      :: packet
type(pbar)        :: bar
type(hit_t)       :: hpoint
type(vector)      :: dir

real(kind=wp), allocatable :: distances(:), image(:,:,:)
real(kind=wp) :: ran, start, time_taken, nscatt
integer       :: i, j, num_args

! mpi/mp variables
integer       :: id, numproc
real(kind=wp) :: nscattGLOBAL!, chance, threshold, absorb

! chance = 1._wp/10._wp
! threshold = 1e-6_wp

dict = toml_table()
num_args = command_argument_count()
if(num_args == 0)then
    allocate(args(1))
    args(1) = "input.toml"
else
    allocate(args(num_args))
    do i = 1, num_args
        call get_command_argument(i, args(i))
    end do
end if


call parse_params("res/"//trim(args(1)), packet, dects, dict)
allocate(image(state%grid%nxg,state%grid%nzg,1))

if(state%tev)then
    !init TEV link
    tev = tevipc()
    call tev%close_image(state%experiment)
    call tev%create_image(state%experiment, state%grid%nxg, state%grid%nzg, ["I", "J", "K"], .true.)
end if

nscatt = 0._wp
! call init_rng(spread(state%iseed+0, 1, 8), fwd=.true.)

call setup_simulation(array, dict)
! render geometry to voxel format for debugging
if(state%render_geom)then
    call render(array, vector(state%grid%xmax, state%grid%ymax, state%grid%zmax), state%render_size, fname=state%renderfile)
end if

allocate(distances(size(array)))

start = get_time()
id = 0
!init MPI
#ifdef MPI
    call MPI_init()
    call MPI_Comm_size(MPI_COMM_WORLD, numproc)
    call MPI_Comm_rank(MPI_COMM_WORLD, id)
#endif


if(id == 0)then
   print*,'# of photons to run',state%nphotons
end if

#ifdef _OPENMP
!is state%seed private, i dont think so...
!$omp parallel default(none) shared(dict, array, numproc, start, state, bar, jmean, tev, dects)&
!$omp& private(ran, id, distances, image, dir, hpoint) reduction(+:nscatt) firstprivate(packet)
    numproc = omp_get_num_threads()
    id = omp_get_thread_num()
    if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
#elif MPI
    !nothing
#else
    numproc = 1
    id = 0
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

    packet%id = id
    distances = 0._wp
    do i = 1, size(distances)
        distances(i) = array(i)%p%evaluate(packet%pos)
        if(distances(i) > 0._wp)distances(i)=-999.0_wp
    end do
    packet%layer=minloc(abs(distances),dim=1)

    ! Find scattering location
    call tauint2(state%grid, packet, array)

    ! dir = vector(packet%nxp, packet%nyp, packet%nzp)
    ! hpoint = hit_t(packet%pos, dir, 1._wp, packet%layer)
    ! do i = 1, size(dects)
    !     call dects(i)%p%record_hit(hpoint)
    ! end do

    ! Photon scatters in grid until it exits (tflag=TRUE)
    do while(.not. packet%tflag)
        ran = ran2()
        if(ran < array(packet%layer)%p%albedo)then!interacts with tissue
            call stokes(packet, array(packet%layer)%p%hgg, array(packet%layer)%p%g2)
            nscatt = nscatt + 1
        else
            packet%tflag = .true.
            exit
        end if
        ! !Find next scattering location
        call tauint2(state%grid, packet, array)
    end do

    dir = vector(packet%nxp, packet%nyp, packet%nzp)
    hpoint = hit_t(packet%pos, dir, 1._wp, packet%layer)
    do i = 1, size(dects)
        call dects(i)%p%record_hit(hpoint)
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
    open(newunit=j,file=trim(fileplace)//"camera.dat",form="unformatted",access="stream")
    associate(x => dects(1)%p)
        select type(x)
            type is(camera)
                write(j)x%data(:x%nbinsX-1,:x%nbinsY-1)
        end select
    end associate
    close(j)
    print*,'write done'
end if
call write_detected_photons(dects)

time_taken = get_time() - start
call print_time(time_taken, id)
#ifdef MPI
    call MPI_Finalize()
#endif
end program mcpolar