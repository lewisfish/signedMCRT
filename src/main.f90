program mcpolar

#ifdef MPI
    use mpi_f08
#endif

#ifdef _OPENMP
    use omp_lib
#endif

!shared data
use photonMod
use iarray
use random,    only : ran2, init_rng, ranu
use constants, only : fileplace, wp

!subroutines
use subs
use inttau2
use stokes_mod
use writer_mod
use vector_class
use sdfs
use utils, only : pbar, str
use parse_mod
use sim_state_mod, only : state

implicit none

type(photon)   :: packet
integer        :: j, i
real(kind=wp)  :: ran, start, time_taken, nscatt
type(pbar)     :: bar
real(kind=wp),   allocatable :: ds(:)
type(container), allocatable :: array(:)

! mpi/mp variables
integer       :: id, numproc
real(kind=wp) :: nscattGLOBAL, optprop(5), focus
type(dict_t)  :: dict

dict = dict_t(4)
call parse_params("res/input.toml", dict)

nscatt = 0._wp
call init_rng(spread(state%iseed+0, 1, 8), fwd=.true.)


call setup_simulation(packet, array, dict)
! render geometry to voxel format for debugging
call render(array, vector(state%grid%xmax, state%grid%ymax, state%grid%zmax), 200, fname=state%renderfile)
allocate(ds(size(array)))

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
!$omp parallel default(none) shared(array, numproc, start, state, bar, dict)&
!$omp& private(ran, id, ds) reduction(+:nscatt) firstprivate(packet)
    numproc = omp_get_num_threads()
    id = omp_get_thread_num()
#elif MPI
    !nothing
#else
    numproc = 1
    id = 0
#endif
if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

!init photon object

! set seed for rnd generator. id to change seed for each process
call init_rng(spread(123456789+id, 1, 8), fwd=.true.)


bar = pbar(state%nphotons/ 10000)
!$OMP do
!loop over photons 
do j = 1, state%nphotons
    if(mod(j, 10000) == 0)call bar%progress()

    ! Release photon from point source 
    call packet%emit(dict)

    packet%id = id
    ds = 0._wp
    do i = 1, size(ds)
        ds(i) = array(i)%p%evaluate(packet%pos)
    end do
    packet%layer=minloc(ds,dim=1)

    ! Find scattering location
    call tauint2(state%grid, packet, array)
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
        
        !Find next scattering location
        call tauint2(state%grid, packet, array)

    end do
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
    call dict%add_entry("grid_data", 'fluence map')
    call dict%add_entry("real_size", str(state%grid%xmax,7)//" "//str(state%grid%ymax,7)//" "//str(state%grid%zmax,7))
    call dict%add_entry("units", "cm")

    jmeanGLOBAL = normalise_fluence(state%grid, jmeanGLOBAL, state%nphotons)
    call write(jmeanGLOBAL, trim(fileplace)//"jmean/"//state%outfile, dict)
    print*,'write done'
end if

time_taken = get_time() - start
call print_time(time_taken, id)

#ifdef MPI
    call MPI_Finalize()
#endif
end program mcpolar