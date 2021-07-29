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
use random, only : ran2, init_rng, ranu

!subroutines
use subs
use gridmod
use inttau2
use stokes_mod
use writer_mod
use vector_class
use sdfs
use utils, only : pbar

implicit none

type(photon)     :: packet
type(cart_grid)  :: grid
integer          :: nphotons, iseed, j
double precision :: nscatt
real             :: ran, start, time_taken
type(pbar)       :: bar

type(cylinder) :: array(3)
! mpi/mp variables
integer :: id, numproc
real    :: nscattGLOBAL

call setup_simulation(nphotons, grid)

array(1) = cylinder(1000., 0.5)
array(2) = cylinder(1000., 0.6)
array(3) = cylinder(1000., 1.75)

start = get_time()
id = 0
!init MPI
#ifdef MPI
    call MPI_init()
    call MPI_Comm_size(MPI_COMM_WORLD, numproc)
    call MPI_Comm_rank(MPI_COMM_WORLD, id)
#endif


if(id == 0)then
   print*,'# of photons to run',nphotons
end if

!loop over photons 
#ifdef _OPENMP
!$omp parallel default(none) shared(array, grid, numproc, start, nphotons, bar)&
!$omp& private(ran, id, packet, iseed, skip) reduction(+:nscatt)
    numproc = omp_get_num_threads()
    id = omp_get_thread_num()
#elif MPI
    !nothing
#else
    numproc = 1
    id = 0
#endif
if(id == 0)print*,'Photons now running on', numproc,' cores.'

nscatt = 0
!init photon object
packet = photon("uniform")

! set seed for rnd generator. id to change seed for each process
call init_rng(123456789)


bar = pbar(nphotons/ 100000)

!$OMP do
do j = 1, nphotons

    if(mod(j, 100000) == 0)call bar%progress()

    ! Release photon from point source 
    call packet%emit(grid, iseed)

    ! Find scattering location
    call tauint2(packet, grid, array)

    ! Photon scatters in grid until it exits (tflag=TRUE) 
    do while(packet%tflag .eqv. .FALSE.)
        ran = ran2()
        if(ran < grid%albedo(packet%layer))then!interacts with tissue
            call stokes(packet, grid)
            nscatt = nscatt + 1        
        else
            packet%tflag = .true.
            exit
        end if
        
        !Find next scattering location
        call tauint2(packet, grid, array)
    
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
    print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons)
#else
    print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
#endif
    !write out files
    call writer()
    print*,'write done'
end if

time_taken = get_time() - start
call print_time(time_taken, id)

#ifdef MPI
    call MPI_Finalize()
#endif
end program mcpolar