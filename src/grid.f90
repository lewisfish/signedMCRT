module gridMod
    
    implicit none

    type :: cart_grid
        integer :: nxg, nyg, nzg ! number of voxels in each cardinal direction for fluence grid
        real    :: xmax, ymax, zmax, delta ! size of each dimension in fluence grid. Delta is the round of for near voxel cell walls
        real, allocatable :: xface(:), yface(:), zface(:) ! position of each cell wall in fluence grid
    end type cart_grid

    interface cart_grid
        module procedure init_grid
    end interface cart_grid

    public  :: cart_grid, init_grid
    private

    contains

    type(cart_grid) function init_grid(nxg, nyg, nzg, xmax, ymax, zmax)
    ! setup grid
    !
    !
        implicit none
        
        integer, intent(IN) :: nxg, nyg, nzg 
        real,    intent(IN) :: xmax, ymax, zmax
        
        integer :: i

        init_grid%nxg = nxg
        init_grid%nyg = nyg
        init_grid%nzg = nzg

        init_grid%xmax = xmax
        init_grid%ymax = ymax
        init_grid%zmax = zmax

        allocate(init_grid%xface(nxg + 1), init_grid%yface(nyg + 1), init_grid%zface(nzg + 2))

        init_grid%xface = 0.
        init_grid%yface = 0.
        init_grid%zface = 0.

        ! Set small distance for use in optical depth integration routines 
        ! for roundoff effects when crossing cell walls
        init_grid%delta = 1.e-8 * min(((2.*xmax)/nxg), ((2.*ymax)/nyg), ((2.*zmax)/nzg))


        do i = 1, nxg + 1
            init_grid%xface(i) = (i - 1) * 2. * xmax/nxg
        end do

        do i = 1, nyg + 1
            init_grid%yface(i) = (i - 1) * 2. * ymax/nyg
        end do

        do i = 1, nzg + 2
            init_grid%zface(i) = (i - 1) * 2. * zmax/nzg
        end do

    end function init_grid
end module gridMod