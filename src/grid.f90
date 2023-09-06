module gridMod
    
    use constants, only : wp

    implicit none

    !! Grid class
    type :: cart_grid
        !> number of voxels in each cardinal direction for fluence grid
        integer       :: nxg, nyg, nzg
        !> half size of each dimension in fluence grid. 
        real(kind=wp) :: xmax, ymax, zmax
        !> Delta is the round off for near voxel cell walls
        real(kind=wp) :: delta
        !> position of each cell wall in fluence grid
        real(kind=wp), allocatable :: xface(:), yface(:), zface(:)
        contains
        procedure :: get_voxel
    end type cart_grid

    interface cart_grid
        module procedure init_grid
    end interface cart_grid

    public  :: cart_grid, init_grid
    private

    contains

    function get_voxel(this, pos) result(res)
        !! get current voxel the photon packet is in
        use vector_class
        
        !> grid class
        class(cart_grid)         :: this
        !> current vector position of photon packet
        type(vector), intent(IN) :: pos
    
        integer :: res(3)

        res(1) = int(this%nxg*(pos%x+this%xmax)/(2._wp*this%xmax))+1
        res(2) = int(this%nyg*(pos%y+this%ymax)/(2._wp*this%ymax))+1
        res(3) = int(this%nzg*(pos%z+this%zmax)/(2._wp*this%zmax))+1

    end function get_voxel

    type(cart_grid) function init_grid(nxg, nyg, nzg, xmax, ymax, zmax)
    !! setup grid
        !> number of voxels in each cardinal direction for fluence grid
        integer,       intent(IN) :: nxg, nyg, nzg
        !> half size of each dimension in fluence grid. 
        real(kind=wp), intent(IN) :: xmax, ymax, zmax
        
        integer :: i

        init_grid%nxg = nxg
        init_grid%nyg = nyg
        init_grid%nzg = nzg

        init_grid%xmax = xmax
        init_grid%ymax = ymax
        init_grid%zmax = zmax

        allocate(init_grid%xface(nxg + 1), init_grid%yface(nyg + 1), init_grid%zface(nzg + 2))

        init_grid%xface = 0._wp
        init_grid%yface = 0._wp
        init_grid%zface = 0._wp

        ! Set small distance for use in optical depth integration routines 
        ! for roundoff effects when crossing cell walls
        init_grid%delta = 1.e-8_wp * min(((2._wp*xmax)/nxg), ((2._wp*ymax)/nyg), ((2._wp*zmax)/nzg))


        do i = 1, nxg + 1
            init_grid%xface(i) = (i - 1) * 2._wp * xmax/nxg
        end do

        do i = 1, nyg + 1
            init_grid%yface(i) = (i - 1) * 2._wp * ymax/nyg
        end do

        do i = 1, nzg + 2
            init_grid%zface(i) = (i - 1) * 2._wp * zmax/nzg
        end do

    end function init_grid
end module gridMod