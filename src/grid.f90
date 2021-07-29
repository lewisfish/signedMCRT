module gridMod
    
    implicit none

    type :: cart_grid
        integer :: nxg, nyg, nzg
        real    :: xmax, ymax, zmax, kappa(3), albedo(3), hgg(3), g2(3), delta
        real, allocatable :: xface(:), yface(:), zface(:)
    end type cart_grid

    interface cart_grid
        module procedure init_grid
    end interface cart_grid

    public  :: cart_grid, init_grid
    private

    contains

    type(cart_grid) function init_grid(nxg, nyg, nzg, xmax, ymax, zmax, n1, n2)

        use ch_opt

        implicit none
        
        integer, intent(IN) :: nxg, nyg, nzg 
        real,    intent(IN) :: xmax, ymax, zmax, n1, n2
        
        real    :: kappa(3), albedo(3), hgg(3), g2(3)
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

        call init_opt1(kappa, albedo, hgg, g2)

        init_grid%hgg = hgg
        init_grid%albedo = albedo
        init_grid%g2 = g2
        init_grid%kappa = kappa

    end function init_grid
end module gridMod