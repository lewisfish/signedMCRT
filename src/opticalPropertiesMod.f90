module opticalProperties

    use constants, only : wp
    use piecewiseMod

    implicit none

    type, abstract :: opticalProp_base
        real(kind=wp) :: mus, mua, hgg, g2, n, kappa, albedo
        contains
            procedure(updateInterface), deferred :: update
    end type opticalProp_base

    abstract interface
        subroutine updateInterface(this, wavelength)
            use constants, only : wp
            use piecewiseMod
            import opticalProp_base
            implicit none
            class(opticalProp_base), intent(inout) :: this
            real(kind=wp),           intent(out)   :: wavelength
        end subroutine updateInterface
    end interface

    type, extends(opticalProp_base) :: opticalProp_t
        class(opticalProp_base), allocatable :: value
        contains
            procedure :: update => update_opticalProp_t
            procedure, private :: opticalProp_t_assign
            generic :: assignment(=) => opticalProp_t_assign
    end type opticalProp_t

    type, extends(opticalProp_base) :: mono
        contains
            procedure :: update => updateMono
    end type mono

    type, extends(opticalProp_base) :: spectral
        type(piecewise1D), private :: mus_a, mua_a, hgg_a, n_a, flux
        contains
            procedure :: update => updateSpectral
    end type spectral

    interface spectral
        module procedure init_spectral
    end interface spectral

    interface mono
        module procedure init_mono
    end interface mono

    private
    public :: spectral, mono, opticalProp_base, opticalProp_t

    contains

    subroutine opticalProp_t_assign(lhs, rhs)

        class(opticalProp_t),    intent(inout) :: lhs
        class(opticalProp_base), intent(in)    :: rhs

        if (allocated(lhs%value))deallocate(lhs%value)
        ! Prevent nested derived type
        select type (rhsT=>rhs)
            class is (opticalProp_t)
                if(allocated(rhsT%value))allocate(lhs%value,source=rhsT%value)
            class default
                allocate(lhs%value,source=rhsT)
        end select

    end subroutine opticalProp_t_assign

    subroutine update_opticalProp_t(this, wavelength)

        class(opticalProp_t), intent(inout) :: this
        real(kind=wp),        intent(out)   :: wavelength

        call this%value%update(wavelength)

    end subroutine update_opticalProp_t

    type(mono) function init_mono(mus, mua, hgg, n) result(res)

        real(kind=wp), intent(in) :: mus, mua, hgg, n

        res%mus = mus
        res%mua = mua

        res%kappa = mus + mua
        if(res%mua < 1e-9_wp)then          
            res%albedo = 1.
        else
            res%albedo = res%mus / res%kappa
        end if

        res%hgg = hgg
        res%g2 = hgg**2
        res%n = n

    end function init_mono

    type(spectral) function init_spectral(mus, mua, hgg, n, flux) result(res)

        real(kind=wp), allocatable, intent(in) :: mus(:, :), mua(:, :), hgg(:, :), n(:, :), flux(:, :)
        
        real(kind=wp) :: tmp

        res%flux = piecewise1D(flux(:, :))
        res%mus_a = piecewise1D(mus(:, :))
        res%mua_a = piecewise1D(mua(:, :))
        res%hgg_a = piecewise1D(hgg(:, :))
        res%n_a = piecewise1D(n(:, :))

        call res%mus_a%sample(res%mus, tmp)
        call res%mua_a%sample(res%mua, tmp)
        call res%hgg_a%sample(res%hgg, tmp)
        res%g2 = res%hgg**2
        call res%n_a%sample(res%n, tmp)

        res%kappa =  res%mus + res%mua
        if(res%mua < 1e-9_wp)then          
            res%albedo = 1.
        else
            res%albedo = res%mus / res%kappa
        end if
    end function init_spectral

    subroutine updateMono(this, wavelength)

        implicit none

        class(Mono),   intent(inout) :: this
        real(kind=wp), intent(out)   :: wavelength

        ! don't do anything as wavelength will not change

    end subroutine updateMono


    subroutine updateSpectral(this, wavelength)

        implicit none

        class(spectral), intent(inout) :: this
        real(kind=wp),   intent(out)   :: wavelength

        real(kind=wp) :: tmp

        !get wavelength
        call this%flux%sample(wavelength, tmp)

        !update mus
        call this%mus_a%sample(this%mus, tmp, wavelength)

        !update mua
        call this%mua_a%sample(this%mua, tmp, wavelength)

        !update hgg
        call this%hgg_a%sample(this%hgg, tmp, wavelength)
        this%g2 = this%hgg**2

        !update n
        call this%n_a%sample(this%n, tmp, wavelength)

        !update kappa and albedo
        this%kappa =  this%mus + this%mua
        this%albedo = this%mus / this%kappa

    end subroutine updateSpectral
end module opticalProperties
! program test

!     use opticalProperties
!     use stdlib_io, only: loadtxt

!     implicit none

!     real(kind=wp), allocatable :: solar(:,:), water(:, :), hgg(:, :), n(:, :), mus(:, :)
!     type(spectral) :: optProp
!     integer :: i

!     call loadtxt("res/solar.dat", solar)
!     call loadtxt("res/water.dat", water)

!     ! allocate(solar(10, 2))
!     ! solar(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900]
!     ! solar(:, 2) = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

!     allocate(n(10, 2))

!     n(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900]
!     n(:, 2) = [1.0, 1.5, 1.5, 1.0, 1.5, 1.8, 1.9, 2.0, 2.1]

!     allocate(mus(10, 2))
!     mus(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900]
!     mus(:, 2) = [0.0, 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.5]

!     allocate(hgg(3, 2))
!     hgg(:, 1) = [100, 450, 900]
!     hgg(:, 2) = [0.9, 0.9, 0.9]

!     optProp = spectral(mus, water, hgg, n, solar)

!     do i = 1, 500
!         call optProp%update()
!     end do

! end program test