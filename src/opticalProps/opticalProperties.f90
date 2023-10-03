module opticalProperties
    !! module implments the optical property abstract type and the types that inheirt from it
    use constants, only : wp
    use piecewiseMod

    implicit none
    !! abstract optical property type
    type, abstract :: opticalProp_base
        !> scattering coeff. \[cm^{-1}\]    
        real(kind=wp) :: mus
        !> absoprtion coeff. \[cm^{-1}\]    
        real(kind=wp) :: mua
        !> g factor
        real(kind=wp) :: hgg
        !> g factor squared
        real(kind=wp) :: g2
        !> refractive index
        real(kind=wp) :: n
        !> \[\kappa = \mu_s + \mu_a\]
        real(kind=wp) :: kappa
        !> \[a = \frac{\mu_s}{\mu_s + \mu_a}\]
        real(kind=wp) :: albedo
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

    interface opticalProp_t
        module procedure opticaProp_new
    end interface

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

    ! optical_property initializer
    type(opticalProp_t) function opticaProp_new(rhs) result(lhs)

        class(opticalProp_base), intent(in) :: rhs
        allocate(lhs%value,source=rhs)

    end function opticaProp_new

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
        
        real(kind=wp) :: wave, tmp

        !setup cdfs
        res%flux = piecewise1D(flux)

        res%mus_a = piecewise1D(mus)
        res%mua_a = piecewise1D(mua)
        res%hgg_a = piecewise1D(hgg)
        res%n_a = piecewise1D(n)

        !sample wavelength so we can sample from other optical properties at the correct points
        call res%flux%sample(wave,  tmp)
        ! sample optical properties
        call res%mus_a%sample(res%mus, wave)
        call res%mua_a%sample(res%mua, wave)
        call res%hgg_a%sample(res%hgg, wave)
        res%g2 = res%hgg**2
        call res%n_a%sample(res%n, wave)

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
        wavelength = 0.0_wp

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