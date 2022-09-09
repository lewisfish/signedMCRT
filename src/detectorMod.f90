module detector_mod
    
    use vector_class
    use constants, only : wp
    
    implicit none
    
    type, public :: hit_t
        type(vector) :: pos, dir
        real(kind=wp) :: value
        integer :: layer
    end type hit_t

    !only needed if using a stack to init with a single null value
    interface hit_t
        module procedure hit_init
    end interface hit_t


    type, abstract :: detector
        ! abstract detector
        ! records non-spatial information, e.g wavelength, time, exit angle etc
        ! currently assumes dectector is in x-y plane
        ! TODO fix this ^
        type(vector)  :: pos
        integer       :: layer, nbins
        real(kind=wp) :: bin_wid
        real(kind=wp), allocatable :: data(:)
        contains
            private
            procedure, public                      :: record_hit => record_hit_sub
            procedure(checkHitInterface), deferred :: check_hit
            procedure                              :: get_bin => get_bin_fn
    end type detector

    abstract interface
        logical function checkHitInterface(this, hitpoint)
            use vector_class
            use constants, only : wp
            import detector, hit_t

            class(detector), intent(inout) :: this
            type(hit_t), intent(IN) :: hitpoint
        end function checkHitInterface
    end interface

    type :: dect_array
        class(detector), pointer :: p => null()
    end type dect_array

    type, extends(detector) :: circle_dect
        real(kind=wp) :: radius
        contains
        procedure :: check_hit  => check_hit_circle
    end type circle_dect

    interface circle_dect
        module procedure init_circle_dect
    end interface circle_dect


    type, extends(detector) :: annulus_dect
        real(kind=wp) :: r1, r2
        contains
        procedure :: check_hit => check_hit_annulus
    end type annulus_dect

    interface annulus_dect
        module procedure init_annulus_dect
    end interface annulus_dect
contains
    
    type(hit_t) function hit_init(val)

        real(kind=wp), intent(in) :: val
        type(vector) :: tmp

        tmp = vector(val, val, val)

        hit_init = hit_t(tmp, tmp, val, int(val))

    end function hit_init

    integer function get_bin_fn(this, value)

        class(detector), intent(IN) :: this
        real(kind=wp),   intent(IN) :: value

        get_bin_fn = nint(value / this%bin_wid) + 1
        get_bin_fn = min(get_bin_fn, this%nbins)

    end function get_bin_fn
    
    subroutine record_hit_sub(this, hitpoint)

        class(detector), intent(INOUT) :: this
        type(hit_t),     intent(IN)    :: hitpoint

        real(kind=wp) :: value
        integer       :: idx

        if(this%check_hit(hitpoint))then
            value = hitpoint%value
            idx = this%get_bin(value)
            !$omp atomic
            this%data(idx) = this%data(idx) + 1
        end if
    end subroutine record_hit_sub
!##############################################################################
!                       CIRCLE DETECTOR
    function init_circle_dect(pos, layer, radius, nbins, maxval) result(out)

        type(vector),  intent(IN) :: pos
        integer,       intent(IN) :: layer, nbins
        real(kind=wp), intent(IN) :: radius, maxval

        type(circle_dect) :: out

        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%radius = radius
        allocate(out%data(out%nbins))
        out%bin_wid = maxval / real(nbins, kind=wp)

    end function init_circle_dect

    logical function check_hit_circle(this, hitpoint)

        class(circle_dect), intent(INOUT) :: this
        type(hit_t),        intent(IN)    :: hitpoint

        real(kind=wp) :: newpos

        check_hit_circle = .false.
        newpos = sqrt((hitpoint%pos%x - this%pos%x)**2 + (hitpoint%pos%y - this%pos%y)**2 + (hitpoint%pos%z - this%pos%z)**2)
        if(newpos <= this%radius)then
            check_hit_circle = .true.
        end if
    end function check_hit_circle
! ##########################################################################
!                       ANNULUS DETECTOR
    function init_annulus_dect(pos, layer, r1, r2, nbins, maxval) result(out)

        type(vector),  intent(IN) :: pos
        integer,       intent(IN) :: layer, nbins
        real(kind=wp), intent(IN) :: r1, r2, maxval

        type(annulus_dect) :: out

        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%r1 = r1
        out%r2 = r2
        allocate(out%data(out%nbins))
        out%bin_wid = maxval / real(nbins, kind=wp)

    end function init_annulus_dect

    logical function check_hit_annulus(this, hitpoint)

        class(annulus_dect), intent(INOUT) :: this
        type(hit_t),         intent(IN)    :: hitpoint

        real(kind=wp) :: newpos

        check_hit_annulus = .false.
        newpos = sqrt((hitpoint%pos%x - this%pos%x)**2 + (hitpoint%pos%y - this%pos%y)**2 + (hitpoint%pos%z - this%pos%z)**2)
        if(newpos >= this%r1 .and. newpos <= this%r2)then
            check_hit_annulus = .true.
        end if

    end function check_hit_annulus
end module detector_mod
! program test
!     use detector_mod
!     use vector_class
!     use constants, only : wp
!     implicit none

!     type(hit_t) :: hit
!     type(vector) :: pos, dir
!     integer :: layer

!     type(circle_dect) :: dect_c
!     type(annulus_dect) :: dect_a

!     dect_c = circle_dect(vector(0._wp, 0._wp, 0._wp), 1, .5_wp, 100, 100._wp)
!     dect_a = annulus_dect(vector(0._wp, 0._wp, 0._wp), 1, .25_wp, .5_wp, 100, 100._wp)

!     layer = 1
!     pos = vector(0._wp, .5_wp, 0._wp)
!     dir = vector(0._wp, 0._wp, 1._wp)
    
!     hit = hit_t(pos, dir, 99._wp, layer)
!     call dect_c%record_hit(hit)
!     print*,sum(dect_c%data)

!     pos = vector(0._wp, .25_wp, 0._wp)
!     dir = vector(0._wp, 0._wp, 1._wp)
    
!     hit = hit_t(pos, dir, 99._wp, layer)
!     call dect_a%record_hit(hit)
!     print*,sum(dect_a%data)
! end program test