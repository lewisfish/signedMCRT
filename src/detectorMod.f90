module detector_mod

    use vector_class
    use constants, only : wp
    
    implicit none
    
    type, public :: hit_t
        type(vector)  :: pos, dir
        real(kind=wp) :: value
        integer       :: layer
    end type hit_t

    !only needed if using a stack to init with a single null value
    interface hit_t
        module procedure hit_init
    end interface hit_t

    type, abstract :: detector
        ! abstract detector
        type(vector)  :: pos, dir
        integer :: layer
        logical :: trackHistory
        contains
            private
            procedure(recordHitInterface), deferred, public :: record_hit
            procedure(checkHitInterface),  deferred :: check_hit
    end type detector

    abstract interface
        logical function checkHitInterface(this, hitpoint)
            use vector_class
            use constants, only : wp
            import detector, hit_t

            class(detector), intent(inout) :: this
            type(hit_t),     intent(in)    :: hitpoint
        end function checkHitInterface

        subroutine recordHitInterface(this, hitpoint, history)
            use constants,     only : wp
            use historyStack,  only : history_stack_t
            use vector_class
            import detector, hit_t

            class(detector),       intent(inout) :: this
            type(hit_t),           intent(in)    :: hitpoint
            type(history_stack_t), intent(inout) :: history
        end subroutine recordHitInterface
    end interface

    type, abstract, extends(detector) :: detector1D
        integer       :: nbins
        real(kind=wp) :: bin_wid
        real(kind=wp), allocatable :: data(:)
        contains
        procedure :: record_hit => record_hit_1D_sub
    end type detector1D

    type, abstract, extends(detector) :: detector2D
        integer       :: nbinsX, nbinsY
        real(kind=wp) :: bin_wid_x, bin_wid_y
        real(kind=wp), allocatable :: data(:,:)
        contains
        procedure :: record_hit => record_hit_2D_sub
    end type detector2D

    type :: dect_array
        class(detector), pointer :: p => null()
    end type dect_array

    type, extends(detector1D) :: circle_dect
        real(kind=wp) :: radius
        contains
        procedure :: check_hit  => check_hit_circle
    end type circle_dect

    interface circle_dect
        module procedure init_circle_dect
    end interface circle_dect


    type, extends(detector1D) :: annulus_dect
        real(kind=wp) :: r1, r2
        contains
        procedure :: check_hit => check_hit_annulus
    end type annulus_dect

    interface annulus_dect
        module procedure init_annulus_dect
    end interface annulus_dect


    type, extends(detector2D) :: camera
        type(vector)  :: n, p2, p3, e1, e2
        real(kind=wp) :: width, height
        contains
        procedure :: check_hit => check_hit_camera
    end type camera

    interface camera
        module procedure init_camera
    end interface camera

    private
    public :: detector, camera, annulus_dect, circle_dect, dect_array

contains
    type(hit_t) function hit_init(val)

        real(kind=wp), intent(in) :: val
        type(vector) :: tmp

        tmp = vector(val, val, val)

        hit_init = hit_t(tmp, tmp, val, int(val))

    end function hit_init
   
    subroutine record_hit_1D_sub(this, hitpoint, history)

        use historyStack, only : history_stack_t

        class(detector1D),     intent(inout) :: this
        type(hit_t),           intent(in)    :: hitpoint
        type(history_stack_t), intent(inout) :: history

        real(kind=wp) :: value
        integer       :: idx

        if(this%check_hit(hitpoint))then
            value = hitpoint%value
            idx = min(nint(value / this%bin_wid) + 1, this%nbins)
            !$omp atomic
            this%data(idx) = this%data(idx) + 1
            if(this%trackHistory)then
                call history%write()
            end if
        end if
        call history%zero()
    end subroutine record_hit_1D_sub

    subroutine record_hit_2D_sub(this, hitpoint, history)

        use historyStack, only : history_stack_t

        class(detector2D),     intent(inout) :: this
        type(hit_t),           intent(in)    :: hitpoint
        type(history_stack_t), intent(inout) :: history

        real(kind=wp), volatile :: x, y
        integer       :: idx, idy

        if(this%check_hit(hitpoint))then
            x = hitpoint%pos%z + this%pos%x
            y = hitpoint%pos%y + this%pos%y
            idx = min(int(x / this%bin_wid_x) + 1, this%nbinsX)
            idy = min(int(y / this%bin_wid_y) + 1, this%nbinsY)
            if(idx < 1)idx = this%nbinsX
            if(idy < 1)idy = this%nbinsY
            !$omp atomic
            this%data(idx, idy) = this%data(idx, idy) + 1
            if(this%trackHistory)then
                call history%write()
            end if
        end if
        call history%zero()
        end subroutine record_hit_2D_sub
!##############################################################################
!                       CIRCLE DETECTOR
    function init_circle_dect(pos, dir, layer, radius, nbins, maxval, trackHistory) result(out)

        type(vector),  intent(in) :: pos, dir
        integer,       intent(in) :: layer, nbins
        real(kind=wp), intent(in) :: radius, maxval
        logical,       intent(in) :: trackHistory

        type(circle_dect) :: out

        out%dir = dir
        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%radius = radius
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = maxval / real(nbins-1, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_circle_dect

    logical function check_hit_circle(this, hitpoint)
        
        use geometry, only : intersectCircle

        class(circle_dect), intent(INOUT) :: this
        type(hit_t),        intent(IN)    :: hitpoint
        
        real(kind=wp) :: t 

        check_hit_circle = .false.
        if(this%layer /= hitpoint%layer)return
        check_hit_circle = intersectCircle(this%dir, this%pos, this%radius, hitpoint%pos, hitpoint%dir, t)
        if(check_hit_circle)then
            if(t > 5e-3_wp)check_hit_circle=.false.
        end if
    end function check_hit_circle
! ##########################################################################
!                       ANNULUS DETECTOR
    function init_annulus_dect(pos, layer, r1, r2, nbins, maxval, trackHistory) result(out)

        type(vector),  intent(IN) :: pos
        integer,       intent(IN) :: layer, nbins
        real(kind=wp), intent(IN) :: r1, r2, maxval
        logical,       intent(in) :: trackHistory

        type(annulus_dect) :: out

        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%r1 = r1
        out%r2 = r2
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = maxval / real(nbins, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_annulus_dect

    logical function check_hit_annulus(this, hitpoint)

        class(annulus_dect), intent(INOUT) :: this
        type(hit_t),         intent(IN)    :: hitpoint

        real(kind=wp) :: newpos

        check_hit_annulus = .false.
        if(this%layer /= hitpoint%layer)return
        newpos = sqrt((hitpoint%pos%x - this%pos%x)**2 + (hitpoint%pos%y - this%pos%y)**2 + (hitpoint%pos%z - this%pos%z)**2)
        if(newpos >= this%r1 .and. newpos <= this%r2)then
            check_hit_annulus = .true.
        end if

    end function check_hit_annulus
! ##########################################################################
!                       CAMERA
    function init_camera(p1, p2, p3, layer, nbins, maxval, trackHistory) result(out)

        type(vector),  intent(in) :: p1, p2, p3
        integer,       intent(in) :: layer, nbins
        real(kind=wp), intent(in) :: maxval
        logical,       intent(in) :: trackHistory
        type(camera) :: out

        out%pos = p1
        out%p2 = p2
        out%p3 = p3
        out%e1 = p2 - p1
        out%e2 = p3 - p1
        out%width = length(out%e1)
        out%height = length(out%e2)
        out%n = out%e2 .cross. out%e1
        out%n = out%n%magnitude()
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbinsX = nbins + 1
        out%nbinsY = nbins + 1
        allocate(out%data(out%nbinsX, out%nbinsY))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid_x = 1._wp
            out%bin_wid_y = 1._wp
        else
            out%bin_wid_x = maxval / real(out%nbinsX, kind=wp)
            out%bin_wid_y = maxval / real(out%nbinsY, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_camera

    logical function check_hit_camera(this, hitpoint)
    ! https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
        class(camera), intent(inout) :: this
        type(hit_t),   intent(in)    :: hitpoint

        real(kind=wp) :: t, proj1, proj2
        type(vector)  :: v

        check_hit_camera = .false.
        if(this%layer /= hitpoint%layer)return

        t = ((this%pos - hitpoint%pos) .dot. this%n) / (hitpoint%dir .dot. this%n)
        if(t >= 0._wp)then
            v = (hitpoint%pos + t * hitpoint%dir) - this%pos
            proj1 = (v .dot. this%e1) / this%width
            proj2 = (v .dot. this%e2) / this%height
            if((proj1 < this%width .and. proj1 > 0._wp) .and. (proj2 < this%height .and. proj2 > 0._wp))then
                check_hit_camera = .true.
            end if
        end if
    end function check_hit_camera
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