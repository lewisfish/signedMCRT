module detector_mod
    !! Module contains photon detector abstract class and the derived types the inherit from it
    !! not fully implmented
    use vector_class
    use constants, only : wp
    
    implicit none
    !> Hit type, which records possible interaction information
    type :: hit_t
        !> Poition of the interaction
        type(vector)  :: pos
        !> Direction the photon came from
        type(vector)  :: dir
        !> Value to deposit
        real(kind=wp) :: value
        !> Layer ID of interaction
        integer       :: layer
    end type hit_t

    !only needed if using a stack to init with a single null value
    interface hit_t
        module procedure hit_init
    end interface hit_t

    !> abstract detector
    type, abstract :: detector
        !> position of the detector
        type(vector)  :: pos
        !> Surface normal of the detector
        type(vector)  :: dir
        !> Layer ID of the detector
        integer :: layer
        !> Boolean, if true store the history of the photon prior to detection.
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
    
    !> 1D detector type. Records linear information
    type, abstract, extends(detector) :: detector1D
        !> Number of bins
        integer       :: nbins
        !> Bin width
        real(kind=wp) :: bin_wid
        !> Bins
        real(kind=wp), allocatable :: data(:)
        contains
        procedure :: record_hit => record_hit_1D_sub
    end type detector1D
    
    !> 2D detecctor type. Records spatial information
    type, abstract, extends(detector) :: detector2D
        !> Number of bins in x dimension (detector space)
        integer       :: nbinsX
        !> Number of bins in y dimension (detector space)
        integer       :: nbinsY
        !> Bin width in the x dimension
        real(kind=wp) :: bin_wid_x
        !> Bin width in the y dimension
        real(kind=wp) :: bin_wid_y
        !> Bins
        real(kind=wp), allocatable :: data(:,:)
        contains
        procedure :: record_hit => record_hit_2D_sub
    end type detector2D

    private
    public :: detector, detector1D, detector2D, hit_t

contains
    type(hit_t) function hit_init(val)

        real(kind=wp), intent(in) :: val
        type(vector) :: tmp

        tmp = vector(val, val, val)

        hit_init = hit_t(tmp, tmp, val, int(val))

    end function hit_init
   
    subroutine record_hit_1D_sub(this, hitpoint, history)
        !! check if a hit is on the detector and record it if so

        use historyStack,  only : history_stack_t
        use sim_state_mod, only : state

        class(detector1D),     intent(inout) :: this
        !> Interaction information
        type(hit_t),           intent(in)    :: hitpoint
        !> Photon packet history
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
        if(state%trackHistory)call history%zero()
    end subroutine record_hit_1D_sub

    subroutine record_hit_2D_sub(this, hitpoint, history)
        !! check if a hit is on the detector and record it if so

        use historyStack, only : history_stack_t
        use sim_state_mod, only : state

        class(detector2D),     intent(inout) :: this
        !> Interaction information
        type(hit_t),           intent(in)    :: hitpoint
        !> Photon packet history
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
        if(state%trackHistory)call history%zero()
        end subroutine record_hit_2D_sub
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