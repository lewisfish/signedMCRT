module hitStack

    use constants,    only : wp
    use detector_mod, only : hit_t

    implicit none
    
    type :: hit_stack_t
        type(hit_t), allocatable :: data(:)
        integer :: size = 0
        contains
            procedure :: pop => hpop_fn
            procedure :: push => hpush_sub
            procedure :: peek => hpeek_fn
            procedure :: empty => hempty_fn
            procedure :: zero => hzero_sub
    end type hit_stack_t

    integer, parameter :: block_size = 32

    private
    public :: hit_stack_t
contains
    type(hit_t) function hpop_fn(this)
        
        class(hit_stack_t) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            hpop_fn = hit_t(-99._wp)
            return
        end if
        
        hpop_fn = this%data(this%size)
        this%size = this%size - 1

    end function hpop_fn

    subroutine hpush_sub(this, val)

        class(hit_stack_t) :: this
        type(hit_t), intent(in) :: val
        
        type(hit_t), allocatable :: tmp(:)

        if(.not. allocated(this%data) .or. size(this%data) == 0)then
            !allocate space if not yet allocated
            allocate(this%data(block_size))
        elseif(this%size == size(this%data))then
            allocate(tmp(size(this%data) + block_size))
            tmp(1:this%size) = this%data
            call move_alloc(tmp, this%data)
        end if

        this%size = this%size + 1
        this%data(this%size) = val

    end subroutine hpush_sub

    type(hit_t) function hpeek_fn(this)

        class(hit_stack_t) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            hpeek_fn = hit_t(-99._wp)
            return
        end if
        hpeek_fn = this%data(this%size)

    end function hpeek_fn

    logical function hempty_fn(this)

        class(hit_stack_t) :: this

        hempty_fn = (this%size == 0 .or. .not. allocated(this%data))

    end function hempty_fn

    subroutine hzero_sub(this)
                
        class(hit_stack_t) :: this

        this%size = 0

    end subroutine hzero_sub
end module hitStack