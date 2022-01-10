module stackMod

    use constants, only : wp

    implicit none

    type :: istack
        integer, allocatable :: data(:) ! istack
        integer              :: size = 0
        contains
            procedure :: pop   => ipop_fn
            procedure :: push  => ipush_sub
            procedure :: peek  => ipeek_fn
            procedure :: empty => iempty_fn
            procedure :: zero  => izero_sub
            procedure :: write => iwrite_sub
            procedure :: write_empty => iwrite_empty_sub
    end type istack


    type :: rstack
        real(kind=wp), allocatable :: data(:) ! istack
        integer              :: size = 0
        contains
            procedure :: pop   => rpop_fn
            procedure :: push  => rpush_sub
            procedure :: peek  => rpeek_fn
            procedure :: empty => rempty_fn
            procedure :: zero  => rzero_sub
            procedure :: write => rwrite_sub
            procedure :: write_empty => rwrite_empty_sub
    end type rstack

    integer, parameter :: block_size = 3

    contains

    subroutine iwrite_empty_sub(this, u)

        implicit none
        
        class(istack) :: this
        integer, intent(IN) :: u

        call this%zero()
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine iwrite_empty_sub

    subroutine iwrite_sub(this, u)

        implicit none
        
        class(istack) :: this
        integer, intent(IN) :: u

        do while(.not. this%empty())
            write(u,"(F10.7,1x)")this%pop()
        end do
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine iwrite_sub

    subroutine izero_sub(this)
        
        implicit none
        
        class(istack) :: this

        integer :: tmp

        do while(.not. this%empty())
            tmp = this%pop()
        end do

    end subroutine izero_sub

    integer function ipop_fn(this)
    ! pop top enrty off istack
        implicit none

        class(istack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            !if nothing in istack send back garbage data
            ipop_fn = -99.
            return
        end if
        ipop_fn = this%data(this%size)
        this%size = this%size - 1

    end function ipop_fn


    integer function ipeek_fn(this)

        implicit none

        class(istack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            ipeek_fn = -99.
            return
        end if
        ipeek_fn = this%data(this%size)

    end function ipeek_fn

    logical function iempty_fn(this)

        implicit none

        class(istack) :: this

        iempty_fn = (this%size == 0 .or. .not. allocated(this%data))

    end function iempty_fn

    subroutine ipush_sub(this, pt)
    ! add pt to istack
        implicit none

        class(istack) :: this

        integer, intent(IN)  :: pt
        integer, allocatable :: tmp(:)

        if(.not. allocated(this%data))then
            ! Allocate space if not yet done
            allocate(this%data(block_size))
        elseif(this%size == size(this%data))then
            ! Grow the allocated space
            allocate(tmp(size(this%data)+block_size))
            tmp(1:this%size) = this%data
            call move_alloc(tmp,this%data)
        end if

        ! Store the data in the istack
        this%size = this%size + 1
        this%data(this%size) = pt
    end subroutine ipush_sub



! *****************************************************
    subroutine rwrite_empty_sub(this, u)

        implicit none
        
        class(rstack) :: this
        integer, intent(IN) :: u

        call this%zero()
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine rwrite_empty_sub

    subroutine rwrite_sub(this, u)

        implicit none
        
        class(rstack) :: this
        integer, intent(IN) :: u

        do while(.not. this%empty())
            write(u,"(F10.7,1x)")this%pop()
        end do
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine rwrite_sub

    subroutine rzero_sub(this)
        
        implicit none
        
        class(rstack) :: this

        real(kind=wp) :: tmp

        do while(.not. this%empty())
            tmp = this%pop()
        end do

    end subroutine rzero_sub

    function rpop_fn(this) result(res)
    ! pop top enrty off rstack
        implicit none

        class(rstack) :: this
        real(kind=wp) :: res

        if(this%size == 0 .or. .not. allocated(this%data))then
            !if nothing in rstack send back garbage data
            res = -99._wp
            return
        end if
        res = this%data(this%size)
        this%size = this%size - 1

    end function rpop_fn


    function rpeek_fn(this) result(res)

        implicit none

        class(rstack) :: this
        real(kind=wp) :: res

        if(this%size == 0 .or. .not. allocated(this%data))then
            res = -99._wp
            return
        end if
        res = this%data(this%size)

    end function rpeek_fn

    logical function rempty_fn(this)

        implicit none

        class(rstack) :: this

        rempty_fn = (this%size == 0 .or. .not. allocated(this%data))

    end function rempty_fn

    subroutine rpush_sub(this, pt)
    ! add pt to rstack
        implicit none

        class(rstack) :: this

        real(kind=wp), intent(IN)  :: pt
        real(kind=wp), allocatable :: tmp(:)

        if(.not. allocated(this%data))then
            ! Allocate space if not yet done
            allocate(this%data(block_size))
        elseif(this%size == size(this%data))then
            ! Grow the allocated space
            allocate(tmp(size(this%data)+block_size))
            tmp(1:this%size) = this%data
            call move_alloc(tmp,this%data)
        end if

        ! Store the data in the rstack
        this%size = this%size + 1
        this%data(this%size) = pt
    end subroutine rpush_sub
end module stackMod