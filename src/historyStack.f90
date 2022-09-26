module historyStack

    use constants,    only : wp
    use vector_class, only : vector

    implicit none

    type :: history_stack_t
        type(vector), allocatable :: data(:)
        integer :: size, vertex_counter, edge_counter
        contains
            procedure :: pop   => histpop_fn
            procedure :: push  => histpush_sub
            procedure :: peek  => histpeek_fn
            procedure :: empty => histempty_fn
            procedure :: zero  => histzero_sub
            procedure :: write => histwrite_sub
            procedure :: finish => histfinish_sub
    end type history_stack_t

    interface history_stack_t
        module procedure init_historyStack
    end interface

    integer, parameter :: block_size = 32

    private
    public ::  history_stack_t

contains

    type(history_stack_t) function init_historyStack()

        init_historyStack%size = 0
        init_historyStack%vertex_counter = 0
        init_historyStack%edge_counter = 0

    end function init_historyStack

    type(vector) function histpop_fn(this)
        
        class(history_stack_t) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            histpop_fn = vector(-99._wp, -99._wp, -99._wp)
            return
        end if
        
        histpop_fn = this%data(this%size)
        this%size = this%size - 1

    end function histpop_fn

    subroutine histpush_sub(this, val)

        class(history_stack_t) :: this
        type(vector), intent(in) :: val
        
        type(vector), allocatable :: tmp(:)

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

    end subroutine histpush_sub

    type(vector) function histpeek_fn(this)

        class(history_stack_t) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            histpeek_fn = vector(-99._wp, -99._wp, -99._wp)
            return
        end if
        histpeek_fn = this%data(this%size)

    end function histpeek_fn

    logical function histempty_fn(this)

        class(history_stack_t) :: this

        histempty_fn = (this%size == 0 .or. .not. allocated(this%data))

    end function histempty_fn

    subroutine histzero_sub(this)
                
        class(history_stack_t) :: this

        this%size = 0

    end subroutine histzero_sub

    subroutine histwrite_sub(this)

        use constants, only : fileplace
        use omp_lib
        use string_utils, only : str

        class(history_stack_t) :: this

        type(vector) :: v
        integer :: u, io, id, counter, i
        logical :: res
        
        id = 0!omp_get_thread_num()
        if(id == 0)then
            inquire(file=trim(fileplace)//"photPos.ply", exist=res)
            if(res)then
                open(newunit=u,file=trim(fileplace)//"photPos.ply", status="old", position="append")
            else
                open(newunit=u,file=trim(fileplace)//"photPos.ply", status="new")
                write(u,"(a)") "ply"//new_line("a")//"format ascii 1.0"//new_line("a")//"element vertex "//str(this%size)
                write(u,"(a)") "property float x"
                write(u,"(a)") "property float y"
                write(u,"(a)") "property float z"
                write(u,"(a)") "element edge"
                write(u,"(a)") "property int vertex1"
                write(u,"(a)") "property int vertex2"
                write(u,"(a)") "end_header"
            end if
            inquire(file=trim(fileplace)//"photPos.ply2", exist=res)
            if(res)then
                open(newunit=io,file=trim(fileplace)//"photPos.ply2", status="old", position="append")
            else
                open(newunit=io,file=trim(fileplace)//"photPos.ply2", status="new")
            end if

            counter = this%vertex_counter
            do i = 1, this%size-1
                write(io, "(2(i0,1x))") counter, counter+1
                counter = counter + 1
                this%edge_counter = this%edge_counter + 1
            end do
            close(io)
            do while(.not. this%empty())
                v = this%pop()
                write(u, "(3(es15.8e2,1x))") v
                this%vertex_counter = this%vertex_counter + 1
            end do
            close(u)

        end if

    end subroutine histwrite_sub

    subroutine histfinish_sub(this)

        use constants,    only : fileplace
        use string_utils, only : str

        class(history_stack_t) :: this

        ! integer :: u
        ! this is the easiest way to edit the vertex count as we don't know how many photons we will track when writing the header.
        ! this saves stroing all photons data in RAM for duration of simulation.
        ! taken from: https://stackoverflow.com/a/11145362
        call execute_command_line("sed -i '3s#.*#element vertex "//str(this%vertex_counter)//"#' "//trim(fileplace)//"photPos.ply")
        call execute_command_line("sed -i '7s#.*#element edge "//str(this%edge_counter)//"#' "//trim(fileplace)//"photPos.ply")
       
        call execute_command_line("cat data/photPos.ply2 >> data/photPos.ply")
        print*,"joined!",this%vertex_counter,this%edge_counter
    end subroutine histfinish_sub
end module historyStack



!json
! inquire(file=trim(fileplace)//"photPos.json", exist=res)
! if(res)then
!     open(newunit=u,file=trim(fileplace)//"photPos.json", status="old", position="append")
!     write(u,"(a)") ","//new_line("a")//'"'//str(counter)//'_'//str(id)//'": '//"["
! else
!     open(newunit=u,file=trim(fileplace)//"photPos.json", status="new")
!     write(u,"(a)") "{"//new_line("a")//'"'//str(counter)//'_'//str(id)//'": '//"["
! end if
! counter = counter + 1

! do while(.not. this%empty())                
!     v = this%pop()
!     if(this%size /= 0)then
!         write(u,"(a,3(es15.8e2,a))")"[",v%x,",",v%y,",",v%z,"],"
!     else
!         write(u,"(a,3(es15.8e2,a))")"[",v%x,",",v%y,",",v%z,"]"
!     end if
! end do
! write(u,"(a)")"]"
! close(u)



! open(newunit=u,file=trim(fileplace)//"photPos.json", status="old", position="append")
! write(u,"(a)") "}"
! close(u)