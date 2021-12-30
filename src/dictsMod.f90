module dict_mod
! module provides a limited dictionary type for storing simulation parameters
!
!
    
    use constants, only : wp

    implicit none

    type :: dict_t
        type(dict_data), allocatable :: dict(:)
        integer :: count
        contains
            procedure :: add_entry
            ! procedure :: get_value
            procedure :: get_value_str
            procedure :: get_value_real
    end type dict_t

    interface dict_t
        module procedure init_dict
    end interface dict_t

    type :: dict_data
        character(len=64) :: key
        class(*),          allocatable :: value
    end type dict_data

    contains
    
    function init_dict(size)

        implicit none

        integer, intent(IN) :: size
        type(dict_t) :: init_dict

        allocate(init_dict%dict(size))
        init_dict%count = 1

    end function init_dict

    subroutine add_entry(this, key, value)
        
        implicit none
        
        class(dict_t) :: this
        class(*),     intent(IN) :: value
        character(*), intent(IN) :: key

        this%dict(this%count)%key = key
        allocate(this%dict(this%count)%value, source=value)
        this%count = min(this%count + 1, size(this%dict))

    end subroutine add_entry


    ! function get_value(this, key)

    !     implicit none

    !     class(dict_t) :: this
    !     character(*), intent(IN) :: key
    !     class(*), allocatable :: get_value

    !     integer :: i, pos

    !     do i = 1, size(this%dict)
    !         pos = index(this%dict(i)%key, key)
    !         if(pos > 0)then
    !             get_value = this%dict(i)%value
    !             return
    !         end if
    !     end do
    !     Error stop "No such key!"

    ! end function get_value

    function get_value_str(this, key)

        use utils, only : str
        use iso_fortran_env, only : sp => real32, dp => real64
        implicit none

        class(dict_t) :: this
        character(*), intent(IN) :: key
        class(*), allocatable :: value
        character(len=:), allocatable :: get_value_str

        integer :: i, pos

        !this is very inefficient, but will suffice for small dicts like intended
        !if a large dict is required then will need to use hash tables I think...
        do i = 1, size(this%dict)
            pos = index(this%dict(i)%key, key)
            if(pos > 0)then
                value = this%dict(i)%value
                    select type (value)
                    type is (character(*))
                        get_value_str = value
                    type is(real(sp))
                        get_value_str = str(real(value, kind=sp))
                    type is(integer)
                        get_value_str = str(int(value))
                    type is(real(dp))
                        get_value_str = str(real(value, kind=dp))
                    type is(logical)
                        get_value_str = str(value)
                    class default
                        continue
                    end select
                return
            end if
        end do
        Error stop "No such key!"

    end function get_value_str

    function get_value_real(this, key)

        implicit none

        class(dict_t) :: this
        character(*), intent(IN) :: key
        class(*), allocatable :: value
        real(kind=wp) :: get_value_real

        integer :: i, pos

        get_value_real = -99._wp

        !this is very inefficient, but will suffice for small dicts like intended
        !if a large dict is required then will need to use hash tables I think...
        do i = 1, size(this%dict)
            pos = index(this%dict(i)%key, key)
            if(pos > 0)then
                value = this%dict(i)%value
                    select type (value)
                    type is (character(*))
                        error stop "cant convert string to real"
                    type is(integer)
                        get_value_real = real(value,kind=wp)
                    type is(real(wp))
                        get_value_real = real(value,kind=wp)
                    type is(logical)
                        error stop "cant convert logical to real"                        
                    class default
                        continue
                    end select
                return
            end if
        end do
        Error stop "No such key!"

    end function get_value_real

end module dict_mod