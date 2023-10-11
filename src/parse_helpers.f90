module parse_HelpersMod
    
    use constants, only : wp
    use tomlf, only : toml_table, toml_context, toml_array, toml_error, get_value, len
    use tomlf_error, only : make_error
    use vector_class
    
    implicit none
    
    private
    public :: get_vector

contains
    
    type(vector) function get_vector(child, key, error, context, default)
    !! Vector helper function for parsing toml

    !> Input Toml entry to read 
    type(toml_table),   pointer,     intent(in)  :: child
    !> Key to read
    character(*),                    intent(in)  :: key
    !> Default value to assign
    type(vector),       optional,    intent(in)  :: default
    !> Context handle for error reporting
    type(toml_context),              intent(in)  :: context
    !> Error Message
    type(toml_error),   allocatable, intent(out) :: error

    type(toml_array), pointer  :: arr => null()
    real(kind=wp) :: tmp(3)
    integer :: j, origin

    call get_value(child, key, arr, origin=origin)
    if (associated(arr))then
        if(len(arr) /= 3)then
            call make_error(error, &
            context%report("Expected vector of size 3 for"//key, origin, "Wrong vector size"), -1)
            return
        end if
        do j = 1, len(arr)
            call get_value(arr, j, tmp(j))
        end do
        get_vector = vector(tmp(1), tmp(2), tmp(3))
    else
        if(present(default))then
            get_vector = default
        else
            call make_error(error, &
            context%report("Expected vector of size 3 for "//key, origin, "Wrong vector size"), -1)
            return
        end if
    end if
    end function get_vector
end module parse_HelpersMod