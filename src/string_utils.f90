module string_utils

    use iso_fortran_env, only : sp => real32, dp => real64

    implicit none
    
    !functions to turn variables into strings
    interface str
        module procedure str_I32
        module procedure str_I64
        module procedure str_Iarray
        module procedure str_R4
        module procedure str_R8
        module procedure str_R8array
        module procedure str_logical
        module procedure str_logicalarray
        module procedure str_vector
    end interface str

    private
    public :: str

contains
    function str_I32(i, len)

        use iso_fortran_env, only : Int32

        implicit none

        integer(int32),    intent(IN)    :: i
        integer, optional, intent(IN) :: len

        character(len=:), allocatable :: str_I32
        character(len=100) :: string
        integer            :: lentmp, lenuse

        write(string,'(I100.1)') I

        if(present(len))then
            lentmp = len_trim(adjustl(string))
            lenuse = len

            if(len >= lentmp)then
                str_I32 = repeat("0", lenuse - lentmp)//trim(adjustl(string))
            else
                str_I32 = trim(adjustl(string))
                str_I32 = trim(adjustl(str_I32(:len)))                        
            end if
        else
            str_I32 = trim(adjustl(string))
        end if
    end function str_I32


    function str_I64(i, len)

        use iso_fortran_env, only : Int64

        implicit none

        integer(int64),    intent(IN)    :: i
        integer, optional, intent(IN) :: len

        character(len=:), allocatable :: str_I64
        character(len=100) :: string
        integer            :: lentmp, lenuse

        write(string,'(I100.1)') I

        if(present(len))then
            lentmp = len_trim(adjustl(string))
            lenuse = len

            if(len >= lentmp)then
                str_I64 = repeat("0", lenuse - lentmp)//trim(adjustl(string))
            else
                str_I64 = trim(adjustl(string))
                str_I64 = trim(adjustl(str_I64(:len)))                        
            end if
        else
            str_I64 = trim(adjustl(string))
        end if
    end function str_I64


    function str_iarray(i)

        implicit none

        integer, intent(IN) :: i(:)

        character(len=:), allocatable :: str_iarray
        character(len=100) :: string
        integer :: k, j, length

        length = 3*size(i)-1
        str_iarray = repeat(" ", length)

        k = 1
        do j = 1, size(i)
            write(string,'(I2.2)') I(j)
            if(j == 1)then
                str_iarray(k:k+2) = trim(adjustl(string))
                k = k + 2
            else
                str_iarray(k:k+2) = ':'//trim(adjustl(string))
                k = k + 3
            end if
        end do

    end function str_iarray


    function str_R4(i, len)

        implicit none

        real(kind=sp), intent(IN) :: i
        integer, optional, intent(IN) :: len

        character(len=:), allocatable :: str_R4
        character(len=100) :: string

        write(string,'(f100.8)') I

        if(present(len))then
            str_R4 = trim(adjustl(string))
            str_R4 = trim(adjustl(str_R4(:len)))
        else
            str_R4 = trim(adjustl(string))
        end if
    end function str_r4

    function str_R8(i, len)

        implicit none

        real(kind=dp),     intent(IN) :: i
        integer, optional, intent(IN) :: len

        character(len=:), allocatable :: str_R8
        character(len=100) :: string

        write(string,'(f100.16)') I

        if(present(len))then
            str_R8 = trim(adjustl(string))
            str_R8 = trim(adjustl(str_R8(:len)))
        else
            str_R8 = trim(adjustl(string))
        end if
    end function str_R8


    function str_R8array(a, width)

        implicit none

        real(kind=dp),     intent(IN) :: a(:)
        integer,           intent(IN) :: width

        character(len=:), allocatable :: str_R8array
        integer :: i, length, lens(size(a)), k

        length = 0
        lens = 0
        do i = 1, size(a)
            lens(i) = len(str_R8(a(i), width))
            length = length + lens(i)
        end do
        length = length + size(a)

        str_R8array = repeat(" ", length)
        k = 1
        do i = 1, size(a)
            str_R8array(k:k+lens(i)) = str_R8(a(i), width)//" "
            k = k + lens(i)+1
        end do
    end function str_R8array


    function str_logical(a)

        implicit none

        logical, intent(IN) :: a

        character(len=:), allocatable :: str_logical
        character(len=100) :: string

        write(string,'(L1)') a
        str_logical = trim(adjustl(string))

    end function str_logical


    function str_logicalarray(a)

        implicit none

        logical, intent(IN) :: a(:)

        character(len=:), allocatable :: str_logicalarray
        character(len=100) :: string
        integer :: i

        do i = 1, size(a)
            write(string,'(L1)') a(i)
            str_logicalarray = str_logicalarray//' '//trim(adjustl(string))
        end do

    end function str_logicalarray

    function str_vector(a)

        use vector_class

        type(vector), intent(IN) :: a

        character(len=:), allocatable :: str_vector
        character(len=100) :: string

        write(string,'(f100.16)') a%x
        str_vector = trim(adjustl(string))//","
        write(string,'(f100.16)') a%y
        str_vector = str_vector//trim(adjustl(string))//","
        write(string,'(f100.16)') a%z
        str_vector = str_vector//trim(adjustl(string))

    end function str_vector
end module string_utils