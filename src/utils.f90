module utils
! module provides various utility functions for simulation
! progress bar, var to str converters
! colours for pretty output
    
    ! use constants, only : wp
    use iso_fortran_env, only : sp => real32, dp => real64

    implicit none

    !foreground colours
    character(len=2), parameter :: black   = '30', &
                                   red     = '31', &
                                   green   = '32', &
                                   yellow  = '33', &
                                   blue    = '34', &
                                   magenta = '35', &
                                   cyan    = '36', &
                                   white   = '37'

    !background colours
    character(len=2), parameter :: black_b   = '40', &
                                   red_b     = '41', &
                                   green_b   = '42', &
                                   yellow_b  = '43', &
                                   blue_b    = '44', &
                                   magenta_b = '45', &
                                   cyan_b    = '46', &
                                   white_b   = '47'

    !styles
    character(len=2), parameter :: bold          = '01', &
                                   italic        = '03', &
                                   underline     = '04', &
                                   inverse       = '07', &
                                   strikethrough = '09'

    !ANSI control characters                               
    character(len=2), parameter :: start = achar(27)//'['
    character(len=3), parameter :: end = '[0m'


    !functions to add colour to output via ANSI colour codes
    interface colour
        module procedure colour_char
        module procedure colour_int
        module procedure colour_real4
        module procedure colour_real8
    end interface

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
    end interface str

    type :: pbar
        integer       :: iters, current_iter, time_remaing(3), time_taken(3), threads
        real(kind=dp) :: percentage, start_t, start_tt, finish_t, average
        logical       :: first
        contains
            procedure :: progress => progress_sub
    end type pbar

    interface pbar
        module procedure :: init_pbar_func
    end interface pbar

    !subroutines to swap variables
    interface swap
        module procedure swap_I
        module procedure swap_R4
        module procedure swap_R8
    end interface swap

    !subroutines to clamp variables
    interface clamp
        module procedure clamp_R4
        module procedure clamp_R8
    end interface clamp

    !subroutines to mix variables
    interface mix
        module procedure mix_R4
        module procedure mix_R8
    end interface mix

    !subroutines to deg2rad variables
    interface deg2rad
        module procedure deg2rad_R4
        module procedure deg2rad_R8
    end interface deg2rad

    !subroutines to rad2deg variables
    interface rad2deg
        module procedure rad2deg_R4
        module procedure rad2deg_R8
    end interface rad2deg


    interface
        function c_chdir(path) bind(C, name="chdir")
                
            use iso_c_binding

            character(kind=c_char), intent(IN) :: path(*)
            integer(kind=C_int) :: c_chdir
        end function c_chdir
    end interface

    private
    public :: str, swap, colour, mem_free, chdir, mix, clamp, rad2deg, deg2rad, lerp, get_time, print_time
    public :: bold, italic, underline, strikethrough, black, red, green, yellow, blue, magenta, cyan, white
    public :: black_b, red_b, green_b, yellow_b, blue_b, magenta_b, cyan_b, white_b, pbar

    contains

        pure function lerp(t, v1, v2) result(res)

            implicit none

            real(kind=dp) :: res
            real(kind=dp), intent(IN) :: v1, v2, t

            res = (1._dp - t) * v1 + t * v2

        end function lerp

        type(pbar) function init_pbar_func(n)

            use omp_lib

            implicit none

            integer, intent(IN) :: n

#ifdef _OPENMP
            init_pbar_func%threads = omp_get_max_threads()
#else
            init_pbar_func%threads = 1
#endif
            init_pbar_func%iters = n
            init_pbar_func%current_iter = 0
            init_pbar_func%time_remaing = 0 
            init_pbar_func%time_taken = 0
            init_pbar_func%percentage = 0.0 
            init_pbar_func%start_t = 0.0
            init_pbar_func%start_tt = 0.0
            init_pbar_func%finish_t = 0.0  
            init_pbar_func%average = 0.0
            init_pbar_func%first = .true.

        end function init_pbar_func


        subroutine progress_sub(this)

            use iso_fortran_env, only : output_unit

            implicit none

            class(pbar) :: this
            integer           :: width
            character(len=52) :: line
            real(kind=dp)     :: time

!$omp critical
            if(.not. this%first)then
                call cpu_time(this%finish_t)
                this%average = this%average + (this%finish_t - this%start_t)
                time = this%average / real(this%threads * this%current_iter)
                time = time * (this%iters - this%current_iter)
                this%time_remaing(1) = floor(time / (60*60))
                this%time_remaing(2) = floor(mod(time / 60, 60._dp))
                this%time_remaing(3) = int(mod(time, 60._dp))

                time = (this%finish_t - this%start_tt) / this%threads
                this%time_taken(1) = floor(time / (60*60))
                this%time_taken(2) = floor(mod(time / 60, 60._dp))
                this%time_taken(3) = int(mod(time, 60._dp))
            else    
                this%first = .false.
                call cpu_time(this%start_tt)
            end if


            this%current_iter = this%current_iter + 1
            if(this%current_iter <= this%iters)then
                this%percentage = 100._dp*real(this%current_iter) / real(this%iters)

                width = int(this%percentage/ 2._dp)
                line = "[" // repeat("#", width) // repeat(" ", 50 - width) // "]"

                write(unit=output_unit,fmt='(A)',advance="no") start//"1000D"//line//" "//str(int(this%percentage),3)//"%  ["//&
                str(this%time_taken)//"<"//str(this%time_remaing)//"]"
                if(this%percentage >= 100._dp)write(unit=output_unit,fmt='(A)')new_line("a")
                flush(output_unit)
            end if
!$omp end critical
            call cpu_time(this%start_t)

        end subroutine progress_sub


        pure function rad2deg_R4(angle) result(res)

            use constants, only : PI

            implicit none

            real(kind=sp), intent(IN) :: angle
            real(kind=sp) :: res

            res = angle/PI*180._sp

        end function rad2deg_R4

        pure function rad2deg_R8(angle) result(res)

            use constants, only : PI

            implicit none

            real(kind=dp), intent(IN) :: angle
            real(kind=dp) :: res

            res = angle/PI*180._dp

        end function rad2deg_R8


        pure function deg2rad_R4(angle) result(res)

            use constants, only : PI

            implicit none

            real(kind=sp), intent(IN) :: angle
            real(kind=sp) :: res

            res = angle*PI/180._sp

        end function deg2rad_R4

        pure function deg2rad_R8(angle) result(res)

            use constants, only : PI

            implicit none

            real(kind=dp), intent(IN) :: angle
            real(kind=dp) :: res

            res = angle*PI/180._dp

        end function deg2rad_R8

        subroutine chdir(path, error)

            use iso_c_binding, only : c_null_char

            implicit none

            character(*), intent(IN)       :: path
            integer, optional, intent(OUT) :: error

            integer :: err

            err = c_chdir(trim(path)//c_null_char)
            if(present(error))error = err
        end subroutine chdir

        pure function clamp_R4(val, lo, hi) result(res)

            implicit none

            real(kind=dp) :: res
            real(kind=dp), intent(IN) :: val, hi, lo

            if(val < lo)then
                res = lo
            elseif(val > hi)then
                res = hi
            else
                res = val
            end if

        end function clamp_R4

        pure function clamp_R8(val, lo, hi) result(res)

            implicit none

            real(kind=sp) :: res
            real(kind=sp), intent(IN) :: val, hi, lo

            if(val < lo)then
                res = lo
            elseif(val > hi)then
                res = hi
            else
                res = val
            end if

        end function clamp_R8

        pure function mix_R4(x, y, a) result(res)

            implicit none

            real(kind=sp) :: res
            real(kind=sp), intent(IN) :: x, y, a

            res = x*(1._sp - a) + y*a

        end function mix_R4

        pure function mix_R8(x, y, a) result(res)

            implicit none

            real(kind=dp) :: res
            real(kind=dp), intent(IN) :: x, y, a

            res = x*(1._dp - a) + y*a

        end function mix_R8


        subroutine swap_I(a, b)

            implicit none

            integer, intent(INOUT) :: a, b
            integer :: tmp

            tmp = a
            a = b
            b = tmp

        end subroutine swap_I


        subroutine swap_R4(a, b)

            implicit none

            real(kind=sp), intent(INOUT) :: a, b
            real(kind=sp) :: tmp

            tmp = a
            a = b
            b = tmp

        end subroutine swap_R4


        subroutine swap_R8(a, b)

            implicit none

            real(kind=dp), intent(INOUT) :: a, b
            real(kind=dp) :: tmp

            tmp = a
            a = b
            b = tmp

        end subroutine swap_R8


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


        function colour_char(string, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            character(*),           intent(IN) :: string
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5
            character(len=:), allocatable      :: colourised

            colourised = string

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_char


        function colour_int(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            integer,                intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(I50.1)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_int


        function colour_real4(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            real(kind=sp),          intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(F50.8)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_real4


        function colour_real8(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            real(kind=dp),          intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(F50.16)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_real8


        function mem_free()

            use iso_fortran_env, only : int64 !as numbers are large

            implicit none

            integer(int64) :: mem_free

            integer(int64)    :: i
            character(len=15) :: tmp
            integer           :: u

            open(newunit=u,file='/proc/meminfo',status='old')

            read(u,*)tmp, i
            read(u,*)tmp, i
            read(u,*)tmp, i
            
            mem_free = i * 1024_int64 !convert from Kib to b 
        end function mem_free

        real function get_time()

#ifdef _OPENMP
            use omp_lib
#endif
            implicit none

#ifdef _OPENMP
                get_time = omp_get_wtime()
#else
                call cpu_time(get_time)
#endif

        end function get_time

        subroutine print_time(time, id)

            use constants, only : wp

            implicit none

            real(kind=wp), intent(IN) :: time
            integer,       intent(IN) :: id

            if(id == 0)then
                if(time >= 60._wp)then
                   print*, floor((time)/60._wp),"mins", mod(time, 60._wp)/100._wp,"s"
                else
                   print*, 'time taken ~',time,'s'
                end if
            end if

        end subroutine print_time

end module utils