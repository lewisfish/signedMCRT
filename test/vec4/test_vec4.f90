module testsVec4Mod

    use vec4_class
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_add", vector_add), &
                new_unittest("Vector_subtract", vector_sub), &
                new_unittest("Vector_multiply", vector_mult), &
                new_unittest("Vector_dot", vector_dot) &
                ]

    end subroutine collect_suite1

    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_add_scal", vector_add_scal), &
                new_unittest("Vector_subtract_scal", vector_sub_scal), &
                new_unittest("Vector_multiply_scal", vector_mult_scal), &
                new_unittest("Vector_div_scal", vector_div_scal) &
                ]

    end subroutine collect_suite2

    subroutine collect_suite3(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_init", vector_init), &
                new_unittest("Vector_sine", vector_sine), &
                new_unittest("Vector_magnitude", vector_mag), &
                new_unittest("Vector_length", vector_length) &
                ]
    end subroutine collect_suite3

    subroutine vector_mag(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, c
        real(kind=wp) :: val

        a = vec4(1._wp, 2._wp, 3._wp, 4._wp)
        c = a%magnitude()

        val = sqrt(30._wp)

        call check(error, c%x, 1._wp / val)
        if(allocated(error))return
        call check(error, c%y, 2._wp / val)
        if(allocated(error))return
        call check(error, c%z, 3._wp / val)
        if(allocated(error))return
        call check(error, c%p, 4._wp / val)
        if(allocated(error))return

    end subroutine vector_mag

    subroutine vector_length(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a
        real(kind=wp) :: c

        a = vec4(1._wp, 2._wp, 3._wp, 4._wp)
        c = a%length()

        call check(error, c, sqrt(30._wp))
        if(allocated(error))return
    end subroutine vector_length

    subroutine vector_sine(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: c
        real(kind=wp) :: pi

        pi = 4._wp*atan(1._wp)
        c = vec4(0._wp, pi/2._wp, pi, (3._wp/2._wp)*pi)
        c = sin(c)

        call check(error, c%x, 0._wp)
        if(allocated(error))return
        call check(error, c%y, 1._wp)
        if(allocated(error))return
        call check(error, c%z, 0._wp)
        if(allocated(error))return
        call check(error, c%p, -1._wp)
        if(allocated(error))return

    end subroutine vector_sine

    subroutine vector_init(error)

        use vector_class, only : vector
        
        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a
        type(vector) :: vec
        real(kind=wp) :: val

        vec = vector(1._wp, 2._wp, 3._wp)
        val = 4._wp
        a = vec4(vec, val)

        call check(error, a%x, 1._wp)
        if(allocated(error))return
        call check(error, a%y, 2._wp)
        if(allocated(error))return
        call check(error, a%z, 3._wp)
        if(allocated(error))return
        call check(error, a%p, 4._wp)
        if(allocated(error))return

    end subroutine vector_init

    subroutine vector_add_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, c
        real(kind=wp) :: b

        a = vec4(1._wp, 1._wp, 1._wp, 1._wp)
        b = 5._wp

        c = a + b
        call check(error, c%x, 6._wp)
        if(allocated(error))return
        call check(error, c%y, 6._wp)
        if(allocated(error))return
        call check(error, c%z, 6._wp)
        if(allocated(error))return
        call check(error, c%p, 6._wp)
        if(allocated(error))return

        c = b + a
        call check(error, c%x, 6._wp)
        if(allocated(error))return
        call check(error, c%y, 6._wp)
        if(allocated(error))return
        call check(error, c%z, 6._wp)
        if(allocated(error))return
        call check(error, c%p, 6._wp)
        if(allocated(error))return

    end subroutine vector_add_scal

    subroutine vector_sub_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, c
        real(kind=wp) :: b

        a = vec4(1._wp, 1._wp, 1._wp, 1._wp)
        b = 5._wp

        c = a - b
        call check(error, c%x, -4._wp)
        if(allocated(error))return
        call check(error, c%y, -4._wp)
        if(allocated(error))return
        call check(error, c%z, -4._wp)
        if(allocated(error))return
        call check(error, c%p, -4._wp)
        if(allocated(error))return

        c = b - a
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return
        call check(error, c%p, 4._wp)
        if(allocated(error))return

    end subroutine vector_sub_scal

    subroutine vector_mult_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, c
        real(kind=wp) :: b

        a = vec4(1._wp, 1._wp, 1._wp, 1._wp)
        b = 5._wp

        c = a * b
        call check(error, c%x, 5._wp)
        if(allocated(error))return
        call check(error, c%y, 5._wp)
        if(allocated(error))return
        call check(error, c%z, 5._wp)
        if(allocated(error))return
        call check(error, c%p, 5._wp)
        if(allocated(error))return

        c = b * a
        call check(error, c%x, 5._wp)
        if(allocated(error))return
        call check(error, c%y, 5._wp)
        if(allocated(error))return
        call check(error, c%z, 5._wp)
        if(allocated(error))return
        call check(error, c%p, 5._wp)
        if(allocated(error))return

    end subroutine vector_mult_scal

    subroutine vector_div_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, c
        real(kind=wp) :: b
        integer :: bi

        a = vec4(10._wp, 10._wp, 10._wp, 10._wp)
        b = 5._wp

        c = a / b
        call check(error, c%x, 2._wp)
        if(allocated(error))return
        call check(error, c%y, 2._wp)
        if(allocated(error))return
        call check(error, c%z, 2._wp)
        if(allocated(error))return
        call check(error, c%p, 2._wp)
        if(allocated(error))return

        bi = 2
        c = a / bi
        call check(error, c%x, 5._wp)
        if(allocated(error))return
        call check(error, c%y, 5._wp)
        if(allocated(error))return
        call check(error, c%z, 5._wp)
        if(allocated(error))return
        call check(error, c%p, 5._wp)
        if(allocated(error))return

    end subroutine vector_div_scal

    subroutine vector_add(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4) :: a, b, c

        a = vec4(1._wp, 1._wp, 1._wp, 1._wp)
        b = vec4(1._wp, 1._wp, 1._wp, 1._wp)

        c = a + b
        call check(error, c%x, 2._wp)
        if(allocated(error))return
        call check(error, c%y, 2._wp)
        if(allocated(error))return
        call check(error, c%z, 2._wp)
        if(allocated(error))return
        call check(error, c%p, 2._wp)
        if(allocated(error))return

    end subroutine vector_add

    subroutine vector_sub(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4) :: a, b, c

        a = vec4(2._wp, 3._wp, 4._wp, 6._wp)
        b = vec4(0.5_wp, 1._wp, 3._wp, 2._wp)

        c = a - b
        call check(error, c%x, 1.5_wp)
        if(allocated(error))return
        call check(error, c%y, 2._wp)
        if(allocated(error))return
        call check(error, c%z, 1._wp)
        if(allocated(error))return
        call check(error, c%p, 4._wp)
        if(allocated(error))return

    end subroutine vector_sub

    subroutine vector_mult(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4) :: a, b, c

        a = vec4(2._wp, 3._wp, 4._wp, 6._wp)
        b = vec4(0.5_wp, 1._wp, 3._wp, 2._wp)

        c = a * b
        call check(error, c%x, 1._wp)
        if(allocated(error))return
        call check(error, c%y, 3._wp)
        if(allocated(error))return
        call check(error, c%z, 12._wp)
        if(allocated(error))return
        call check(error, c%p, 12._wp)
        if(allocated(error))return

    end subroutine vector_mult

    subroutine vector_dot(error)

        type(error_type), allocatable, intent(out) :: error

        type(vec4)    :: a, b
        real(kind=wp) :: c

        a = vec4(2._wp, 3._wp, 4._wp, 6._wp)
        b = vec4(0.5_wp, 1._wp, 3._wp, 2._wp)

        c = a .dot. b
        call check(error, c, 28._wp)
        if(allocated(error))return
    end subroutine vector_dot
end module testsVec4Mod

program test_vector4

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testsVec4Mod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("Suite: Vector .op. vector", collect_suite1), &
                  new_testsuite("Suite: Vector .op. scalar", collect_suite2), &
                  new_testsuite("Suite: Vector functions", collect_suite3) &
                  ]

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat)
    end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_vector4