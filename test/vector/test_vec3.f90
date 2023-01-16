module testsVecMod

    use vector_class, only : vector, max, abs, nint, min
    use testdrive, only : new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type
    use constants, only : wp

    implicit none

    private
    public :: Vector_suite

    contains

    subroutine Vector_suite(testsuites)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)

        testsuites = [new_testsuite("Vector .op. vector", collect_suite1),&
                      new_testsuite("Vector .op scalar", collect_suite2),&
                      new_testsuite("Vector functions", collect_suite3),&
                      new_testsuite("Vector .op. matrix", collect_suite4) &
                     ]

    end subroutine Vector_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_add", vector_add), &
                new_unittest("Vector_subtract", vector_sub), &
                new_unittest("Vector_multiply", vector_mult), &
                new_unittest("Vector_dot", vector_dot), &
                new_unittest("Vector_cross", vector_cross), &
                new_unittest("Vector_equal", vector_equal) &
                ]

    end subroutine collect_suite1

    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_add_scalar", vector_add_scal), &
                new_unittest("Vector_subtract_scalar", vector_sub_scal), &
                new_unittest("Vector_multiply_scalar", vector_mult_scal), &
                new_unittest("Vector_div_scalar", vector_div_scal), &
                new_unittest("Vector_div_scalar_int", vector_div_scal_int) &
                ]

    end subroutine collect_suite2

    subroutine collect_suite3(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_magnitude", vector_mag), &
                new_unittest("Vector_length", vector_length), &
                new_unittest("Vector_nint", vector_nint), &
                new_unittest("Vector_abs", vector_abs), &
                new_unittest("Vector_min", vector_min), &
                new_unittest("Vector_max", vector_max) &
                ]

    end subroutine collect_suite3

    subroutine collect_suite4(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Vector_dot_mat", vector_dot_mat) &
                ]

    end subroutine collect_suite4

    subroutine vector_add(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b, c

        a = vector(1._wp, 1._wp, 1._wp)
        b = vector(1._wp, 1._wp, 1._wp)

        c = a + b
        call check(error, c%x, 2._wp)
        if(allocated(error))return
        call check(error, c%y, 2._wp)
        if(allocated(error))return
        call check(error, c%z, 2._wp)
        if(allocated(error))return

    end subroutine vector_add

    subroutine vector_sub(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b, c

        a = vector(2._wp, 2._wp, 2._wp)
        b = vector(1._wp, 1._wp, 1._wp)

        c = a - b
        call check(error, c%x, 1._wp)
        if(allocated(error))return
        call check(error, c%y, 1._wp)
        if(allocated(error))return
        call check(error, c%z, 1._wp)
        if(allocated(error))return

    end subroutine vector_sub

    subroutine vector_mult(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b, c

        a = vector(2._wp, 2._wp, 2._wp)
        b = vector(2._wp, 2._wp, 2._wp)

        c = a * b
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return

    end subroutine vector_mult

    subroutine vector_dot(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b
        real(kind=wp) :: c

        a = vector(1._wp, 2._wp, 3._wp)
        b = vector(6._wp, 5._wp, 4._wp)

        c = a .dot. b
        call check(error, c, 28._wp)
        if(allocated(error))return

    end subroutine vector_dot

    subroutine vector_cross(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b, c

        a = vector(1._wp, 2._wp, 3._wp)
        b = vector(4._wp, 5._wp, 6._wp)

        c = a .cross. b
        call check(error, c%x, -3._wp)
        if(allocated(error))return
        call check(error, c%y, 6._wp)
        if(allocated(error))return
        call check(error, c%z, -3._wp)
        if(allocated(error))return

    end subroutine vector_cross


    subroutine vector_equal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, b
        logical :: c

        a = vector(1._wp, 2._wp, 3._wp)
        b = vector(1._wp, 2._wp, 3._wp)

        c = (a == b)
        call check(error, c, .true.)
        if(allocated(error))return

        b = vector(2._wp, 3._wp, 4._wp)
        c = (a == b)
        call check(error, c, .false.)
        if(allocated(error))return
    end subroutine vector_equal

    subroutine vector_add_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c
        real(kind=wp) :: b

        a = vector(1._wp, 1._wp, 1._wp)
        b = 3._wp

        c = a + b
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return

        c = b + a
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return

    end subroutine vector_add_scal

    subroutine vector_sub_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c
        real(kind=wp) :: b
        
        a = vector(2._wp, 2._wp, 2._wp)
        b = 1._wp

        c = a - b
        call check(error, c%x, 1._wp)
        if(allocated(error))return
        call check(error, c%y, 1._wp)
        if(allocated(error))return
        call check(error, c%z, 1._wp)
        if(allocated(error))return

        c = b - a
        call check(error, c%x, -1._wp)
        if(allocated(error))return
        call check(error, c%y, -1._wp)
        if(allocated(error))return
        call check(error, c%z, -1._wp)
        if(allocated(error))return

    end subroutine vector_sub_scal

    subroutine vector_mult_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c
        real(kind=wp) :: b

        a = vector(2._wp, 2._wp, 2._wp)
        b = 2._wp

        c = a * b
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return

        c = b * a
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return

    end subroutine vector_mult_scal

    subroutine vector_div_scal(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c
        real(kind=wp) :: b

        a = vector(8._wp, 8._wp, 8._wp)
        b = 2._wp

        c = a / b
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return
    end subroutine vector_div_scal

    subroutine vector_div_scal_int(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c
        integer      :: b

        a = vector(8._wp, 8._wp, 8._wp)
        b = 2

        c = a / b
        call check(error, c%x, 4._wp)
        if(allocated(error))return
        call check(error, c%y, 4._wp)
        if(allocated(error))return
        call check(error, c%z, 4._wp)
        if(allocated(error))return
    end subroutine vector_div_scal_int

    subroutine vector_mag(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c

        a = vector(3._wp, 3._wp, 3._wp)
        c = a%magnitude()

        call check(error, c%x, 3._wp/sqrt(27._wp))
        if(allocated(error))return
        call check(error, c%y, 3._wp/sqrt(27._wp))
        if(allocated(error))return
        call check(error, c%z, 3._wp/sqrt(27._wp))
        if(allocated(error))return
    end subroutine vector_mag

    subroutine vector_length(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a
        real(kind=wp) :: b

        a = vector(2._wp, 2._wp, 2._wp)
        b = a%length()

        call check(error, b, sqrt(12._wp))
        if(allocated(error))return
    end subroutine vector_length

    subroutine vector_abs(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c

        a = vector(-8._wp, -.0009_wp, -80000._wp)
        c = abs(a)

        call check(error, c%x, 8._wp)
        if(allocated(error))return
        call check(error, c%y, 0.0009_wp)
        if(allocated(error))return
        call check(error, c%z, 80000._wp)
        if(allocated(error))return
    end subroutine vector_abs

    subroutine vector_nint(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: a, c

        a = vector(8.123213_wp, 8000.721321_wp, 0.000321_wp)
        c = nint(a)

        call check(error, c%x, 8._wp)
        if(allocated(error))return
        call check(error, c%y, 8001._wp)
        if(allocated(error))return
        call check(error, c%z, 0.0_wp)
        if(allocated(error))return
    end subroutine vector_nint

    subroutine vector_max(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector)  :: a, c
        real(kind=wp) :: b

        a = vector(8._wp, -1._wp, 2.2_wp)
        b = 2._wp

        c = max(a, b)
        call check(error, c%x, 8._wp)
        if(allocated(error))return
        call check(error, c%y, 2._wp)
        if(allocated(error))return
        call check(error, c%z, 2.2_wp)
        if(allocated(error))return
    end subroutine vector_max

    subroutine vector_min(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector)  :: a, c
        real(kind=wp) :: b

        a = vector(8._wp, -1._wp, 2.2_wp)
        b = 2._wp

        c = min(a, b)
        call check(error, c%x, 2._wp)
        if(allocated(error))return
        call check(error, c%y, -1._wp)
        if(allocated(error))return
        call check(error, c%z, 2._wp)
        if(allocated(error))return
    end subroutine vector_min

    subroutine vector_dot_mat(error)

        type(error_type), allocatable, intent(out) :: error

        type(vector)  :: a, c
        real(kind=wp) :: b(4, 4)

        a = vector(8._wp, -1._wp, 2.2_wp)
        b(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
        b(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
        b(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
        b(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

        c = a .dot. b
        call check(error, c%x, 8._wp)
        if(allocated(error))return
        call check(error, c%y, -1._wp)
        if(allocated(error))return
        call check(error, c%z, 2.2_wp)
        if(allocated(error))return
    end subroutine vector_dot_mat
end module testsVecMod