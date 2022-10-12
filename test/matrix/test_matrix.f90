module testsMatrixMod

    use mat_class
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Matrix_add", Matrix_add) &
                ! ! new_unittest("Vector_subtract", vector_sub), &
                ! ! new_unittest("Vector_multiply", vector_mult), &
                ! ! new_unittest("Vector_dot", vector_dot) &
                ]

    end subroutine collect_suite1

    subroutine matrix_add(error)

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: a(4,4), b, c(4,4), val
        integer :: i, j

        a(:, 1) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a(:, 2) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a(:, 3) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a(:, 4) = [1._wp, 2._wp, 3._wp, 4._wp] 
        b = 5._wp

        c = a + b
        val = 6._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c(i,j), val)
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do
    end subroutine matrix_add
end module testsMatrixMod
program test_matrix

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testsMatrixMod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("Suite: Matrix .op. scalar", collect_suite1) &
    !               new_testsuite("Suite: Vector .op. scalar", collect_suite2), &
    !               new_testsuite("Suite: Vector functions", collect_suite3), &
    !               new_testsuite("Suite: Vector .op. matrix", collect_suite4) &
                  ]

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat)
    end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_matrix