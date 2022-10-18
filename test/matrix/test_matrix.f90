module testsMatrixMod

    use mat_class
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Matrix_add", Matrix_add), &
                new_unittest("Matrix_subtract", matrix_sub), &
                new_unittest("Matrix_multiply", matrix_mult), &
                new_unittest("Matrix_div", matrix_div) &
                ]

    end subroutine collect_suite1


    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Matrix_invert", Matrix_inv), &
                new_unittest("Matrix_init", matrix_init) &
                ]

    end subroutine collect_suite2

    subroutine matrix_add(error)

        type(error_type), allocatable, intent(out) :: error

        type(mat) :: a, c
        real(kind=wp) :: b, val
        integer :: i, j

        a%vals(:, 1) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 2) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 3) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 4) = [1._wp, 2._wp, 3._wp, 4._wp] 
        b = 5._wp

        c = a + b
        val = 6._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do

        c = b + a
        val = 6._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do
    end subroutine matrix_add

    subroutine matrix_sub(error)

        type(error_type), allocatable, intent(out) :: error

        type(mat) :: a, c
        real(kind=wp) :: b, val
        integer :: i, j

        a%vals(:, 1) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 2) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 3) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 4) = [1._wp, 2._wp, 3._wp, 4._wp] 
        b = 1._wp

        c = a - b
        val = 0._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do
    end subroutine matrix_sub

    subroutine matrix_mult(error)

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: b, val
        type(mat) :: a, c
        integer :: i, j

        a%vals(:, 1) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 2) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 3) = [1._wp, 2._wp, 3._wp, 4._wp] 
        a%vals(:, 4) = [1._wp, 2._wp, 3._wp, 4._wp] 
        b = 2._wp

        c = a * b
        val = 2._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 2._wp
        end do

        c = b * a
        val = 2._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 2._wp
        end do
    end subroutine matrix_mult

    subroutine matrix_div(error)

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: b, val
        type(mat) :: a, c
        integer :: i, j

        a%vals(:, 1) = [2._wp, 4._wp, 6._wp, 8._wp] 
        a%vals(:, 2) = [2._wp, 4._wp, 6._wp, 8._wp] 
        a%vals(:, 3) = [2._wp, 4._wp, 6._wp, 8._wp] 
        a%vals(:, 4) = [2._wp, 4._wp, 6._wp, 8._wp] 
        b = 2._wp

        c = a / b
        val = 1._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, c%vals(i,j), val)
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do
    end subroutine matrix_div

    subroutine matrix_inv(error)

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: a(4,4), c(4,4), b(4,4), val
        integer :: i, j

        a(:, 1) = [4._wp, 0._wp, 2._wp, 1._wp] 
        a(:, 2) = [0._wp, 0._wp, 2._wp, 0._wp] 
        a(:, 3) = [0._wp, 1._wp, 2._wp, 0._wp] 
        a(:, 4) = [1._wp, 0._wp, 0._wp, 1._wp] 

        val = 1._wp/3._wp
        b(:, 1) = [val, -1._wp*val, 0._wp, -1._wp*val]
        val = 1._wp - epsilon(1._wp)
        b(:, 2) = [0._wp, -1._wp*val, val, 0._wp]
        val = 0.5_wp - epsilon(1._wp)
        b(:, 3) = [0.0_wp, val, 0._wp, 0._wp] 
        val = 1._wp/3._wp
        b(:, 4) = [-1._wp*val, val, 0._wp, 1._wp + val] 

        c = invert(a)
        do i = 1, 4
            do j = 1, 4
                call check(error, c(i,j), b(i,j))
                if(allocated(error))return
            end do
            val = val + 1._wp
        end do
    end subroutine matrix_inv

    subroutine matrix_init(error)

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: val, b(16)
        type(mat) :: a
        integer :: i, j

        b = [1._wp, 5._wp, 9._wp, 13._wp, &
        2._wp, 6._wp, 10._wp, 14._wp, &
        3._wp, 7._wp, 11._wp, 15._wp, &
        4._wp, 8._wp, 12._wp, 16._wp]

        a = mat(b)

        val = 1._wp
        do i = 1, 4
            do j = 1, 4
                call check(error, a%vals(i,j), val)
                if(allocated(error))return
                val = val + 1._wp
            end do
        end do
    end subroutine matrix_init

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

    testsuites = [new_testsuite("Suite: Matrix .op. scalar", collect_suite1), &
                  new_testsuite("Suite: Matrix Other", collect_suite2) &
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