module testScatterMod

    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none
    
    contains
    
    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Scatter_test", Scatter_test)]

    end subroutine collect_suite1

    subroutine Scatter_test(error)
    ! run the test case of a sphere radius 1, tau=10cm^-1
    ! expect to follow theory of N_scatt \approx tau^2 / 2 + tau
        use kernels, only : pathlength_scatter

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: val
        integer :: u

        call pathlength_scatter("scat_test.toml")
        open(newunit=u,file="nscatt.dat")
        read(u,*)val
        close(u)

        call check(error, val, 57.5_wp, thr=0.5_wp)
        if(allocated(error))return

    end subroutine Scatter_test
end module testScatterMod
program test_scat

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testScatterMod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("Suite: Scatter test", collect_suite1)]

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat)
    end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_scat