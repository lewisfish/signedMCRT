module testsPhotonMod

    use photonMod
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        ! testsuite = [ &
                ! new_unittest("Vector_add", vector_add) &
                ! ! new_unittest("Vector_subtract", vector_sub), &
                ! ! new_unittest("Vector_multiply", vector_mult), &
                ! ! new_unittest("Vector_dot", vector_dot) &
                ! ]

    end subroutine collect_suite1

end module testsPhotonMod
program test_photon

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testsPhotonMod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    ! testsuites = [new_testsuite("Suite: Vector .op. vector", collect_suite1), &
    !               new_testsuite("Suite: Vector .op. scalar", collect_suite2), &
    !               new_testsuite("Suite: Vector functions", collect_suite3), &
    !               new_testsuite("Suite: Vector .op. matrix", collect_suite4) &
    !               ]

    ! do i = 1, size(testsuites)
    !     write(error_unit, fmt) "Testing:", testsuites(i)%name
    !     call run_testsuite(testsuites(i)%collect, error_unit, stat)
    ! end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_photon