module testScatterMod

    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none
    
    private
    public :: End_to_End_suite

    contains
    
    subroutine End_to_End_suite(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Scatter_test", Scatter_test)]

    end subroutine End_to_End_suite

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