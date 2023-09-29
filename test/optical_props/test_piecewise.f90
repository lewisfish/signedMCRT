module testsPiecewiseMod

    use piecewiseMod
    use testdrive, only : new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type, context_t, test_failed
    use constants, only : wp

    implicit none

    private
    public :: Piecewise_suite

    contains

    subroutine Piecewise_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test Piecewise", collect_suite1, context)&
                     ]

    end subroutine Piecewise_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Piecewise1D", test_piecewise1D)&
                ]

    end subroutine collect_suite1

    subroutine test_piecewise1D(error)

        use random, only : init_rng
        use utils,  only : str

        type(error_type), allocatable, intent(out) :: error

        integer ::  u, idx, i
        integer, parameter :: n=376
        real(kind=wp) :: bin_wid, xr, yr, bins(n), data1d(n, 2), diff_sum
        type(piecewise1D) :: obj1D
    
        !set seed
        call init_rng(spread(123456789, 1, 8), .true.)
        
        data1d = 0.
        bin_wid = (1000.-250.)/n
        bins = 0.
    
        open(newunit=u, file="test/optical_props/blood.dat")
        do i = 1, 376
            read(u,*) data1d(i, :)
        end do
    
        obj1D = piecewise1D(data1d)
    
        do i = 1, 1000000
            call obj1D%sample(xr, yr)
            idx = nint((xr-250) / bin_wid) + 1
            if(idx > 0 .and. idx < 377)bins(idx) = bins(idx) + 1
        end do
        
        data1d(:, 2) = data1d(:, 2) / maxval(data1d(:, 2))
        bins = bins / maxval(bins)

        diff_sum = sum(abs(data1d(:, 2) - bins))

        if(diff_sum > 2.0)then
            call test_failed(error, "Piecewise check failed!", "Expected a value less than 1.0! Got "//str(diff_sum, 5))
        end if
        if(allocated(error))return

    end subroutine test_piecewise1D
end module testsPiecewiseMod