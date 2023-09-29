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
                new_unittest("Piecewise1D", test_piecewise1D),&
                new_unittest("Piecewise2D", test_piecewise2D)&
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
    
    subroutine test_piecewise2D(error)

        use random, only : init_rng, rang
        use utils,  only : str

        type(error_type), allocatable, intent(out) :: error

        integer ::  ix, iy, u
        integer, parameter :: n=200
        real(kind=wp) :: xx, yy, bin_wid
        integer(kind=int64) :: i
        real(kind=wp) :: xr, yr, diff_sum
        real(kind=wp) :: data2d(n, n), bins(n,n)
        type(piecewise2D) :: obj2D
    
        !set seed
        call init_rng(spread(123456789, 1, 8), .true.)
        
        ! generate image
        data2d = 0.
        bins = 0.
        bin_wid = 2. / n
    
        do i = 1, 10000000
            call rang(xx, yy, 1._wp, 0.1_wp)
            ix = nint(xx / bin_wid)-10
            iy = nint(yy / bin_wid)-30
            if(ix > n .or. ix < 1)cycle
            if(iy > n .or. iy < 1)cycle
            data2d(ix, iy) = data2d(ix, iy) + 1.
            call rang(xx, yy, 1._wp, 0.1_wp)
            ix = nint(xx / bin_wid)+50
            iy = nint(yy / bin_wid)+50
            if(ix > n .or. ix < 1)cycle
            if(iy > n .or. iy < 1)cycle
            data2d(ix, iy) = data2d(ix, iy) + 1.
        end do
        
        data2d = data2d / maxval(data2d)
    
        obj2D = piecewise2D(0.5_wp, 0.5_wp, data2d)
    
        do i = 1, 1000000
            call obj2D%sample(xr, yr)
            ix = nint(xr)+2
            iy = nint(yr)+2
            if(ix > n .or. ix < 1)cycle
            if(iy > n .or. iy < 1)cycle
            bins(ix, iy) = bins(ix, iy) + 1.
        end do

        bins = bins / maxval(bins)

        diff_sum = sum(abs(data2d - bins)) / (n**2)
        if(diff_sum > 1e-2_wp)then
            call test_failed(error, "Piecewise2D check failed!", "Expected a value less than 1e-2! Got "//str(diff_sum, 5))
        end if
        if(allocated(error))return

    end subroutine test_piecewise2D
end module testsPiecewiseMod