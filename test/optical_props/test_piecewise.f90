module testsPiecewiseMod

    use piecewiseMod
    use testdrive, only : new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type, context_t
    use constants, only : wp

    implicit none

    private
    public :: Piecewise_suite

    contains

    subroutine Piecewise_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test Sources", collect_suite1, context)&
                     ]

    end subroutine Piecewise_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Init", test_piecewise1D_init)&
                ! new_unittest("Sample", test_piecewise1D_Sample),&
                ! new_unittest("CDF", test_piecewise1d_CDF)&
                ]

    end subroutine collect_suite1

    subroutine test_piecewise1D_init(error)

        type(error_type), allocatable, intent(out) :: error

        integer :: u
        integer, parameter :: n=129
        real(kind=wp) :: bin_wid
        integer(kind=int64) :: i
        real(kind=wp) :: xr, yr
        real(kind=wp) :: data1d(n, 2)
        type(piecewise1D) :: obj1D
        
        data1d = 0.
        bin_wid = 20. / n

        do i = 1, n
            xr = i*bin_wid
            yr = exp(-(xr-10.)**2 / (2*(0.9)**2))
            data1d(i, 1) = xr
            data1d(i, 2) = yr
        end do

        ! open(newunit=u, file="../res/spectrum.dat")
        ! do i = 1, n
        !     write(u,"(es24.16e3,1x,es24.16e3)") data1d(i, 1),  data1d(i,2)
        ! end do
        ! close(u)
        
        data1d(:,2) = data1d(:,2) / maxval(data1d,2)
        obj1D = piecewise1D(data1d)

        ! open(newunit=u,file="cdf.dat")
        ! do i = 1, size(obj1D%cdf)
        !     write(u,*) obj1D%cdf(i)
        ! end do
        ! close(u)

        open(newunit=u, file="sampled.dat")
        do i = 1, 10000
            call obj1D%sample(xr, yr)
            write(u,*) xr
        end do
        close(u)

        ! call check(error, c%x, 2._wp)
        ! if(allocated(error))return
        ! call check(error, c%y, 2._wp)
        ! if(allocated(error))return
        ! call check(error, c%z, 2._wp)
        ! if(allocated(error))return

    end subroutine test_piecewise1D_init
end module testsPiecewiseMod