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
                new_unittest("Scatter_test1", Scatter_test1),&
                new_unittest("Scatter_test2", Scatter_test2)]

    end subroutine End_to_End_suite

    subroutine Scatter_test1(error)
    ! run the test case of a sphere radius 1, tau=10cm^-1
    ! expect to follow theory of N_scatt \approx tau^2 / 2 + tau
        use kernels, only : test_kernel

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: val
        integer :: u

        call test_kernel("scat_test.toml", .false.)
        open(newunit=u,file="nscatt.dat")
        read(u,*)val
        close(u)

        call check(error, val, 57.5_wp, thr=0.5_wp)
        if(allocated(error))return

    end subroutine Scatter_test1

    subroutine Scatter_test2(error)
        ! run the test case of pencil beam in +z in ~infinite medium mua=0, mus=1cm^-1, hgg=0.9
            use kernels, only : test_kernel
            use vector_class, only : vector
            type(error_type), allocatable, intent(out) :: error
    
            type(vector) :: pos 
            integer :: i, u
            real(kind=wp) :: vals(8, 3)

            ! values taken from table 7 of Two-step verification method for Monte Carlo codes in biomedical optics applications
            ! 1st moments
            vals(1, :) = [0.0, 0.0, 1.0]
            vals(2, :) = [0.0, 0.0, 1.9]
            vals(3, :) = [0.0, 0.0, 2.71]
            vals(4, :) = [0.0, 0.0, 3.349]
            !2nd moments
            vals(5, :) = [0.0, 0.0, 2.0]
            vals(6, :) = [0.1266666666, 0.1266666666, 5.5466666666]
            vals(7, :) = [0.469933, 0.469933, 10.28013]
            vals(8, :) = [1.091246, 1.091246, 15.91551]

            call test_kernel("scat_test2.toml", .true.)
            open(newunit=u,file="positions.dat")
            do i = 1, 8
                read(u,*)pos%x, pos%y, pos%z
                call check(error, pos%x, vals(i, 1), thr=0.1_wp)
                if(allocated(error))then
                    print*,"Error in x, scatter order:", i
                    return
                end if
                call check(error, pos%y, vals(i, 2), thr=0.1_wp)
                if(allocated(error))then
                    print*,"Error in y, scatter order:", i
                    return
                end if
                call check(error, pos%z, vals(i, 3), thr=0.143_wp)
                if(allocated(error))then
                    print*,"Error in z, scatter order:", i
                    return
                end if
            end do
            close(u)
        end subroutine Scatter_test2
end module testScatterMod