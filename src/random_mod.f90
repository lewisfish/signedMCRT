module random
! module provides an interface to call random_numbers and various other random distributions
!
!
    use vector_class
    use constants, only : wp

    implicit none

    private
    public  :: ran2, ranu, rang, randint, init_rng

    contains

        subroutine init_rng(input_seed, fwd)
        ! initiate RNG state with reproducable state

            implicit none
            
            integer, optional, intent(IN) :: input_seed(:)
            logical, optional, intent(IN) :: fwd

            integer, allocatable :: seed(:)
            integer              :: n, i
            logical              :: ffwd
            real(kind=wp)        :: a


            call random_seed(size=n)
            allocate(seed(n))

            if(present(input_seed))then
                seed = 0
                seed = input_seed
            else
                seed = 1234567
            end if

            if(present(fwd))then
                ffwd = fwd
            else
                ffwd = .true.
            end if

            call random_seed(put=seed)

            !fast forward rng state 100 times to avoid any potential bad seeds
            if(ffwd)then
                call random_seed(get=seed)
                do i = 1, 100
                    a = ran2()
                    call random_seed(get=seed)
                end do
            end if
        end subroutine init_rng

        function ran2() result(res)
        !wrapper for call random number

            implicit none

            real(kind=wp) :: res

            call random_number(res)

        end function ran2

        function ranu(a, b) result(res)
        !uniformly sample in range[a, b)

            implicit none

            real(kind=wp) :: res
            real(kind=wp), intent(IN) :: a, b

            res = a + ran2() * (b - a)

        end function ranu

        subroutine rang(x, y, avg, sigma)
        ! sample a 2D Guassian distribution

            implicit none

            real(kind=wp), intent(IN)  :: avg, sigma
            real(kind=wp), intent(OUT) :: x,y
            
            real(kind=wp) :: s, tmp

            s = 1._wp

            do while(s >= 1._wp)
                x = ranu(-1._wp, 1._wp)
                y = ranu(-1._wp, 1._wp)
                s = y**2 + x**2
            end do

            tmp = x*sqrt(-2._wp*log(s)/s)
            x = avg + sigma*tmp

            tmp = y*sqrt(-2._wp*log(s)/s)
            y = avg + sigma*tmp

        end subroutine rang

        integer function randint(a, b)

            implicit none

            integer, intent(IN) :: a, b

            randint = a + floor((b + 1 - a)*ran2()) 

        end function randint

end module random

! Program test
    
!     use random, only : randint

!     implicit none
    
!     integer :: i

!     do i = 1, 100
!         print*,randint(0, 5)
!     end do

! end program test