module random

    implicit none

    private
    public  :: ran2, ranu, rang, init_rng

    contains

        subroutine init_rng(input_seed)
        ! initiate RNG state with reproducable state

            implicit none
            
            integer, optional, intent(IN) :: input_seed

            integer, allocatable :: seed(:)
            integer              :: n, i
            real                 :: a

            call random_seed(size=n)
            allocate(seed(n))
            
            if(present(input_seed))then
                seed = input_seed
            else
                seed = 1234567
            end if

            call random_seed(put=seed)

            !fast forward rng state 100 times to avoid any potential bad seeds
            do i = 1, 100
                a = ran2()
            end do

        end subroutine init_rng

        real function ran2()
        !wrapper for call random number

            implicit none

            call random_number(ran2)

        end function ran2

        real function ranu(a, b)
        !uniformly sample in range[a, b)

            implicit none

            real,    intent(IN)    :: a, b

            ranu = a + ran2() * (b - a)

        end function ranu

        subroutine rang(x, y, avg, sigma)
        ! sample a 2D Guassian distribution

            implicit none

            real,    intent(IN)    :: avg, sigma
            real,    intent(OUT)   :: x,y
            
            real :: s, tmp

            s = 1.

            do while(s >= 1.)
                x = ranu(-1., 1.)
                y = ranu(-1., 1.)
                s = y**2 + x**2
            end do

            tmp = x*sqrt(-2.*log(s)/s)
            x = avg + sigma*tmp

            tmp = y*sqrt(-2.*log(s)/s)
            y = avg + sigma*tmp

        end subroutine rang
end module random