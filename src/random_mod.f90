module random
!! module provides an interface to call random_numbers and various other random distributions=======    !!This module defines a set of functions that return random numbers in different distributions.    !!- ran2. Returns a single float uniformly in the range [0, 1)    !!- ranu. Return a single float uniformly in the range [a, b)    !!- randint. Returns a single integer uniformly in the range [a, b)    !!- rang. Returns a single float from a Gaussian distribution with mean *avg* and std *sigma*.    !!- init_rng. Seeds the internal random number generator with a reproducible seed.

    use vector_class
    use constants, only : wp

    implicit none
    !> Sequence type for quasi-monte carlo
    type :: seq
        !> Current index to get value for.
        integer :: index
        !> Base from which to calculate radical inverse from.
        integer :: base
        contains
            procedure :: next
    end type seq
    
    private
    public  :: ran2, ranu, rang, randint, init_rng, seq

    contains

        real(kind=wp) function next(this) result(res)

            class(seq) :: this

            real(kind=wp) :: fraction
            integer :: i

            fraction = 1.
            res = 0.
            i = this%index

            do while(i > 0)
                fraction = fraction / this%base
                res = res + (fraction * mod(i, this%base))
                i = floor(i / real(this%base, kind=wp))
            end do

            this%index = this%index + 1

        end function next

        subroutine init_rng(input_seed, fwd)
        !! initiate RNG state with reproducible state
            !> input seed
            integer, optional, intent(IN) :: input_seed(:)
            !> boolean that if True runs the generator for 100 steps before returning
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
                ffwd = .false.
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
        !! wrapper for call random number

            real(kind=wp) :: res

            call random_number(res)

        end function ran2

        function ranu(a, b) result(res)
        !! uniformly sample in range[a, b)

            real(kind=wp) :: res
            !> lower bound
            real(kind=wp), intent(IN) :: a
            !> upper bound
            real(kind=wp), intent(IN) :: b

            res = a + ran2() * (b - a)

        end function ranu

        subroutine rang(x, y, avg, sigma)
        !! sample a 2D Guassian distribution

            !> mean of the gaussian to sample from
            real(kind=wp), intent(IN)  :: avg
            !> \(\sigma\) of the guassian to sample from.
            real(kind=wp), intent(IN)  :: sigma
            !> first value to return
            real(kind=wp), intent(OUT) :: x
            !> 2nd value to return
            real(kind=wp), intent(OUT) :: y            
            
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
        !! sample a random integer between [a, b]
            !> lower bound
            integer, intent(IN) :: a
            !> higher bound
            integer, intent(IN) :: b

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