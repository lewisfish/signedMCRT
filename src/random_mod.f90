module random

    use vector_class

    implicit none

    private
    public  :: ran2, ranu, rang, init_rng, hemi

    contains

        subroutine init_rng(input_seed, fwd)
        ! initiate RNG state with reproducable state

            implicit none
            
            integer, optional, intent(IN) :: input_seed(:)
            logical, optional, intent(IN) :: fwd

            integer, allocatable :: seed(:)
            integer              :: n, i
            logical              :: ffwd
            real                 :: a


            call random_seed(size=n)
            allocate(seed(n))

            if(present(input_seed))then
                seed = 0
                seed(:) = input_seed
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


        type(vector) function hemi(cosThetaMax, x, y, z)

            use constants, only: twopi
            use utils, only : lerp

            implicit none

            ! real :: z, phi, r
            type(vector) :: x, y, z
            real :: ran, cosTheta, sinTheta, phi, cosThetaMax

            ! z = ran2()
            ! r = sqrt(max(0., 1.-z**2))
            ! phi = twopi*ran2()

            ! hemi = vector(r * cos(phi), r * sin(phi), z)
            ! ran = ran2()
            ! cosTheta = 1. - ran + ran*cosThetaMax
            ! sinTheta = sqrt(1. - cosTheta**2)
            ! phi = ran2()*twopi
            ! hemi = vector(cos(phi)*sinTheta, sin(phi)*sinTheta, cosTheta)

            cosTheta = lerp(ran2(), cosThetaMax, 1.)
            sinTheta = sqrt(max(0., 1. - cosTheta**2))
            phi = ran2()*twopi
            hemi = cos(phi)*sinTheta*x + sin(phi)*sinTheta*y + cosTheta*z

        end function hemi


end module random