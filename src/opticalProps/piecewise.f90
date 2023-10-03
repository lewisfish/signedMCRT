module piecewiseMod
!! This file contains the piecewise abstract type, for sampling from constants, 1D or 2D arrays. Inspired by [PBRT](https://www.pbr-book.org/) piecewise class.
!! Currently, the following public types are defined:
    
!! - Constant. Used in the case where there is only one value.
!! - 1D. Used in the case where there is a spectrum
!! - 2D. Used in the case where SLM or other image based source types are needed.
    
!! The piecewise type ensures that there is a method (sample) that can be called on all inherited types, e.g
!! call 2Dimage%p%sample(x, y)
!! will return a position (x,y) from where to release a photon.
!! This class can be used to have multi-spectral or single valued wavelength, or used as a 2D image input source i.e SLMs.
!! NOTE: optical properties are not currently adjusted on wavelength change.
    
    use iso_fortran_env, only : int32, int64
    use constants,       only : wp

    implicit none

    !> Abstract spectrum base type. 
    type, abstract :: piecewise
        contains
            !> Deferred procdure. Used to generate a sample from spectrum or get constant value etc.
            procedure(sampleInterface), deferred :: sample
    end type piecewise

    abstract interface
        subroutine sampleInterface(this, x, y, value)
            use constants, only : wp
            import piecewise
            implicit none
            class(piecewise), intent(in) :: this
            real(kind=wp), intent(out) :: x, y
            real(kind=wp), intent(in), optional :: value
        end subroutine sampleInterface
    end interface
    
    !> Spectrum_t type. Used as a container type
    type :: spectrum_t
        class(piecewise), pointer :: p => null()
    end type spectrum_t
    
    !> Constant piecewise type. i.e a piecewise function that does not change value
    type, extends(piecewise) :: constant
        !> The constant value
        real(kind=wp) :: value
        contains
            !> Sampling routine
            procedure :: sample => getValue
    end type constant

    !> 1D piecewise type. Used for the spectral type
    type, extends(piecewise) :: piecewise1D
        !> Input array to sample from. Should be size(n, 2). 1st column is x-axis, 2nd column is y-axis
        real(kind=wp), allocatable :: array(:, :)
        !> cumulative distribution function (CDF) of array.
        real(kind=wp), allocatable :: cdf(:)
        contains
            !> Overloaded sampling function
            procedure :: sample => sample1D
    end type piecewise1D
    
    !> 2D piecewise type. Used for images
    type, extends(piecewise) :: piecewise2D
        !> Height of each cell
        real(kind=wp) :: cell_height
        !> Width of each cell
        real(kind=wp) :: cell_width
        !>cumulative distribution function (CDF) of array.
        real(kind=wp), allocatable :: cdf(:)
        !> Offsets
        integer, private :: xoffset, yoffset
        contains
            !> Overloaded sampling function
            procedure :: sample => sample2D
    end type piecewise2D

    interface piecewise1D
        !> Initalise piecewise1D
        module procedure init_piecewise1D
    end interface piecewise1D

    interface piecewise2D
        !> Initalise piecewise2D
        module procedure init_piecewise2D
    end interface piecewise2D

    ! private
    public :: spectrum_t, piecewise, piecewise1D, piecewise2D, constant

    contains

    subroutine getValue(this, x, y, value)
        !! The constant version of sample
        
        class(constant), intent(in) :: this
        !> Output value
        real(kind=wp),   intent(out) :: x
        !> Not used. Kept to keep interface the same for constant, piecewise1D and piecewise2D
        real(kind=wp),   intent(out) :: y
        !> Not used. Kept to keep interface the same for constant, piecewise1D and piecewise2D
        real(kind=wp), intent(in), optional :: value

        x = this%value
        y = -9999._wp

    end subroutine getValue

    subroutine sample1D(this, x, y, value)
        !! Randomly sample from 1D array
        use random, only : ran2, ranu

        class(piecewise1D), intent(in)  :: this
        !> Return value
        real(kind=wp),      intent(out) :: x
        !> Not used, but here so we can have same interface as 2D sample routine.
        real(kind=wp),      intent(out) :: y
        !> Optional x value. If not present we generate a random one in the range [0., 1.] 
        real(kind=wp), intent(in), optional :: value

        integer(kind=int64) :: idx
        real(kind=wp)       :: val

        if(.not. present(value))then
            !get random x coordinate then get corresponding y
            val = ran2()
            call search_1D(this%cdf, idx, val)
            
            x = this%array(idx, 1) + &
            ((val - this%cdf(idx))*(this%array(idx + 1, 1) - this%array(idx, 1))) / (this%cdf(idx + 1) - this%cdf(idx))

        else
            !already have x so get y
            call search_2D(this%array, idx, value)
            x = this%array(idx, 2) + (this%array(idx+1, 2) - this%array(idx, 2)) * &
                   ((value - this%array(idx, 1))/(this%array(idx+1, 1) - this%array(idx, 1)))
        end if

    end subroutine sample1D


    type(piecewise1D) function init_piecewise1D(array) result(res)
        !! initalise the piecewise1D type with an array size (n, 2). Calculates the CDF of this array.
        !> Input array

        use stdlib_quadrature, only: trapz_weights

        real(kind=wp), intent(in) :: array(:, :)
        
        integer :: i, length
        real(kind=wp) :: weights(size(array, 1)), sumer

        if(size(array, 2) /= 2)error stop "Array must be size (n, 2)"

        res%array = array
        length = size(array, 1)
        allocate(res%cdf(length))
        res%cdf = 0.
        ! Generate CDF array from PDF array via Trapezoidal rule
        weights = trapz_weights(array(:, 1))
        sumer = 0.
        do i = 2, length
            sumer = sumer + weights(i)*array(i,2)
            res%cdf(i) = sumer
        end do
         ! normalise
        res%cdf=res%cdf/res%cdf(length)
    end function init_piecewise1D


    subroutine sample2D(this, x, y, value)
        ! TODO cite where you got this from...
        use random, only : ran2, ranu

        class(piecewise2D), intent(in)  :: this
        real(kind=wp),      intent(out) :: x, y
        real(kind=wp), intent(in), optional :: value

        integer(kind=int32) :: xr, yr
        integer(kind=int64) :: idx
        real(kind=wp)       :: val

        val = ran2()
        call search_1D(this%cdf, idx, val)
        call decode(idx, xr, yr)

        x = real(xr - this%xoffset, kind=wp) + ranu(-this%cell_width, this%cell_width)
        y = real(yr - this%yoffset, kind=wp) + ranu(-this%cell_height, this%cell_height)

    end subroutine sample2D

    
    type(piecewise2D) function init_piecewise2D(cell_width, cell_height, image)
        !! Initalise the piecewise2D type with a given cell_width, cell_height and input image

        !> Input cell width
        real(kind=wp), intent(in) :: cell_width
        !> Input cell height
        real(kind=wp), intent(in) :: cell_height
        !> Input image
        real(kind=wp), intent(in) :: image(:,:)
        
        real(kind=wp), allocatable :: HC1D(:), imagenew(:,:)
        
        integer :: width, height, w2, h2
        integer(kind=int64) :: i
        integer(kind=int32) :: x, y
        
        width = size(image, 1)
        height = size(image, 2)
        
        ! need to pad image for z-order to work...
        w2 = nextpwr2(width)
        h2 = nextpwr2(height)
        
        allocate(imagenew(w2, h2))
        imagenew = 0.
        
        init_piecewise2D%xoffset = (h2 - height)/2
        init_piecewise2D%yoffset = (w2 - width)/2
        
        imagenew(init_piecewise2D%xoffset:init_piecewise2D%xoffset+width-1, &
        init_piecewise2D%yoffset:init_piecewise2D%yoffset+height-1) = image
        
        allocate(init_piecewise2D%cdf(w2 * h2))
        allocate(HC1D(w2 * h2))
        
        HC1D = 0.
        
        do i = 0, (h2*w2)-1
            call decode(i, x, y)
            HC1D(i+1) = imagenew(x+1, y+1)
        end do
        
        init_piecewise2D%cdf(1) = HC1D(1)
        do i = 2, size(HC1D)
            init_piecewise2D%cdf(i) = init_piecewise2D%cdf(i-1) + HC1D(i)
        end do
        
        init_piecewise2D%cell_height = cell_height
        init_piecewise2D%cell_width = cell_width
        init_piecewise2D%cdf = init_piecewise2D%cdf/ init_piecewise2D%cdf(size(init_piecewise2D%cdf))
        
    end function init_piecewise2D
    
    integer function nextpwr2(v) result(res)
    !! Get the next power of 2. i.e given 5 will return 8 (4^2)
    !! only works on 32bit ints
    !! [ref](https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2)
        integer, intent(in) :: v

        res = v - 1
        res = ior(res, rshift(res, 1))
        res = ior(res, rshift(res, 2))
        res = ior(res, rshift(res, 4))
        res = ior(res, rshift(res, 8))
        res = ior(res, rshift(res, 16))
        res = res + 1

    end function nextpwr2

    subroutine search_1D(array, nlow, value)
        !! search by bisection for 1D array
        
        !> Array to search
        real(kind=wp),       intent(in)  :: array(:)
        !> index of found value
        integer(kind=int64), intent(out) :: nlow
        !> value to find in 1D array
        real(kind=wp),       intent(in)  :: value
        
        integer :: nup, middle
        
        nup = size(array)
        nlow = 1
        middle = int((nup+nlow)/2.)

        do while((nup - nlow) > 1)
            middle = int((nup + nlow)/2.)
            if(value > array(middle))then
                nlow = middle
            else
                nup = middle   
            end if
        end do
    end subroutine search_1D

    subroutine search_2D(array, nlow, value)
        !! search by bisection for 1D array
        
        !> 2D array to search. Only searches 1st column
        real(kind=wp),       intent(in)  :: array(:, :)
        !> Index of found index
        integer(kind=int64), intent(out) :: nlow
        !> Value to find in the array.
        real(kind=wp),       intent(in)  :: value
        
        integer :: nup, middle
        
        nup = size(array, 1)
        nlow = 1
        middle = int((nup+nlow)/2.)

        do while((nup - nlow) > 1)
            middle = int((nup + nlow)/2.)
            if(value > array(middle, 1))then
                nlow = middle
            else
                nup = middle   
            end if
        end do
    end subroutine search_2D

    integer(kind=int64) function pack_bits(z) result(x)
    !! Reverse the split function. I.e go from 0a0b0c0d to abcd
    !! Adapted from archer2 cpp [course](https://github.com/EPCCed/archer2-cpp/tree/main/exercises/morton-order)
        !> Input interleaved integer
        integer(kind=int64), intent(in) :: z

        x = z

        x = iand(x, 6148914691236517205_int64)
        x = ior(rshift(x, 1), x)
        x = iand(x, 3689348814741910323_int64)
        x = ior(rshift(x, 2), x)
        x = iand(x, 1085102592571150095_int64)
        x = ior(rshift(x, 4), x)
        x = iand(x, 71777214294589695_int64)
        x = ior(rshift(x, 8), x)
        x = iand(x, 281470681808895_int64)
        x = ior(rshift(x, 16), x)

    end function pack_bits


    subroutine decode(z, x, y)
        !! Compute the 2 indices from a Morton index
        !! Adapted from archer2 cpp [course](https://github.com/EPCCed/archer2-cpp/tree/main/exercises/morton-order)
        
        !> Morton Index
        integer(kind=int64), intent(in) :: z
        !> The computed indices
        integer(kind=int32), intent(out) :: x, y

        integer(kind=int64) :: i, j

        i = z
        x = pack_bits(i)
        j = rshift(z, 1)
        y = pack_bits(j)

    end subroutine decode
end module piecewiseMod