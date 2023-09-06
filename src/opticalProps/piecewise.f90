module piecewiseMod

    use iso_fortran_env, only : int32, int64
    use constants,       only : wp

    implicit none

    type, abstract :: piecewise
        contains
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

    type :: spectrum_t
        class(piecewise), pointer :: p => null()
    end type spectrum_t

    type, extends(piecewise) :: constant
        real(kind=wp) :: value
        contains
            procedure :: sample => getValue
    end type constant

    type, extends(piecewise) :: piecewise1D
        real(kind=wp), allocatable :: array(:, :), cdf(:)
        contains
            procedure :: sample => sample1D
    end type piecewise1D

    type, extends(piecewise) :: piecewise2D
        real(kind=wp) :: cell_width, cell_height
        real(kind=wp), allocatable :: cdf(:)
        integer, private :: xoffset, yoffset
        contains
            procedure :: sample => sample2D
    end type piecewise2D

    interface piecewise1D
        module procedure init_piecewise1D
    end interface piecewise1D

    interface piecewise2D
        module procedure init_piecewise2D
    end interface piecewise2D

    private
    public :: spectrum_t, piecewise, piecewise1D, piecewise2D, constant

    contains

    subroutine getValue(this, x, y, value)

        class(constant), intent(in) :: this
        real(kind=wp),   intent(out) :: x, y
        real(kind=wp), intent(in), optional :: value

        x = this%value
        y = -9999._wp

    end subroutine getValue

    subroutine sample1D(this, x, y, value)
        
        use random, only : ran2, ranu

        class(piecewise1D), intent(in)  :: this
        real(kind=wp),      intent(out) :: x, y
        real(kind=wp), intent(in), optional :: value

        integer(kind=int64) :: idx
        real(kind=wp)       :: val

        if(.not. present(value))then
            !get random x coordinate then get corresponding y
            val = ran2()
            call search_1D(this%cdf, idx, val)
            ! linear interpolation
            ! offset by +1 if not get wrong values. e.g miss out on end values and get more near start.
            x = this%array(idx+1, 1) + (this%array(idx+2, 1) - this%array(idx+1, 1)) * &
                   ((val - this%cdf(idx))/(this%cdf(idx+1) - this%cdf(idx)))
        else
            !already have x so get y
            call search_2D(this%array, idx, value)
            x = this%array(idx, 2) + (this%array(idx+1, 2) - this%array(idx, 2)) * &
                   ((value - this%array(idx, 1))/(this%array(idx+1, 1) - this%array(idx, 1)))
        end if

    end subroutine sample1D


    type(piecewise1D) function init_piecewise1D(array) result(res)

        real(kind=wp), intent(in) :: array(:, :)
        
        integer :: i, j
        real :: summ

        res%array = array
        allocate(res%cdf(size(array,1)-1))

        do j = 1 , size(array, 1)-1
        summ = 0.
            do i = 1, j
                summ = summ + 0.5 * (array(i+1, 2) + array(i, 2)) * (array(i+1, 1) - array(i, 1))
            end do
            res%cdf(j) = summ
        end do
        res%cdf = res%cdf / maxval(res%cdf)

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
    
        real(kind=wp), intent(in) :: cell_width, cell_height
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
    ! https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
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
        !
        !  search by bisection for 1D array
        !
        
        integer(kind=int64), intent(out) :: nlow
        real(kind=wp),       intent(in)  :: array(:)
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
        !
        !  search by bisection for 1D array
        !
        
        integer(kind=int64), intent(out) :: nlow
        real(kind=wp),       intent(in)  :: array(:, :)
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



    integer(kind=int64) function split(a) result(x)
    ! taken from archer2 cpp course
        integer(kind=int32) :: a

        x = a
        x = iand(ior(x, lshift(x, 16)), 281470681808895_int64)
        x = iand(ior(x, lshift(x, 8)), 71777214294589695_int64)
        x = iand(ior(x, lshift(x, 4)), 1085102592571150095_int64)
        x = iand(ior(x, lshift(x, 2)), 3689348814741910323_int64)
        x = iand(ior(x, lshift(x, 1)), 6148914691236517205_int64)

    end function split

    integer(kind=int64) function pack_bits(z) result(x)
    ! taken from archer2 cpp course

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

    integer(kind=int64) function encode(x, y)
    ! taken from archer2 cpp course

        integer(kind=int32), intent(in) :: x, y
        
        encode = ior(split(x), lshift(split(y),1))
        
    end function encode

    subroutine decode(z, x, y)
    ! taken from archer2 cpp course

        integer(kind=int64), intent(in) :: z
        integer(kind=int32), intent(out) :: x, y

        integer(kind=int64) :: i, j

        i = z
        x = pack_bits(i)
        j = rshift(z, 1)
        y = pack_bits(j)

    end subroutine decode
end module piecewiseMod
! program test

!     use, intrinsic :: iso_c_binding
!     use iso_fortran_env, only : int64

!     use piecewiseMod
!     use random,    only : rang, init_rng
!     use constants, only : wp

!     implicit none

!     integer ::  ix, iy, u
!     integer, parameter :: n=129
!     real(kind=wp) :: xx, yy, bin_wid
!     integer(kind=int64) :: i
!     real(kind=wp) :: xr, yr
!     real :: data2d(n, n), data1d(n, 2)
!     type(piecewise1D) :: obj
    
!     !set seed
!     call init_rng(spread(123456789, 1, 8), .true.)
    
!     data1d = 0.
!     bin_wid = 20. / n

!     do i = 1, n
!         xr = i*bin_wid
!         yr = exp(-(xr-10.)**2 / (2*(0.9)**2))
!         data1d(i, 1) = xr
!         data1d(i, 2) = yr
!     end do

!     open(newunit=u, file="../res/spectrum.dat")
!     do i = 1, n
!         write(u,"(es24.16e3,1x,es24.16e3)") data1d(i, 1),  data1d(i,2)
!     end do
!     close(u)
    
    
!     data1d(:,2) = data1d(:,2) / maxval(data1d,2)
!     obj = piecewise1D(data1d)

!     open(newunit=u,file="cdf.dat")
!     do i = 1, size(obj%cdf)
!         write(u,*) obj%cdf(i)
!     end do
!     close(u)

!     open(newunit=u, file="sampled.dat")
!     do i = 1, 1000
!         call obj%sample(xr, yr)
!         write(u,*) xr
!     end do
!     close(u)


    !generate image
    ! data2d = 0.
    ! data1d = 0.
    ! bin_wid = 2. / n
    ! do i = 1, 100000
    !     call rang(xx, yy, 1._wp, 0.1_wp)
    !     ix = nint(xx / bin_wid)-10
    !     iy = nint(yy / bin_wid)-30
    !     if(ix > n .or. ix < 1)cycle
    !     if(iy > n .or. iy < 1)cycle
    !     data2d(ix, iy) = data2d(ix, iy) + 1.
    !     call rang(xx, yy, 1._wp, 0.1_wp)
    !     ix = nint(xx / bin_wid)+30
    !     iy = nint(yy / bin_wid)+30
    !     if(ix > n .or. ix < 1)cycle
    !     if(iy > n .or. iy < 1)cycle
    !     data2d(ix, iy) = data2d(ix, iy) + 1.
    ! end do

    ! open(newunit=u, file="image.dat")
    ! do i = 1, n
    !     write(u,*) data2d(:, i)
    ! end do
    ! close(u)

    ! data2d = data2d / maxval(data2d)

    ! obj = piecewise2D(0.5_wp, 0.5_wp, data2d)

    ! open(newunit=u,file="cdf.dat")
    ! do i = 1, (n**2)
    !     write(u,*) obj%cdf(i)
    ! end do
    ! close(u)

    ! open(newunit=u, file="sampled.dat")
    ! do i = 1, 1000
    !     call obj%sample(xr, yr)
    !     write(u,*) xr, yr
    ! end do
    ! close(u)
! end program