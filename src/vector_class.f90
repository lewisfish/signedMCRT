Module vector_class
! module contains full vector class implmentation
!
!
    use constants, only : wp

    type :: vector
        real(kind=wp)  :: x, y, z
        contains

        procedure :: magnitude       => magnitude_fn
        procedure :: print           => print_sub
        procedure :: length          => length
        generic   :: operator(.dot.) => vec_dot_vec, vec_dot_mat
        generic   :: operator(/)     => vec_div_scal_r4, vec_div_scal_int
        generic   :: operator(*)     => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)     => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)     => vec_minus_vec, vec_minus_scal, scal_minus_vec

        procedure, pass(a), private :: vec_dot_vec
        procedure, pass(a), private :: vec_dot_mat

        procedure, pass(a), private :: vec_div_scal_r4
        procedure, pass(a), private :: vec_div_scal_int

        procedure, pass(a), private :: vec_mult_vec
        procedure, pass(a), private :: vec_mult_scal
        procedure, pass(b), private :: scal_mult_vec

        procedure, pass(a), private :: vec_add_vec
        procedure, pass(a), private :: vec_add_scal
        procedure, pass(b), private :: scal_add_vec

        procedure, pass(a), private :: vec_minus_vec
        procedure, pass(a), private :: vec_minus_scal
        procedure, pass(b), private :: scal_minus_vec

    end type vector

    private
    public :: magnitude, vector, print, abs, length, max, invert, nint, clamp_vec

    interface nint
        module procedure nint_vec
    end interface nint

    interface abs
        module procedure abs_vec
    end interface abs

    interface max
        module procedure max_vec
    end interface max

    contains

        type(vector) function clamp_vec(this, lo, hi)

            use utils, only : clamp

            implicit none

            type(vector), intent(IN) :: this, lo, hi

            real(kind=wp) :: x, y, z

            x = clamp(this%x, lo%x, hi%x)
            y = clamp(this%y, lo%y, hi%y)
            z = clamp(this%z, lo%z, hi%z)

            clamp_vec = vector(x, y, z)

        end function clamp_vec


        type(vector) function nint_vec(this)

            implicit none

            type(vector), intent(IN) :: this

            nint_vec = vector(real(nint(this%x), kind=wp), real(nint(this%y), kind=wp), real(nint(this%z), kind=wp))

        end function nint_vec

        type(vector) function abs_vec(this)

            implicit none

            type(vector), intent(IN) :: this

            abs_vec = vector(abs(this%x), abs(this%y), abs(this%z))

        end function abs_vec

        type(vector) function max_vec(this, val)

            implicit none

            type(vector),  intent(IN) :: this
            real(kind=wp), intent(IN) :: val

            max_vec = vector(max(this%x, val), max(this%y, val), max(this%z, val))

        end function max_vec

        type(vector) function min_fn(this, val)

            implicit none

            type(vector),  intent(IN) :: this
            real(kind=wp), intent(IN) :: val

            min_fn = vector(min(this%x, val), min(this%y, val), min(this%z, val))

        end function min_fn

        type(vector) function vec_minus_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_minus_vec = vector(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus_vec


        type(vector) function vec_add_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_add_scal = vector(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal


        type(vector) function scal_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_add_vec = vector(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vector) function vec_minus_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_minus_scal = vector(a%x - b, a%y - b, a%z - b)

        end function vec_minus_scal


        type(vector) function scal_minus_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_minus_vec = vector(b%x - a, b%y - a, b%z - a)

        end function scal_minus_vec


        type(vector) function vec_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_add_vec = vector(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        elemental function vec_dot_vec(a, b) result (dot)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b
            real(kind=wp) :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vec_dot_vec

        function vec_dot_mat(a, b) result (dot)

            implicit none

            class(vector), intent(IN) :: a
            real(kind=wp), intent(IN) :: b(4, 4)
            type(vector) :: dot

            dot%x = b(1, 1)*a%x + b(2, 1)*a%y + b(3, 1)*a%z + b(4, 1)*1.
            dot%y = b(1, 2)*a%x + b(2, 2)*a%y + b(3, 2)*a%z + b(4, 2)*1.
            dot%z = b(1, 3)*a%x + b(2, 3)*a%y + b(3, 3)*a%z + b(4, 3)*1.

        end function vec_dot_mat

        type(vector) function vec_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_mult_vec = vector(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vector) function vec_mult_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_mult_scal = vector(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vector) function scal_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_mult_vec = vector(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vector) function vec_div_scal_r4(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_div_scal_r4 = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal_r4

        type(vector) function vec_div_scal_int(a, b)

            implicit none

            class(vector), intent(IN) :: a
            integer,       intent(IN) :: b

            vec_div_scal_int = vector(a%x / real(b, kind=wp), a%y / real(b, kind=wp), a%z / real(b, kind=wp))

        end function vec_div_scal_int


        type(vector) function magnitude_fn(this)

            implicit none

            class(vector) :: this

            real(kind=wp) :: tmp

            tmp = sqrt(this%x**2 + this%y**2 + this%z**2)
            magnitude_fn = this / tmp

        end function magnitude_fn


        real(kind=wp) function length(this)

            implicit none

            class(vector) :: this

            length = sqrt(this%x**2 + this%y**2 + this%z**2)

        end function length


        subroutine print_sub(this)

            implicit none

            class(vector) :: this

            print*,this%x, this%y, this%z

        end subroutine

    pure function invert(A) result(B)
        !! from http://fortranwiki.org/fortran/show/Matrix+inversion
        !! Performs a direct calculation of the inverse of a 4??4 matrix.
        real(kind=wp), intent(in) :: A(4,4)   !! Matrix

        real(kind=wp) :: B(4,4)   !! Inverse matrix
        real(kind=wp) :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1./(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
        - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
        + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
        - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end function invert
end Module vector_class