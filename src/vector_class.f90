Module vector_class
! module contains full vector class implmentation
!
!
    use constants, only : wp

    type :: vector
        real(kind=wp)  :: x, y, z
        contains

        procedure :: magnitude         => magnitude_fn
        procedure :: print             => print_sub
        procedure :: length            => length
        generic   :: operator(.dot.)   => vec_dot_vec, vec_dot_mat
        generic   :: operator(.cross.) => vec_cross_vec
        generic   :: operator(/)       => vec_div_scal_r4, vec_div_scal_int
        generic   :: operator(*)       => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)       => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)       => vec_minus_vec, vec_minus_scal, scal_minus_vec

        procedure, pass(a), private :: vec_dot_vec
        procedure, pass(a), private :: vec_dot_mat

        procedure, pass(a), private :: vec_cross_vec

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
    public :: magnitude, vector, print, abs, length, max, nint, min

    interface nint
        module procedure nint_vec
    end interface nint

    interface abs
        module procedure abs_vec
    end interface abs

    interface max
        module procedure max_vec
        module procedure maxval_vec
    end interface max

    interface min
        module procedure min_vec
        module procedure minval_vec
    end interface min

    contains

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

        type(vector) function min_vec(this, val)

            implicit none

            type(vector),  intent(IN) :: this
            real(kind=wp), intent(IN) :: val

            min_vec = vector(min(this%x, val), min(this%y, val), min(this%z, val))

        end function min_vec


        real(kind=wp) function maxval_vec(this)

            implicit none

            type(vector),  intent(IN) :: this

            maxval_vec = max(this%x, this%y, this%z)

        end function maxval_vec

        real(kind=wp) function minval_vec(this)

            implicit none

            type(vector),  intent(IN) :: this

            minval_vec = min(this%x, this%y, this%z)

        end function minval_vec

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

            scal_minus_vec = vector(a - b%x, a - b%y, a - b%z)

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

        function vec_cross_vec(a, b) result(cross)

            class(vector), intent(in) :: a
            type(vector),  intent(in) :: b
            type(vector) :: cross

            cross%x = a%y*b%z - a%z*b%y
            cross%y = a%x*b%z - a%z*b%x
            cross%z = a%x*b%y - a%y*b%x

        end function vec_cross_vec

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

            tmp = this%length()
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
end Module vector_class
