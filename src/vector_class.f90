Module vector_class

    type :: vector
        real :: x, y, z
        contains

        procedure :: magnitude       => magnitude_fn
        procedure :: print           => print_sub
        procedure :: length          => length
        generic   :: operator(.dot.) => vec_dot
        generic   :: operator(/)     => vec_div_scal
        generic   :: operator(*)     => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)     => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)     => vec_minus_vec

        procedure, pass(a), private :: vec_dot

        procedure, pass(a), private :: vec_div_scal

        procedure, pass(a), private :: vec_mult_vec
        procedure, pass(a), private :: vec_mult_scal
        procedure, pass(b), private :: scal_mult_vec

        procedure, pass(a), private :: vec_add_vec
        procedure, pass(a), private :: vec_add_scal
        procedure, pass(b), private :: scal_add_vec

        procedure, pass(a), private :: vec_minus_vec

    end type vector

    private
    public :: magnitude, vector, print, abs, length, max


    interface abs
        module procedure abs_vec
    end interface abs

    interface max
        module procedure max_vec
    end interface max

    contains

        type(vector) function abs_vec(this)

            implicit none

            type(vector), intent(IN) :: this

            abs_vec = vector(abs(this%x), abs(this%y), abs(this%z))

        end function abs_vec

        type(vector) function max_vec(this, val)

            implicit none

            type(vector), intent(IN) :: this
            real, intent(IN) :: val

            max_vec = vector(max(this%x, val), max(this%y, val), max(this%z, val))

        end function max_vec

        type(vector) function min_fn(this, val)

            implicit none

            type(vector), intent(IN) :: this
            real, intent(IN) :: val

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
            real,          intent(IN) :: b

            vec_add_scal = vector(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal


        type(vector) function scal_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real,          intent(IN) :: a

            scal_add_vec = vector(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vector) function vec_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_add_vec = vector(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        elemental function vec_dot(a, b) result (dot)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b
            real :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vec_dot


        type(vector) function vec_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_mult_vec = vector(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vector) function vec_mult_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real,          intent(IN) :: b

            vec_mult_scal = vector(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vector) function scal_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real,          intent(IN) :: a

            scal_mult_vec = vector(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vector) function vec_div_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vec_div_scal = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal


        type(vector) function magnitude_fn(this)

            implicit none

            class(vector) :: this

            real :: tmp

            tmp = sqrt(this%x**2 + this%y**2 + this%z**2)
            magnitude_fn = this / tmp

        end function magnitude_fn


        real function length(this)

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