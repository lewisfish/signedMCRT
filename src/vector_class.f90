module vector_class
!! Vector class module. Defines a vector type (x, y, z) and associated operations on vectors and other types.
  
    use constants, only : wp

    implicit none

    !> Vector class
    type :: vector
        !> vector components
        real(kind=wp)  :: x, y, z
        contains

        procedure :: magnitude         => magnitude
        procedure :: print             => print
        procedure :: length            => length
        !> .dot. operator. Dot product
        generic   :: operator(.dot.)   => vec_dot_vec, vec_dot_mat
        !> .cross. operator. Cross product
        generic   :: operator(.cross.) => vec_cross_vec
        !> Overloads the Division operator for vec3
        generic   :: operator(/)       => vec_div_scal_r4, vec_div_scal_r8, vec_div_scal_int
        !> Overloads the Multiplication operator for vec3
        generic   :: operator(*)       => vec_mult_vec, vec_mult_scal, scal_mult_vec
        !> Overloads the exponential operator for vec3
        generic   :: operator(**)      => vec_mult_exp_scal_int, vec_mult_exp_scal_r4, vec_mult_exp_scal_r8
        !> Overloads the Addition operator for vec3
        generic   :: operator(+)       => vec_add_vec, vec_add_scal, scal_add_vec
        !> Overloads the Subtraction operator for vec3
        generic   :: operator(-)       => vec_minus_vec, vec_minus_scal, scal_minus_vec
        !> Overloads the equal operator for vec3
        generic   :: operator(==)      => vec_equal_vec

        procedure, pass(a), private :: vec_dot_vec
        procedure, pass(a), private :: vec_dot_mat

        procedure, pass(a), private :: vec_cross_vec

        procedure, pass(a), private :: vec_div_scal_r4
        procedure, pass(a), private :: vec_div_scal_r8
        procedure, pass(a), private :: vec_div_scal_int

        procedure, pass(a), private :: vec_mult_vec
        procedure, pass(a), private :: vec_mult_scal
        procedure, pass(b), private :: scal_mult_vec
        
        procedure, pass(a), private :: vec_mult_exp_scal_int
        procedure, pass(a), private :: vec_mult_exp_scal_r4
        procedure, pass(a), private :: vec_mult_exp_scal_r8

        procedure, pass(a), private :: vec_add_vec
        procedure, pass(a), private :: vec_add_scal
        procedure, pass(b), private :: scal_add_vec

        procedure, pass(a), private :: vec_minus_vec
        procedure, pass(a), private :: vec_minus_scal
        procedure, pass(b), private :: scal_minus_vec

        procedure, pass(a), private :: vec_equal_vec

    end type vector

    private
    public :: magnitude, vector, print, abs, length, max, nint, min

    interface nint
        !! Overload of the nint intrinsic for a vec3
        module procedure nint_vec
    end interface nint

    interface abs
        !! Overload of the abs intrinsic for a vec3
        module procedure abs_vec
    end interface abs

    interface max
        !! Overload of the max intrinsic for a vec3
        module procedure max_vec
        module procedure maxval_vec
    end interface max

    interface min
        !! Overload of the min intrinsic for a vec3
        module procedure min_vec
        module procedure minval_vec
    end interface min

    contains

        type(vector) pure elemental function vec_mult_exp_scal_int(a, b)
            !! vec3**scalar for integer scalar
            !> Input Vector
            class(vector), intent(in) :: a
            !> Input scalar
            integer,       intent(in) :: b

            vec_mult_exp_scal_int = vector(a%x**b, a%y**b, a%z**b)

        end function vec_mult_exp_scal_int

        type(vector) pure elemental function vec_mult_exp_scal_r4(a, b)
            !! vec3**scalar for 32-bit float scalar

            use constants, only : sp
            
            !> Input Vector
            class(vector), intent(in) :: a
            !> Input scalar
            real(kind=sp), intent(in) :: b

            vec_mult_exp_scal_r4 = vector(a%x**b, a%y**b, a%z**b)

        end function vec_mult_exp_scal_r4

        type(vector) pure elemental function vec_mult_exp_scal_r8(a, b)
            !! vec3**scalar for 64-bit float scalar

            use constants, only : dp
                
            !> Input Vector
            class(vector), intent(in) :: a
            !> Input scalar
            real(kind=dp), intent(in) :: b

            vec_mult_exp_scal_r8 = vector(a%x**b, a%y**b, a%z**b)

        end function vec_mult_exp_scal_r8


        logical pure elemental function vec_equal_vec(a, b)
            !! vec3 == vec3
            !> Input vec3s
            class(vector), intent(in) :: a, b

            vec_equal_vec = .false.

            if(a%x == b%x)then
                if(a%y == b%y)then
                    if(a%z == b%z)then
                        vec_equal_vec = .true.
                    end if
                end if
            end if           

        end function vec_equal_vec

        type(vector) pure elemental function nint_vec(this)
            !! Overload the nint intrinsic for a vec3 elementwise
            !> Input vector
            type(vector), intent(IN) :: this

            nint_vec = vector(real(nint(this%x), kind=wp), real(nint(this%y), kind=wp), real(nint(this%z), kind=wp))

        end function nint_vec

        type(vector) pure elemental function abs_vec(this)
        !! Calculate the absoulte of a vector elementwise
            !> Input vector
            type(vector), intent(IN) :: this

            abs_vec = vector(abs(this%x), abs(this%y), abs(this%z))

        end function abs_vec

        type(vector) pure elemental function max_vec(this, val)
            !! Get the max value elementwise between a vec3 and a scalar
            !> Input vector
            type(vector),  intent(IN) :: this
            !> Input max value
            real(kind=wp), intent(IN) :: val

            max_vec = vector(max(this%x, val), max(this%y, val), max(this%z, val))

        end function max_vec

        type(vector) pure elemental function min_vec(this, val)
            !! Get the min value elementwise between a vec3 and a scalar
            !> Input vector
            type(vector),  intent(IN) :: this
            !> Input minimum value
            real(kind=wp), intent(IN) :: val

            min_vec = vector(min(this%x, val), min(this%y, val), min(this%z, val))

        end function min_vec


        real(kind=wp) pure elemental function maxval_vec(this)
            !! Get the max value in a vec3
            !> Input vector
            type(vector),  intent(IN) :: this

            maxval_vec = max(this%x, this%y, this%z)

        end function maxval_vec

        real(kind=wp) pure elemental function minval_vec(this)
            !! Get the min value in a vec3
            !> Input vector
            type(vector),  intent(IN) :: this

            minval_vec = min(this%x, this%y, this%z)

        end function minval_vec

        type(vector) pure elemental function vec_minus_vec(a, b)
            !! vec3 - vec3
            !> Input vector
            class(vector), intent(IN) :: a
            !> vec3 to subtract
            type(vector),  intent(IN) :: b

            vec_minus_vec = vector(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus_vec


        type(vector) pure elemental function vec_add_scal(a, b)
            !! vec3 + scalar
            !> Input vector
            class(vector), intent(IN) :: a
            !> Scalar to add
            real(kind=wp), intent(IN) :: b

            vec_add_scal = vector(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal


        type(vector) pure elemental function scal_add_vec(a, b)
            !! vec3 + scalar
            !> Input vector
            class(vector), intent(IN) :: b
            !> Scalar to add
            real(kind=wp), intent(IN) :: a

            scal_add_vec = vector(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vector) pure elemental function vec_minus_scal(a, b)
            !! vec3 - scalar
            !> Input vector
            class(vector), intent(IN) :: a
            !> Scalar to subtract
            real(kind=wp), intent(IN) :: b

            vec_minus_scal = vector(a%x - b, a%y - b, a%z - b)

        end function vec_minus_scal


        type(vector) pure elemental function scal_minus_vec(a, b)
            !! scalar - vec3
            !> Input vector
            class(vector), intent(IN) :: b
            !> Scalar to subtract from
            real(kind=wp), intent(IN) :: a

            scal_minus_vec = vector(a - b%x, a - b%y, a - b%z)

        end function scal_minus_vec


        type(vector) pure elemental function vec_add_vec(a, b)
            !! vec3 + vec3
            !> Input vector
            class(vector), intent(IN) :: a
            !> Vec3 to add
            type(vector),  intent(IN) :: b

            vec_add_vec = vector(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        pure elemental function vec_dot_vec(a, b) result (dot)
            !! vec3 . vec3
            !> Input vec3
            class(vector), intent(IN) :: a
            !> vec3 to dot
            type(vector),  intent(IN) :: b
            
            real(kind=wp) :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vec_dot_vec

        pure function vec_dot_mat(a, b) result (dot)
            !! vec3 . matrix
            !> Input vec3
            class(vector), intent(IN) :: a
            !> Matrix to dot with
            real(kind=wp), intent(IN) :: b(4, 4)
            type(vector) :: dot

            dot%x = b(1, 1)*a%x + b(2, 1)*a%y + b(3, 1)*a%z + b(4, 1)*1.
            dot%y = b(1, 2)*a%x + b(2, 2)*a%y + b(3, 2)*a%z + b(4, 2)*1.
            dot%z = b(1, 3)*a%x + b(2, 3)*a%y + b(3, 3)*a%z + b(4, 3)*1.

        end function vec_dot_mat

        pure elemental function vec_cross_vec(a, b) result(cross)
            !! vec3 x vec3
            !> Input vector
            class(vector), intent(in) :: a
            !> vec3 to cross with
            type(vector),  intent(in) :: b
            type(vector) :: cross

            cross%x = a%y*b%z - a%z*b%y
            cross%y = -a%x*b%z + a%z*b%x
            cross%z = a%x*b%y - a%y*b%x

        end function vec_cross_vec

        type(vector) pure elemental function vec_mult_vec(a, b)
            !! vec3 * vec3 elementwise
            !> input vec3
            class(vector), intent(IN) :: a
            !> vec3 to multiply by
            type(vector),  intent(IN) :: b

            vec_mult_vec = vector(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vector) pure elemental function vec_mult_scal(a, b)
            !! vec3 * scalar elementwise
            !> input vec3
            class(vector), intent(IN) :: a
            !> Scalar to multiply by
            real(kind=wp), intent(IN) :: b

            vec_mult_scal = vector(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vector) pure elemental function scal_mult_vec(a, b)
            !! Scalar * vec3 elementwise
            !> input vec3
            class(vector), intent(IN) :: b
            !> Scalar to multiply by
            real(kind=wp), intent(IN) :: a

            scal_mult_vec = vector(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vector) pure elemental function vec_div_scal_r4(a, b)
            !! vec3 / scalar elementwise. Scalar is a 32-bit float
            use constants, only : sp
            !> input vec3
            class(vector), intent(IN) :: a
            !> Scalar to divide by
            real(kind=sp), intent(IN) :: b

            vec_div_scal_r4 = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal_r4

        type(vector) pure elemental function vec_div_scal_r8(a, b)
            !! vec3 / scalar elementwise. Scalar is a 64-bit float
            use constants, only : dp
            !> input vec3
            class(vector), intent(IN) :: a
            !> Scalar to divide by
            real(kind=dp), intent(IN) :: b

            vec_div_scal_r8 = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal_r8

        type(vector) pure elemental function vec_div_scal_int(a, b)
            !! vec3 / scalar elementwise. Scalar is an integer
            !> input vec3
            class(vector), intent(IN) :: a
            !> Scalar to divide by
            integer,       intent(IN) :: b

            vec_div_scal_int = vector(a%x / real(b, kind=wp), a%y / real(b, kind=wp), a%z / real(b, kind=wp))

        end function vec_div_scal_int


        type(vector) pure elemental function magnitude(this)
            !! Returns the magnitude of a vec3

            class(vector), intent(in) :: this

            real(kind=wp) :: tmp

            tmp = this%length()
            magnitude = this / tmp

        end function magnitude


        real(kind=wp) pure elemental function length(this)
            !! Returns the length of a vec3
            class(vector), intent(in) :: this

            length = sqrt(this%x**2 + this%y**2 + this%z**2)

        end function length

        subroutine print(this)

            class(vector) :: this

            print*,this%x, this%y, this%z

        end subroutine print
end Module vector_class
