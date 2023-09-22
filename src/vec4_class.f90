Module vec4_class
!! Vector4 class module. Defines a vector4 type (x, y, z, p) and associated operations on vectors and other types.
    use constants, only : wp

    implicit none

    !> not fully implmented vec4 class
    type :: vec4
        !> vec4 components
        real(kind=wp) :: x, y, z, p
        contains

        procedure :: magnitude         => magnitude_fn
        procedure :: length            => length
        !> .dot. operator
        generic   :: operator(.dot.) => vec_dot_vec
        !> Overloaded Division operator
        generic   :: operator(/)     => vec_div_scal_r4, vec_div_scal_r8, vec_div_scal_int
        !> Overloaded Mulitiplication operator
        generic   :: operator(*)     => vec_mult_vec, vec_mult_scal, scal_mult_vec
        !> Overloaded Addition operator
        generic   :: operator(+)     => vec_add_vec, vec_add_scal, scal_add_vec
        !> Overloaded Subtraction operator
        generic   :: operator(-)     => vec_minus_vec, vec_minus_scal, scal_minus_vec

        procedure, pass(a), private :: vec_dot_vec

        procedure, pass(a), private :: vec_div_scal_r4
        procedure, pass(a), private :: vec_div_scal_r8
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

    end type vec4

    interface sin
        !! Vec4 overload of the sin intrinsic 
        module procedure sin_vec
    end interface sin

    interface vec4
        !! Initalise a vec4 from a vec3 and a scalar
        module procedure init_vec4_vector_real
    end interface vec4

    private
    public :: vec4, sin

    contains
        type(vec4) function init_vec4_vector_real(vec, val) result(out)
            !! Initalise vec4 from a vec3 and Scalar
            !! e.g vec4 = [vec3%x, vec3%y, vec3%z, scalar]
            use vector_class
            !> Input vec3
            type(vector),  intent(in) :: vec
            !> Input Scalar
            real(kind=wp), intent(in) :: val

            out%x = vec%x
            out%y = vec%y
            out%z = vec%z
            out%p = val

        end function init_vec4_vector_real

        type(vec4) pure elemental function sin_vec(p)
            !! Sine of a vec4, elementwise
            !> Input vec4
            type(vec4), intent(IN) :: p

            sin_vec = vec4(sin(p%x), sin(p%y), sin(p%z), sin(p%p))

        end function sin_vec

        type(vec4) pure elemental function vec_minus_vec(a, b)
            !! Elementwise vec4 - vec4
            !> Input vec4
            class(vec4), intent(IN) :: a
            !> vec4 to subtract
            type(vec4),  intent(IN) :: b

            vec_minus_vec = vec4(a%x - b%x, a%y - b%y, a%z - b%z, a%p - b%p)

        end function vec_minus_vec

        type(vec4) pure elemental function vec_add_scal(a, b)
            !! Elementwise vec4 + scalar
            !> Input vec4
            class(vec4),     intent(IN) :: a
            !> Scalar to add
            real(kind=wp),   intent(IN) :: b

            vec_add_scal = vec4(a%x + b, a%y + b, a%z + b, a%p + b)

        end function vec_add_scal

        type(vec4) pure elemental function scal_add_vec(a, b)
            !! Elementwise scalar + vec4
            !> Input vec4
            class(vec4),   intent(IN) :: b
            !> Scalar to add
            real(kind=wp), intent(IN) :: a

            scal_add_vec = vec4(b%x + a, b%y + a, b%z + a,  b%p + a)

        end function scal_add_vec

        type(vec4) pure elemental function vec_minus_scal(a, b)
            !! Elementwise vec4 - scalar
            !> Input vec4
            class(vec4),   intent(IN) :: a
            !> Scalar to subtract
            real(kind=wp), intent(IN) :: b

            vec_minus_scal = vec4(a%x - b, a%y - b, a%z - b, a%p - b)

        end function vec_minus_scal

        type(vec4) pure elemental function scal_minus_vec(a, b)
            !! Elementwise Scalar - vec4
            !> Input vec4
            class(vec4),   intent(IN) :: b
            !> Scalar to subtract
            real(kind=wp), intent(IN) :: a

            scal_minus_vec = vec4(a - b%x, a - b%y, a - b%z, a - b%p)

        end function scal_minus_vec

        type(vec4) pure elemental function vec_add_vec(a, b)
            !! Elementwise vec4 + vec4
            !> Input vec4
            class(vec4), intent(IN) :: a
            !> vec4 to add
            type(vec4),  intent(IN) :: b

            vec_add_vec = vec4(a%x + b%x, a%y + b%y, a%z + b%z, a%p + b%p)

        end function vec_add_vec

        pure elemental function vec_dot_vec(a, b) result (dot)
            !! dot product between two vec4s
            !> Input vec4
            class(vec4), intent(IN) :: a
            !> vec4 to dot with
            type(vec4),  intent(IN) :: b

            real(kind=wp) :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z) + (a%p * b%p)

        end function vec_dot_vec

        type(vec4) pure elemental function vec_mult_vec(a, b)
            !! Elementwise vec4 * vec4
            !> Input vec4
            class(vec4), intent(IN) :: a
            !> vec4 to multiply by
            type(vec4),  intent(IN) :: b

            vec_mult_vec = vec4(a%x * b%x, a%y * b%y, a%z * b%z, a%p * b%p)

        end function vec_mult_vec

        type(vec4) pure elemental function vec_mult_scal(a, b)
            !! Elementwise vec4 * Scalar
            !> Input vec4
            class(vec4),   intent(IN) :: a
            !> Scalar to multiply by
            real(kind=wp), intent(IN) :: b

            vec_mult_scal = vec4(a%x * b, a%y * b, a%z * b, a%p * b)

        end function vec_mult_scal

        type(vec4) pure elemental function scal_mult_vec(a, b)
            !! Elementwise Scalar * vec4
            !> Input vec4
            class(vec4),   intent(IN) :: b
            !> Scalar to multiply by
            real(kind=wp), intent(IN) :: a

            scal_mult_vec = vec4(a * b%x, a * b%y, a * b%z, a * b%p)

        end function scal_mult_vec

        type(vec4) pure elemental function vec_div_scal_r4(a, b)
            !! Elementwise vec4 / Scalar. Scalar is 32-bit float
            use constants, only : sp
        
            !> Input vec4
            class(vec4),   intent(IN) :: a
            !> Scalar to divide by
            real(kind=sp), intent(IN) :: b

            vec_div_scal_r4 = vec4(a%x / b, a%y / b, a%z / b, a%p / b)

        end function vec_div_scal_r4

        type(vec4) pure elemental function vec_div_scal_r8(a, b)
            !! Elementwise vec4 / Scalar. Scalar is 32-bit float
            use constants, only : dp
        
            !> Input vec4
            class(vec4),   intent(IN) :: a
            !> Scalar to divide by
            real(kind=dp), intent(IN) :: b

            vec_div_scal_r8 = vec4(a%x / b, a%y / b, a%z / b, a%p / b)

        end function vec_div_scal_r8

        type(vec4) pure elemental function vec_div_scal_int(a, b)
            !! Elementwise vec4 / Scalar. Scalar is an integer

            !> Input vec4
            class(vec4),   intent(IN) :: a
            !> Scalar to divide by
            integer,     intent(IN) :: b

            vec_div_scal_int = vec4(a%x / real(b, kind=wp), a%y / real(b, kind=wp), a%z / real(b, kind=wp), a%p / real(b, kind=wp))

        end function vec_div_scal_int

        type(vec4) pure elemental function magnitude_fn(this)
            !! Returns the magnitude of a vec4
            class(vec4), intent(in) :: this

            magnitude_fn = this / this%length()

        end function magnitude_fn

        real(kind=wp) pure elemental function length(this)
            !! Returns the length of a vec4
            class(vec4), intent(in) :: this

            length = sqrt(this%x**2 + this%y**2 + this%z**2 + this%p**2)

        end function length
end Module vec4_class