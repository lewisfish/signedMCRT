Module vec4_class

    use constants, only : wp

    implicit none

    !not fully implmented vec4 class
    type :: vec4
        real(kind=wp) :: x, y, z, p
        contains

        generic   :: operator(.dot.) => vec_dot_vec
        generic   :: operator(/)     => vec_div_scal_r4, vec_div_scal_int
        generic   :: operator(*)     => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)     => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)     => vec_minus_vec, vec_minus_scal, scal_minus_vec

        procedure, pass(a), private :: vec_dot_vec

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

    end type vec4

    interface sin
        module procedure sin_vec
    end interface sin

    
    interface vec4
        module procedure init_vec4_vector_real
    end interface vec4


    private
    public :: vec4, sin
    contains
        type(vec4) function init_vec4_vector_real(vec, val) result(out)

            use vector_class

            type(vector),  intent(in) :: vec
            real(kind=wp), intent(in) :: val

            out%x = vec%x
            out%y = vec%y
            out%z = vec%z
            out%p = val

        end function init_vec4_vector_real

        type(vec4) pure elemental function sin_vec(p)

            type(vec4), intent(IN) :: p


            sin_vec = vec4(sin(p%x), sin(p%y), sin(p%z), sin(p%p))

        end function sin_vec

        type(vec4) function vec_minus_vec(a, b)

            class(vec4), intent(IN) :: a
            type(vec4),  intent(IN) :: b

            vec_minus_vec = vec4(a%x - b%x, a%y - b%y, a%z - b%z, a%p - b%p)

        end function vec_minus_vec

        type(vec4) function vec_add_scal(a, b)

            class(vec4),     intent(IN) :: a
            real(kind=wp),   intent(IN) :: b

            vec_add_scal = vec4(a%x + b, a%y + b, a%z + b, a%p + b)

        end function vec_add_scal

        type(vec4) function scal_add_vec(a, b)

            class(vec4),   intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_add_vec = vec4(b%x + a, b%y + a, b%z + a,  b%p + a)

        end function scal_add_vec

        type(vec4) function vec_minus_scal(a, b)

            class(vec4),   intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_minus_scal = vec4(a%x - b, a%y - b, a%z - b, a%p - b)

        end function vec_minus_scal

        type(vec4) function scal_minus_vec(a, b)

            class(vec4),   intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_minus_vec = vec4(b%x - a, b%y - a, b%z - a, b%p - a)

        end function scal_minus_vec

        type(vec4) function vec_add_vec(a, b)

            class(vec4), intent(IN) :: a
            type(vec4),  intent(IN) :: b

            vec_add_vec = vec4(a%x + b%x, a%y + b%y, a%z + b%z, a%p + b%p)

        end function vec_add_vec

        elemental function vec_dot_vec(a, b) result (dot)

            class(vec4), intent(IN) :: a
            type(vec4),  intent(IN) :: b
            real(kind=wp) :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z) + (a%p * b%p)

        end function vec_dot_vec

        type(vec4) function vec_mult_vec(a, b)

            class(vec4), intent(IN) :: a
            type(vec4),  intent(IN) :: b

            vec_mult_vec = vec4(a%x * b%x, a%y * b%y, a%z * b%z, a%p * b%p)

        end function vec_mult_vec

        type(vec4) function vec_mult_scal(a, b)

            class(vec4),   intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_mult_scal = vec4(a%x * b, a%y * b, a%z * b, a%p * b)

        end function vec_mult_scal

        type(vec4) function scal_mult_vec(a, b)

            class(vec4),   intent(IN) :: b
            real(kind=wp), intent(IN) :: a

            scal_mult_vec = vec4(a * b%x, a * b%y, a * b%z, a * b%p)

        end function scal_mult_vec

        type(vec4) function vec_div_scal_r4(a, b)

            class(vec4),   intent(IN) :: a
            real(kind=wp), intent(IN) :: b

            vec_div_scal_r4 = vec4(a%x / b, a%y / b, a%z / b, a%p / b)

        end function vec_div_scal_r4

        type(vec4) function vec_div_scal_int(a, b)

            class(vec4), intent(IN) :: a
            integer,       intent(IN) :: b

            vec_div_scal_int = vec4(a%x / real(b, kind=wp), a%y / real(b, kind=wp), a%z / real(b, kind=wp), a%p / real(b, kind=wp))

        end function vec_div_scal_int

        type(vec4) function magnitude_fn(this)

            class(vec4) :: this

            real(kind=wp) :: tmp

            tmp = sqrt(this%x**2 + this%y**2 + this%z**2 + this%p**2)
            magnitude_fn = this / tmp

        end function magnitude_fn

        real(kind=wp) function length(this)

            class(vec4) :: this

            length = sqrt(this%x**2 + this%y**2 + this%z**2 + this%p**2)

        end function length
end Module vec4_class