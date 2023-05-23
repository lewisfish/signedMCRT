module sdfModifiers
    
    use constants,   only : wp
    use sdf_baseMod, only : sdf_base, primitive
    use sdfHelpers,  only : identity
    use vector_class

    implicit none
    
    type, extends(sdf_base) :: revolution
        real(kind=wp) :: o
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_revolution
    end type revolution

    type, extends(sdf_base) :: extrude
        real(kind=wp) :: h
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_extrude
    end type extrude

    type, extends(sdf_base) :: onion
        real(kind=wp) :: thickness
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_onion
    end type onion

    type, extends(sdf_base) :: twist
        real(kind=wp) :: k
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_twist
    end type twist

    type, extends(sdf_base) :: displacement
        procedure(primitive), nopass, pointer :: func
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_disp
    end type displacement

    type, extends(sdf_base) :: bend
        real(kind=wp) :: k
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_bend
    end type bend

    type, extends(sdf_base) :: elongate
        type(vector) :: size
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_elongate
    end type elongate

    type, extends(sdf_base) :: repeat
        real(kind=wp) :: c
        type(vector) :: la, lb
        class(sdf_base), pointer :: prim
        contains
        procedure :: evaluate => eval_repeat
    end type repeat

    interface revolution
        module procedure revolution_init
    end interface revolution

    interface extrude
        module procedure extrude_init
    end interface extrude

    interface onion
        module procedure onion_init
    end interface onion

    interface twist
        module procedure twist_init
    end interface twist

    interface displacement
        module procedure displacement_init
    end interface displacement

    interface bend
        module procedure bend_init
    end interface bend

    interface elongate
        module procedure elongate_init
    end interface elongate

    interface repeat
        module procedure repeat_init
    end interface repeat

    private
    public :: onion, extrude, twist, displacement, bend, elongate, repeat, revolution
    public :: union, SmoothUnion, intersection, subtraction

contains
    
    type(twist) function twist_init(prim, k) result(out)

        class(sdf_base), target :: prim
        real :: k

        out%k = k
        out%prim => prim
        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function twist_init

    type(extrude) function extrude_init(prim, h) result(out)

        class(sdf_base), target :: prim
        real(kind=wp), intent(IN)   :: h

        out%h = h
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function extrude_init

    type(elongate) function elongate_init(prim, size) result(out)

        class(sdf_base), target :: prim
        type(vector), intent(IN) :: size

        out%size = size
        out%prim => prim

        out%optProps = prim%optProps
        out%layer = prim%layer
        out%transform = identity()

    end function elongate_init

    type(displacement) function displacement_init(prim, func) result(out)

        class(sdf_base), target :: prim
        procedure(primitive) :: func

        out%func => func
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function displacement_init

    type(bend) function bend_init(prim, k) result(out)

        class(sdf_base), target :: prim
        real(kind=wp), intent(IN) :: k

        out%k = k
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function bend_init

    type(repeat) function repeat_init(prim, c, la, lb) result(out)

        class(sdf_base), target :: prim
        type(vector),  intent(IN) :: la, lb
        real(kind=wp), intent(IN) :: c

        out%c = c
        out%la = la
        out%lb = lb
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function repeat_init

    pure elemental function eval_extrude(this, pos) result(res)

        class(extrude), intent(in) :: this
        type(vector),   intent(IN) :: pos
        real(kind=wp) :: res

        type(vector)  :: w
        real(kind=wp) :: d

        d = this%prim%evaluate(pos)
        w = vector(d, abs(pos%z) - this%h, 0._wp)
        res = min(max(w%x, w%y), 0._wp) + length(max(w, 0._wp))

    end function eval_extrude

    type(revolution) function revolution_init(prim, o) result(out)

        class(sdf_base), target :: prim
        real(kind=wp), intent(IN) :: o

        out%o = o
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function revolution_init

    pure elemental function eval_revolution(this, pos) result(res)

        class(revolution), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: pxz, q

        pxz = vector(pos%x, pos%z, 0._wp)

        q = vector(length(pxz) - this%o, pos%y, 0._wp)
        res = this%prim%evaluate(q)
    
    end function eval_revolution

    type(onion) function onion_init(prim, thickness) result(out)

        class(sdf_base), target :: prim
        real(kind=wp), intent(IN) :: thickness

        out%thickness = thickness
        out%prim => prim

        out%optProps = prim%optProps

        out%layer = prim%layer
        out%transform = identity()

    end function onion_init

    pure elemental function eval_onion(this, pos) result(res)

        class(onion), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        res = abs(this%prim%evaluate(pos)) - this%thickness

    end function eval_onion

    pure elemental function eval_elongate(this, pos) result(res)

        class(elongate), intent(in) :: this
        type(vector),    intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: w
        type(vector) :: q

        q = abs(pos) - this%size
        w = min(max(q%x, max(q%y, q%z)), 0._wp)

        res = this%prim%evaluate(max(q, 0._wp)) + w

    end function eval_elongate

    pure elemental function eval_twist(this, pos) result(res)

        class(twist), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: c, s, x2, y2, z2

        c = cos(this%k * pos%z)
        s = sin(this%k * pos%z)
        x2 = c*pos%x - s*pos%y
        y2 = s*pos%x + c*pos%y
        z2 = pos%z

        res = this%prim%evaluate(vector(x2, y2, z2))

    end function eval_twist

    pure elemental function eval_bend(this, pos) result(res)

        class(bend),  intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: c, s, x2, y2, z2

        c = cos(this%k * pos%x)
        s = sin(this%k * pos%x)
        x2 = c * pos%x - s * pos%y
        y2 = s * pos%x + c * pos%y
        z2 = pos%z

        res = this%prim%evaluate(vector(x2, y2, z2))

    end function eval_bend

    pure elemental function eval_disp(this, pos) result(res)

        class(displacement), intent(in) :: this
        type(vector),        intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: d1, d2

        d1 = this%prim%evaluate(pos)
        d2 = this%func(pos)

        res = d1 + d2

    end function eval_disp

    pure elemental function eval_repeat(this, pos) result(res)
        
        ! use utils, only : clamp

        class(repeat), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: q

        error stop "Not implmented as no vector dependacny in utils yet!"
        ! q = pos - this%c*clamp(nint(pos/this%c), this%la, this%lb)
        res = this%prim%evaluate(q)

    end function eval_repeat

    pure function union(d1, d2, k) result(res)

        real(kind=wp), intent(IN) :: d1, d2, k
        real(kind=wp) :: res

        res = min(d1, d2)
    end function union

    pure function SmoothUnion(d1, d2, k) result(res)

        real(kind=wp), intent(IN) :: d1, d2, k
        real(kind=wp) :: res

        real(kind=wp) :: h

        h = max(k - abs(d1 - d2), 0._wp) / k
        res = min(d1, d2) - h*h*h*k*(1._wp/6._wp)

    end function SmoothUnion

    pure function subtraction(d1, d2, k) result(res)

        real(kind=wp), intent(IN) :: d1, d2, k
        real(kind=wp) :: res

        res = max(-d1, d2)

    end function subtraction

    pure function intersection(d1, d2, k) result(res)

        real(kind=wp), intent(IN) :: d1, d2, k
        real(kind=wp) :: res

        res = max(d1, d2)

    end function intersection
end module sdfModifiers