!   Module provides signed distance functions (SDFs) for various shapes 
!   and some operations to modify them
!   All SDF functions are adapted from Inigo Quilex exhaustive list at:
!   https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm!
!   API based upon https://fortran-lang.discourse.group/t/attempting-type-erasure-in-fortran/4402/2
module sdfNew

    use constants,         only : wp
    use opticalProperties, only : opticalProp_t
    use sdfHelpers,        only : identity
    use vector_class

    ! displacement, bend, twist, elongate, repeat

    implicit none

    type, abstract :: sdf_base
        type(opticalProp_t) :: optProps
        real(kind=wp) :: transform(4, 4)
        integer :: layer
        contains
            procedure(evalInterface), deferred :: evaluate
    end type sdf_base
        
    type, public, extends(sdf_base) :: sdf
        class(sdf_base), allocatable :: value
        contains
            procedure :: getKappa
            procedure :: getAlbedo
            procedure :: getMua, gethgg, getG2, getN
            procedure :: evaluate => sdf_evaluate
            procedure, private :: sdf_assign
            generic :: assignment(=) => sdf_assign
    end type sdf

!    shapes
!################################################################
    type, public, extends(sdf_base) :: box
        type(vector) :: lengths
        contains
        procedure :: evaluate => evaluate_box
    end type box

    type, public, extends(sdf_base) :: sphere
        real(kind=wp) :: radius
        contains
        procedure :: evaluate => evaluate_sphere
    end type sphere

    type, public, extends(sdf_base) :: cylinder
        real(kind=wp) :: radius
        type(vector) :: a, b
        contains
        procedure :: evaluate => evaluate_cylinder
    end type cylinder

    type, public, extends(sdf_base) :: torus
        real(kind=wp) :: oradius, iradius
        contains
        procedure :: evaluate => evaluate_torus
    end type torus

    type, extends(sdf_base) :: triprism
        real(kind=wp) :: h1, h2
        contains
        procedure :: evaluate => evaluate_triprism
    end type triprism

    type, extends(sdf_base) :: cone
        type(vector)  :: a, b
        real(kind=wp) :: ra, rb
        contains
        procedure :: evaluate => evaluate_cone
    end type cone

    type, extends(sdf_base) :: capsule
        type(vector)  :: a, b
        real(kind=wp) :: r
        contains
        procedure :: evaluate => evaluate_capsule
    end type capsule

    type, extends(sdf_base) :: plane
        type(vector) :: a
        contains
        procedure :: evaluate => evaluate_plane
    end type plane

    type, extends(sdf_base) :: segment
        type(vector) :: a, b
        contains
        procedure :: evaluate => evaluate_segment
    end type segment

    type, extends(sdf_base) :: egg
        real(kind=wp) :: r1, r2, h
        contains
        procedure :: evaluate => evaluate_egg
    end type egg
!####################################################################
!   modifiers
!####################################################################
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
!####################################################################
    
    abstract interface
        pure elemental function evalInterface(this, pos) result(res)
            use vector_class
            use constants, only : wp
            import sdf_base
            class(sdf_base), intent(in) :: this
            type(vector),    intent(in) :: pos
            real(kind=wp) :: res
        end function
    end interface

    abstract interface
        pure function primitive(pos) result(res)
            use vector_class, only : vector
            use constants,    only : wp
            implicit none
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res
        end function primitive
    end interface

    interface sdf
        module procedure sdf_new
    end interface
    
    interface revolution
        module procedure revolution_init
    end interface revolution

    interface extrude
        module procedure extrude_init
    end interface extrude

    interface onion
        module procedure onion_init
    end interface onion
    
    interface sphere
        module procedure sphere_init
    end interface sphere

    interface box
        module procedure box_init
    end interface box

    interface torus
        module procedure torus_init
    end interface torus

    interface cylinder
        module procedure cylinder_init
    end interface cylinder

    interface triprism
        module procedure triprism_init
    end interface triprism

    interface egg
        module procedure egg_init
    end interface egg

    interface segment
        module procedure segment_init
    end interface segment

    interface cone
        module procedure cone_init
    end interface cone

    interface capsule
        module procedure capsule_init
    end interface capsule

    interface plane
        module procedure plane_init
    end interface plane

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

    function segment_init(a, b, optProp, layer, transform) result(out)

        type(segment) :: out
        type(opticalProp_t), intent(in) :: optProp
        type(vector),            intent(IN) :: a, b
        ! real(kind=wp),           intent(IN) :: mus, mua, hgg, n
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%a = a
        out%b = b

        out%layer = layer
        out%transform = t
        
        out%optProps = optProp

    end function segment_init

    function egg_init(r1, r2, h, optProp, layer, transform) result(out)
        ! makes a Moss egg. https://www.shadertoy.com/view/WsjfRt
        ! R1 controls "fatness" of the egg. Actually controls the base circle radius.
        ! R2 contorls the pointiness of the egg. Actually controls radius of top circle.
        ! h controls the height of the egg. Actually controls y position of top circle.
        type(egg) :: out
        
        type(opticalProp_t),      intent(in) :: optProp
        real(kind=wp),            intent(IN) :: r1, r2, h
        integer,                  intent(IN) :: layer
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%h = h
        out%r1 = r1
        out%r2 = r2
        out%layer = layer
        out%transform = t
        out%optProps = optProp

    end function egg_init

    function plane_init(a, optProp, layer, transform) result(out)
        
        type(plane) :: out
        
        type(opticalProp_t),      intent(in) :: optProp
        type(vector),             intent(IN) :: a
        integer,                  intent(IN) :: layer
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%a = a
        out%layer = layer
        out%transform = t

        out%optProps = optProp

    end function plane_init

    function capsule_init(a, b, r, optProp, layer, transform) result(out)
        
        type(capsule) :: out
        
        type(vector),            intent(IN) :: a, b
        type(opticalProp_t),     intent(in) :: optProp
        real(kind=wp),           intent(IN) :: r
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%a = a
        out%b = b
        out%r = r
        out%layer = layer
        out%transform = t

        out%optProps = optProp

    end function capsule_init

    function triprism_init(h1, h2, optProp, layer, transform) result(out)
        !h1 is height
        !h2 is length
        !        
        type(triprism) :: out
        
        type(opticalProp_t),     intent(in) :: optProp
        real(kind=wp),           intent(IN) :: h1, h2
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%h1 = h1
        out%h2 = h2
        out%layer = layer
        out%transform = t
        out%optProps = optProp

    end function triprism_init

    function cone_init(a, b, ra, rb, optProp, layer, transform) result(out)
        
        type(cone) :: out
        
        type(opticalProp_t),     intent(in) :: optProp
        type(vector),            intent(IN) :: a, b
        real(kind=wp),           intent(IN) :: ra, rb
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%a = a
        out%b = b
        out%ra = ra
        out%rb = rb
        out%layer = layer
        out%transform = t
        out%optProps = optProp

    end function cone_init

    function cylinder_init(a, b, radius, optProp, layer, transform) result(out)
                
        type(cylinder) :: out
        
        type(opticalProp_t),     intent(in) :: optProp
        real(kind=wp),           intent(in) :: radius
        integer,                 intent(IN) :: layer
        type(vector),            intent(IN) :: a, b
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%a = a
        out%b = b
        out%radius = radius
        out%layer = layer
        out%transform = t

        out%optProps = optProp

    end function cylinder_init

    function torus_init(oradius, iradius, optProp, layer, transform) result(out)
        
        type(torus) :: out
        
        type(opticalProp_t),     intent(in) :: optProp
        real(kind=wp),           intent(IN) :: oradius, iradius
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%oradius = oradius
        out%iradius = iradius
        out%layer = layer
        out%transform = t

        out%optProps = optProp

    end function torus_init

    function box_init(lengths, optProp, layer, transform) result(out)
                
        type(box) :: out
        
        type(vector),            intent(IN) :: lengths
        type(opticalProp_t),     intent(in) :: optProp
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(in) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%lengths = .5_wp*lengths! as only half lengths
        out%layer = layer
        out%transform = t

        out%optProps = optProp

    end function box_init
    
    function sphere_init(radius, optProp, layer, transform) result(out)
                
        type(sphere) :: out
        
        type(opticalProp_t),     intent(in) :: optProp
        real(kind=wp),           intent(IN) :: radius
        integer,                 intent(IN) :: layer
        real(kind=wp), optional, intent(IN) :: transform(4, 4)

        real(kind=wp) :: t(4, 4)

        if(present(transform))then
            t = transform
        else
            t = identity()
        end if

        out%radius = radius
        out%layer = layer

        out%transform = t

        out%optProps = optProp

    end function sphere_init

    pure elemental function evaluate_sphere(this, pos) result(res)
       class(sphere), intent(in) :: this
       type(vector), intent(in) :: pos
       real(kind=wp) :: res

        type(vector) :: p

        p = pos .dot. this%transform
        res = sqrt(p%x**2+p%y**2+p%z**2) - this%radius

    end function evaluate_sphere
 
    pure elemental function evaluate_box(this, pos) result(res)
       class(box), intent(in) :: this
       type(vector), intent(in) :: pos
       real(kind=wp) :: res

        type(vector) :: p, q

        p = pos .dot. this%transform
        q = abs(p) - this%lengths
        res = length(max(q, 0._wp)) + min(max(q%x, max(q%y, q%z)), 0._wp)

    end function evaluate_box

    pure elemental function evaluate_torus(this, pos) result(res)
       class(torus), intent(in) :: this
       type(vector), intent(in) :: pos
       real(kind=wp) :: res

        type(vector) :: p, q

        p = pos .dot. this%transform
        q = vector(length(vector(p%x, 0._wp, p%z)) - this%oradius, p%y, 0._wp)
        res = length(q) - this%iradius

    end function evaluate_torus

    pure elemental function evaluate_cylinder(this, pos) result(res)
       class(cylinder), intent(in) :: this
       type(vector), intent(in) :: pos
       real(kind=wp) :: res

        type(vector)  :: p, ba, pa
        real(kind=wp) :: x, y, x2, y2, d, baba, paba

        p = pos .dot. this%transform

        ba = this%b - this%a
        pa = p - this%a
        baba = ba .dot. ba
        paba = pa .dot. ba
        x = length(pa * baba - ba*paba) - this%radius*baba
        y = abs(paba - baba*.5_wp) - baba*.5_wp
        x2 = x**2
        y2 = (y**2)*baba
        if(max(x, y) < 0._wp)then
            d = -min(x2, y2)
        else
            if(x > 0._wp .and. y > 0._wp)then
                d = x2 + y2
            elseif(x > 0._wp)then
                d = x2
            elseif(y > 0._wp)then
                d = y2
            else
                d = 0._wp
            end if
        end if

        res = sign(sqrt(abs(d))/baba, d)

    end function evaluate_cylinder

    pure elemental function evaluate_triprism(this, pos) result(res)

        class(triprism), intent(in) :: this
        type(vector),  intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: q, p

        p = pos .dot. this%transform
        q = abs(p)
        res = max(q%z - this%h2, max(q%x*.866025_wp + p%y*.5_wp, -p%y) - this%h1*.5_wp) 

    end function evaluate_triprism

    pure elemental function evaluate_segment(this, pos) result(res)
        !p = pos
        !a = pt1
        !b = pt2
        !draws segment along the axis between 2 points a and b

        use utils, only : clamp

        class(segment), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        type(vector)  :: pa, ba, p
        real(kind=wp) :: h
       
        p = pos .dot. this%transform

        pa = p - this%a
        ba = this%b - this%a
        h = clamp((pa .dot. ba) / (ba .dot. ba), 0.0_wp, 1.0_wp)

        res = length(pa - ba*h) - 0.1_wp

    end function evaluate_segment

    pure elemental function evaluate_capsule(this, pos) result(res)

        use utils, only : clamp

        class(capsule), intent(in) :: this
        type(vector), intent(in) :: pos
        real(kind=wp) :: res

        type(vector) :: pa, ba, p
        real(kind=wp) :: h

        p = pos .dot. this%transform

        pa = p - this%a
        ba = this%b - this%a
        h = clamp((pa .dot. ba) / (ba .dot. ba), 0._wp, 1._wp)
        res = length(pa - ba*h) - this%r

    end function evaluate_capsule

    pure elemental function evaluate_cone(this, pos) result(res)

        use utils, only : clamp

        class(cone), intent(in) :: this
        type(vector),  intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: rba, baba, papa, paba, x, cax, cay, k, f, cbx, cby, s
        type(vector) :: p

        p = pos .dot. this%transform

        rba = this%rb - this%ra
        baba = (this%b-this%a) .dot. (this%b-this%a)
        papa = (p-this%a) .dot. (p-this%a)
        paba =  ((p-this%a) .dot. (this%b-this%a))/ baba
        x = sqrt(papa - baba*paba**2)
        if(paba < 0.5_wp)then
            cax = max(0._wp, x - this%ra)
        else
            cax = max(0._wp, x - this%rb)
        end if
        cay = abs(paba - 0.5_wp) - .5_wp
        k = rba**2 + baba
        f = clamp((rba * (x - this%ra) + paba*baba) / k, 0._wp, 1._wp)
        cbx = x - this%ra - f*rba
        cby = paba - f
        if(cbx < 0._wp .and. cay < 0._wp)then
            s = -1._wp
        else
            s = 1._wp
        end if 
        res = s * sqrt(min(cax**2 + baba*cay**2, cbx**2 + baba*cby**2)) 

    end function evaluate_cone

    pure elemental function evaluate_egg(this, pos) result(res)
        !https://www.shadertoy.com/view/WsjfRt

        class(egg),   intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        real(kind=wp) :: r, l, h_in
        type(vector) :: p_in, p

        p = pos .dot. this%transform

        p_in = p

        p_in%x = abs(p%x)
        r = this%r1 - this%r2
        h_in = this%h + r
        l = (h_in**2 - r**2) / (2._wp * r)

        if(p_in%y <= 0._wp)then
            res = length(p_in) - this%r1
        else
            if((p_in%y - h_in) * l > p_in%x*h_in)then
                res = length(p_in - vector(0._wp, h_in, 0._wp)) - ((this%r1+l) - length(vector(h_in,l, 0._wp)))
            else
                res = length(p_in + vector(l, 0._wp, 0._wp)) - (this%r1+l)
            end if
        end if
    end function evaluate_egg

    pure elemental function evaluate_plane(this, pos) result(res)

        class(plane), intent(in) :: this
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: p

        p = pos .dot. this%transform

        !a must be normalised
        res = (p .dot. this%a)

    end function evaluate_plane

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

!#############################################################
!           Helpers
!#############################################################
    type(vector) function calcNormal(p, obj)

        type(vector), intent(IN) :: p
        class(sdf_base) :: obj

        real(kind=wp) :: h
        type(vector) :: xyy, yyx, yxy, xxx

        h = 1e-6_wp
        xyy = vector(1._wp, -1._wp, -1._wp)
        yyx = vector(-1._wp, -1._wp, 1._wp)
        yxy = vector(-1._wp, 1._wp, -1._wp)
        xxx = vector(1._wp, 1._wp, 1._wp)

        calcNormal = xyy*obj%evaluate(p + xyy*h) +  &
                    yyx*obj%evaluate(p + yyx*h) +  &
                    yxy*obj%evaluate(p + yxy*h) +  &
                    xxx*obj%evaluate(p + xxx*h)

        calcNormal = calcNormal%magnitude()

    end function calcNormal

    function getKappa(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%kappa

    end function getKappa

    function getMua(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%mua

    end function getMua

    function gethgg(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%hgg

    end function gethgg


    function getg2(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%g2

    end function getg2
    function getN(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%n

    end function getN

    function getAlbedo(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%albedo

    end function getAlbedo
!#########################################################################
!       SDF bound procedures
!#########################################################################
    pure elemental function sdf_evaluate(this, pos) result(res)

        class(sdf), intent(in) :: this
        type(vector), intent(in) :: pos
        real(kind=wp) :: res

        res = this%value%evaluate(pos)

    end function sdf_evaluate
 
    ! sdf initializer
    subroutine sdf_assign(lhs, rhs)

        class(sdf),      intent(inout) :: lhs
        class(sdf_base), intent(in)    :: rhs

        if (allocated(lhs%value))deallocate(lhs%value)
        ! Prevent nested derived type
        select type (rhsT=>rhs)
            class is (sdf)
                if(allocated(rhsT%value))allocate(lhs%value,source=rhsT%value)
            class default
                allocate(lhs%value,source=rhsT)
        end select
    end subroutine sdf_assign
 
    ! sdf initializer
    type(sdf) function sdf_new(rhs) result(lhs)

        class(sdf_base), intent(in) :: rhs
        allocate(lhs%value,source=rhs)

    end function sdf_new
end module sdfNew