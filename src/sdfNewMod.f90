!https://fortran-lang.discourse.group/t/attempting-type-erasure-in-fortran/4402/2
module sdfNew

    use vector_class
    use constants, only : wp
    use opticalProperties, only : opticalProp_t
    use sdfHelpers, only : identity

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

    type, extends(sdf) :: triprism
        real(kind=wp) :: h1, h2
        contains
        procedure :: evaluate => evaluate_triprism
    end type triprism

    type, extends(sdf) :: cone
        type(vector)  :: a, b
        real(kind=wp) :: ra, rb
        contains
        procedure :: evaluate => evaluate_cone
    end type cone

    type, extends(sdf) :: capsule
        type(vector)  :: a, b
        real(kind=wp) :: r
        contains
        procedure :: evaluate => evaluate_capsule
    end type capsule

    type, extends(sdf) :: plane
        type(vector) :: a
        contains
        procedure :: evaluate => evaluate_plane
    end type plane

    type, extends(sdf) :: segment
        type(vector) :: a, b
        contains
        procedure :: evaluate => evaluate_segment
    end type segment

    type, extends(sdf) :: egg
        real(kind=wp) :: r1, r2, h
        contains
        procedure :: evaluate => evaluate_egg
    end type egg

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
 
    interface sdf
       module procedure sdf_new
    end interface

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

    contains

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

        res = this%value%optProps%p%kappa

    end function getKappa

    function getMua(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%p%mua

    end function getMua

    function gethgg(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%p%hgg

    end function gethgg


    function getg2(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%p%g2

    end function getg2
    function getN(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%p%n

    end function getN

    function getAlbedo(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%p%albedo

    end function getAlbedo

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

            class(egg), intent(in) :: this
            type(vector),  intent(IN) :: pos
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

    pure elemental function sdf_evaluate(this, pos) result(res)

       class(sdf), intent(in) :: this
       type(vector), intent(in) :: pos
       real(kind=wp) :: res

       res = this%value%evaluate(pos)

    end function
 
    ! sdf initializer
    subroutine sdf_assign(lhs, rhs)

        class(sdf),      intent(inout) :: lhs
        class(sdf_base), intent(in)    :: rhs

        if (allocated(lhs%value))deallocate(lhs%value)
        ! Prevent nested derived type
        select type (rhsT=>rhs)
            class is (sdf)
                if (allocated(rhsT%value)) allocate(lhs%value,source=rhsT%value)
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