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

    contains


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