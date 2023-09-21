module sdfs
    
    !! This module defines the signed distance function (SDF) abstract type and all types that inherit from it.
    !! The SDF abstract type defines the optical properties of an SDF (mus, mua, kappa, albedo, hgg, g2,and n), as well as a transform (4x4 matrix), and the layer ID code of the SDF.
    !! The SDF abstract type also provides an abstract interface (evaluate) which each inheriting function must implement. This evaluate function is the heart of the SDF implementation.
    !! Each individual evaluate is the direct implementation of that SDF, e.g. that function defines the mathematical SDF.
    !! For more information on SDFs, check out Inigo Quilez's [website](https://iquilezles.org/articles/) from which most of the below SDFs and transforms have been taken.
    
    !! - cylinder
    !! - sphere
    !! - box
    !! - torus
    !! - cone
    !! - triprism (triangular prism)
    !! - capsule
    !! - plane
    !! - segment
    !! - egg
    
    !! **This is the module the user should import to other module not sdf_base!**

    use constants,         only : wp
    use opticalProperties, only : opticalProp_t
    use sdf_baseMod,       only : sdf, sdf_base, model, calcNormal, render
    use sdfHelpers,        only : identity
    use vector_class

    implicit none
    
    !> Box SDF
    type, extends(sdf_base) :: box
        !> Length of each dimension of the box
        type(vector) :: lengths
        contains
        procedure :: evaluate => evaluate_box
    end type box

    !> Sphere SDF
    type, extends(sdf_base) :: sphere
        real(kind=wp) :: radius
        contains
        procedure :: evaluate => evaluate_sphere
    end type sphere

    !> Cylinder SDF
    type, extends(sdf_base) :: cylinder
        real(kind=wp) :: radius
        type(vector) :: a, b
        contains
        procedure :: evaluate => evaluate_cylinder
    end type cylinder

    !> Torus SDF
    type, extends(sdf_base) :: torus
        real(kind=wp) :: oradius, iradius
        contains
        procedure :: evaluate => evaluate_torus
    end type torus

    !> Triprisim SDF
    type, extends(sdf_base) :: triprism
        real(kind=wp) :: h1, h2
        contains
        procedure :: evaluate => evaluate_triprism
    end type triprism

    !> Cone SDF
    type, extends(sdf_base) :: cone
        type(vector)  :: a, b
        real(kind=wp) :: ra, rb
        contains
        procedure :: evaluate => evaluate_cone
    end type cone

    !> Capsule SDF
    type, extends(sdf_base) :: capsule
        type(vector)  :: a, b
        real(kind=wp) :: r
        contains
        procedure :: evaluate => evaluate_capsule
    end type capsule

    !> Plane SDF
    type, extends(sdf_base) :: plane
        type(vector) :: a
        contains
        procedure :: evaluate => evaluate_plane
    end type plane

    !> Segment SDF (2D)
    type, extends(sdf_base) :: segment
        type(vector) :: a, b
        contains
        procedure :: evaluate => evaluate_segment
    end type segment

    !> Egg SDF
    type, extends(sdf_base) :: egg
        real(kind=wp) :: r1, r2, h
        contains
        procedure :: evaluate => evaluate_egg
    end type egg

    interface sphere
        module procedure sphere_init
    end interface sphere

    interface box
        !! Interface to box SDF initialising function
        module procedure box_init
    end interface box

    interface torus
        !! Interface to torus SDF initialising function
        module procedure torus_init
    end interface torus

    interface cylinder
        !! Interface to cylinder SDF initialising function
        module procedure cylinder_init
    end interface cylinder

    interface triprism
        !! Interface to triprisim SDF initialising function
        module procedure triprism_init
    end interface triprism

    interface egg
        !! Interface to egg SDF initialising function
        module procedure egg_init
    end interface egg

    interface segment
        !! Interface to segment SDF initialising function
        module procedure segment_init
    end interface segment

    interface cone
        !! Interface to cone SDF initialising function
        module procedure cone_init
    end interface cone

    interface capsule
        !! Interface to capsule SDF initialising function
        module procedure capsule_init
    end interface capsule

    interface plane
        !! Interface to plane SDF initialising function
        module procedure plane_init
    end interface plane

    private
    public :: plane, capsule, cone, segment, egg, triprism, cylinder, torus, box, sphere, sdf, model, calcNormal, render

    contains
    
    function segment_init(a, b, optProp, layer, transform) result(out)
        !! Initalising function for segment SDF.
        !! Note this is a 2D function

        type(segment) :: out

        !> Optical properties of the SDF
        type(opticalProp_t),     intent(in) :: optProp
        !> segment start point
        type(vector),            intent(IN) :: a
        !> segment end point
        type(vector),            intent(IN) :: b
        !> ID number of sdf
        integer,                 intent(IN) :: layer
        !> Optional transform to apply to SDF
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
        !! Initalising function for egg SDF.
        !! makes a Moss egg. [ref](https://www.shadertoy.com/view/WsjfRt).

        type(egg) :: out
        
        !> R1 controls "fatness" of the egg. Actually controls the base circle radius.
        real(kind=wp),            intent(IN) :: r1
        !> R2 contorls the pointiness of the egg. Actually controls radius of top circle.
        real(kind=wp),            intent(in) :: r2
        !> h controls the height of the egg. Actually controls y position of top circle.
        real(kind=wp),            intent(in) :: h
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for plane SDF.

        type(plane) :: out
        
        !> Plane normal. must be normalised
        type(vector),             intent(IN) :: a
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for capsule SDF.

        type(capsule) :: out
        
        !> Capsule startpoint
        type(vector),            intent(IN) :: a
        !> Capsule endpoint
        type(vector),            intent(IN) :: b
        !> Capsule radius
        real(kind=wp),           intent(IN) :: r
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for triprisim SDF.
        
        type(triprism) :: out
        
        !> Height of triprisim
        real(kind=wp),            intent(IN) :: h1
        !> length of triprisim
        real(kind=wp),            intent(IN) :: h2
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for Capped Cone SDF.

        type(cone) :: out
        
        !> Centre of base of Cone
        type(vector),             intent(IN) :: a
        !> Tip of cone
        type(vector),             intent(IN) :: b
        !> Radius of Cones base
        real(kind=wp),            intent(IN) :: ra
        !> Radius of Cones tip. For rb = 0.0 get normal uncapped cone.
        real(kind=wp),            intent(in) :: rb
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for Cylinder SDF.

        type(cylinder) :: out
        
        !> Radius of cylinder
        real(kind=wp),           intent(in) :: radius
        !> Vector position at centre of the bottom circle
        type(vector),            intent(IN) :: a
        !> Vector position at centre of the top circle
        type(vector),            intent(IN) :: b
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for Torus SDF.

        type(torus) :: out
        !> Outer radius of Torus
        real(kind=wp),           intent(IN) :: oradius
        !> Inner radius of Torus
        real(kind=wp),           intent(IN) :: iradius
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for Box SDF.

        type(box) :: out
        
        !> Lengths of each dimension of the box
        type(vector),             intent(IN) :: lengths
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Initalising function for Sphere SDF.

        type(sphere) :: out
        
        !> radius of the Sphere
        real(kind=wp),            intent(IN) :: radius
        !> ID number of sdf
        integer,                  intent(IN) :: layer
        !> Optional transform to apply to SDF
        real(kind=wp),  optional, intent(IN) :: transform(4, 4)
        !> Optical properties of the SDF
        type(opticalProp_t),      intent(in) :: optProp

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
        !! Evaluation function for Sphere SDF.

        class(sphere), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector),  intent(in) :: pos

        real(kind=wp) :: res

        type(vector) :: p

        p = pos .dot. this%transform
        res = sqrt(p%x**2+p%y**2+p%z**2) - this%radius

    end function evaluate_sphere
 
    pure elemental function evaluate_box(this, pos) result(res)
        !! Evaluation function for Box SDF.

        class(box),   intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector), intent(in) :: pos

        real(kind=wp) :: res

        type(vector) :: p, q

        p = pos .dot. this%transform
        q = abs(p) - this%lengths
        res = length(max(q, 0._wp)) + min(max(q%x, max(q%y, q%z)), 0._wp)

    end function evaluate_box

    pure elemental function evaluate_torus(this, pos) result(res)
        !! Evaluation function for Torus SDF.

        class(torus), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector), intent(in) :: pos

        real(kind=wp) :: res

        type(vector) :: p, q

        p = pos .dot. this%transform
        q = vector(length(vector(p%x, 0._wp, p%z)) - this%oradius, p%y, 0._wp)
        res = length(q) - this%iradius

    end function evaluate_torus

    pure elemental function evaluate_cylinder(this, pos) result(res)
        !! Evaluation function for Cylinder SDF.

        class(cylinder), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector),    intent(in) :: pos
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
        !! Evaluation function for Triprisim SDF.

        class(triprism), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector),    intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: q, p

        p = pos .dot. this%transform
        q = abs(p)
        res = max(q%z - this%h2, max(q%x*.866025_wp + p%y*.5_wp, -p%y) - this%h1*.5_wp) 

    end function evaluate_triprism

    pure elemental function evaluate_segment(this, pos) result(res)
        !! Evaluation function for Segment SDF.

        !p = pos
        !a = pt1
        !b = pt2
        !draws segment along the axis between 2 points a and b

        use utils, only : clamp

        class(segment), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector),   intent(IN) :: pos

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
        !! Evaluation function for Capsule SDF.

        use utils, only : clamp

        class(capsule), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector),   intent(in) :: pos
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
        !! Evaluation function for Cone SDF.

        use utils, only : clamp

        class(cone),  intent(in) :: this
        type(vector), intent(IN) :: pos
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
        !! Evaluation function for Egg SDF.
        !! [ref](https://www.shadertoy.com/view/WsjfRt)

        class(egg),   intent(in) :: this
        !> vector position to evaluate SDF at
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
        !! Evaluation function for Plane SDF.

        class(plane), intent(in) :: this
        !> vector position to evaluate SDF at
        type(vector), intent(IN) :: pos
        real(kind=wp) :: res

        type(vector) :: p

        p = pos .dot. this%transform

        !a must be normalised
        res = (p .dot. this%a)

    end function evaluate_plane
end module sdfs