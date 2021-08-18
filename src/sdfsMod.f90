module sdfs
!   Module provides signed distance functions (SDFs) for various shapes 
!   and some operations to adjust, move, rotate etc them
!   All SDF functions are adapted from Inigo Quilex exhaustive list at:
!   https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
!
    use vector_class

    implicit none

    type, abstract :: sdf
    ! base class.
    ! provides deferred evaluation func to get distance to shape
    ! and gives the transform and optical prperties.
    ! Layer is an important bookeeping integer.
        real :: mus, mua, kappa, albedo, hgg, g2, n
        real :: transform(4, 4)
        integer :: layer
        contains
        procedure(evalInterface), deferred :: evaluate
    end type sdf

    abstract interface
        real function evalInterface(this, pos)
            use vector_class
            import sdf
            implicit none
            class(sdf) :: this
            type(vector), intent(IN) :: pos
        end function evalInterface
    end interface



    type, extends(sdf) :: sphere
        real :: radius
        contains
        procedure :: evaluate => eval_sphere
    end type sphere

    interface sphere
        module procedure sphere_init
    end interface sphere


    type, extends(sdf) :: cylinder
        real         :: radius
        type(vector) :: a, b
        contains
        procedure :: evaluate => eval_cylinder
    end type cylinder

    interface cylinder
        module procedure cylinder_init
    end interface cylinder


    type, extends(sdf) :: box
        type(vector) :: lengths
        contains
        procedure :: evaluate => eval_box
    end type box

    interface box
        module procedure box_init_vec
        module procedure box_init_scal
    end interface box


    type, extends(sdf) :: torus
        real :: oradius, iradius
        contains
        procedure :: evaluate => eval_torus
    end type torus

    interface torus
        module procedure torus_init
    end interface torus


    type, extends(sdf) :: triprisim
        real :: h1, h2
        contains
        procedure :: evaluate => eval_triprisim
    end type triprisim

    interface triprisim
        module procedure triprisim_init
    end interface triprisim


    type, extends(sdf) :: cone
        type(vector) :: a, b
        real         :: ra, rb
        contains
        procedure :: evaluate => eval_cone
    end type cone

    interface cone
        module procedure cone_init
    end interface cone


    type, extends(sdf) :: model
        type(container), allocatable  :: array(:)
        procedure(op),nopass, pointer :: func
        contains
        procedure :: evaluate => eval_model
    end type model

    interface model
        module procedure model_init
    end interface model


    type :: container
        class(sdf), pointer :: p => null()
    end type container




    type, extends(sdf) :: twist
        real :: k
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_twist
    end type twist

    interface twist
        module procedure twist_init
    end interface twist

    type, extends(sdf) :: displacement
        procedure(primitive), nopass, pointer :: func
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_disp
    end type displacement

    interface displacement
        module procedure displacement_init
    end interface displacement

    type, extends(sdf) :: bend
        real :: k
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_bend
    end type bend

    interface bend
        module procedure bend_init
    end interface bend

    type, extends(sdf) :: elongate
        type(vector) :: size
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_elongate
    end type elongate

    interface elongate
        module procedure elongate_init
    end interface elongate



    abstract interface
        real function op(d1, d2)
            implicit none
            real, intent(IN) :: d1, d2
        end function op
    end interface

    abstract interface
        real function primitive(pos)
            use vector_class, only : vector
            implicit none
            type(vector), intent(IN) :: pos
        end function primitive
    end interface

    private
    ! shapes
    public :: sdf, cylinder, sphere, box, torus, cone, triprisim
    ! meta
    public :: model, container
    ! boolean ops
    public :: union, intersection, subtraction, SmoothUnion, op
    ! move ops
    public :: rotate_x, rotate_y, rotate_z, identity, translate
    ! deform ops
    public :: displacement, bend, twist, elongate
    ! utility funcs
    public :: calcNormal, model_init, render

    contains

        type(vector) function calcNormal(p, obj)

            implicit none

            type(vector), intent(IN) :: p
            class(sdf) :: obj

            real :: h
            type(vector) :: xyy, yyx, yxy, xxx

            h = 1d-10
            xyy = vector(1., -1., -1.)
            yyx = vector(-1., -1., 1.)
            yxy = vector(-1., 1., -1.)
            xxx = vector(1., 1., 1.)

            calcNormal = xyy*obj%evaluate(p + xyy*h) +  &
                         yyx*obj%evaluate(p + yyx*h) +  &
                         yxy*obj%evaluate(p + yxy*h) +  &
                         xxx*obj%evaluate(p + xxx*h)

            calcNormal = calcNormal%magnitude()

        end function calcNormal

        function model_init(array, func) result(out)
        !TODO make sure optical properties are same in all inputs            
            implicit none

            type(model) :: out

            procedure(op) :: func
            type(container), intent(IN) :: array(:)

            out%array = array
            out%func => func

            associate(member => array(1)%p)
                out%mus = member%mus
                out%mua = member%mua
                out%kappa = member%kappa
                out%albedo = member%albedo
                out%g2 = member%g2
                out%hgg = member%hgg
                out%n = member%n
                out%layer = member%layer
            end associate

        end function model_init

        real function eval_model(this, pos)

            implicit none

            class(model) :: this
            type(vector), intent(IN) :: pos

            integer :: i

            eval_model = this%array(1)%p%evaluate(pos)
            do i = 2, size(this%array)
                eval_model = this%func(eval_model, this%array(i)%p%evaluate(pos))
            end do

        end function eval_model


        function cylinder_init(a, b, radius, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(cylinder) :: out
            
            real,           intent(IN) :: radius, mus, mua, hgg, n
            integer,        intent(IN) :: layer
            type(vector),   intent(IN) :: a, b
            real, optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function cylinder_init

        real function eval_cylinder(this, pos)

            implicit none

            class(cylinder) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_cylinder = cylinder_fn(p, this%a, this%b, this%radius)

        end function eval_cylinder

        real function cylinder_fn(p, a, b, r)
            !p = pos
            !a = pt1
            !b = pt2
            !r = radius
            !draws cylinder along the axis between 2 points a and b
            implicit none

            type(vector), intent(IN) :: p, a, b
            real,         intent(IN) :: r

            type(vector) :: ba, pa
            real :: x, y, x2, y2, d, baba, paba

            ba = b - a
            pa = p - a
            baba = ba .dot. ba
            paba = pa .dot. ba
            x = length(pa * baba - ba*paba) - r*baba
            y = abs(paba - baba*.5) - baba*.5
            x2 = x**2
            y2 = (y**2)*baba
            if(max(x, y) < 0.)then
                d = -min(x2, y2)
            else
                if(x > 0. .and. y > 0.)then
                    d = x2 + y2
                elseif(x > 0.)then
                    d = x2
                elseif(y > 0)then
                    d = y2
                else
                    d = 0.
                end if
            end if

            cylinder_fn = sign(sqrt(abs(d))/baba, d)

        end function cylinder_fn


        function box_init(lengths, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            type(vector),   intent(IN) :: lengths
            real,           intent(IN) :: mus, mua, hgg, n
            integer,        intent(IN) :: layer
            real, optional, intent(in) :: transform(4, 4)

            real :: t(4, 4)


            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%lengths = .5*lengths! as only half lengths
            out%layer = layer
            out%transform = t

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function box_init

        function box_init_vec(lengths, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            type(vector),   intent(IN) :: lengths
            real,           intent(IN) :: mus, mua, hgg, n
            integer,        intent(IN) :: layer
            real, optional, intent(in) :: transform(4, 4)

            out = box_init(lengths, mus, mua, hgg, n, layer, transform)

        end function box_init_vec

        function box_init_scal(length, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            real,           intent(IN) :: length, mus, mua, hgg, n
            integer,        intent(IN) :: layer
            real, optional, intent(in) :: transform(4, 4)

            type(vector) :: lengths

            lengths = vector(length, length, length)

            out = box_init(lengths, mus, mua, hgg, n, layer, transform)

        end function box_init_scal

        real function eval_box(this, pos)

            implicit none

            class(box) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_box = box_fn(p, this%lengths)

        end function eval_box

        real function box_fn(p, b)

            implicit none

            type(vector), intent(IN) :: p, b

            type(vector) :: q

            q = abs(p) - b
            box_fn = length(max(q, 0.)) + min(max(q%x, max(q%y, q%z)), 0.)

        end function box_fn



        function sphere_init(radius, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(sphere) :: out
            
            real,            intent(IN) :: radius, mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%radius = radius
            out%layer = layer

            out%transform = t

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function sphere_init

        real function eval_sphere(this, pos)

            implicit none

            class(sphere) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_sphere = sphere_fn(p, this%radius)

        end function eval_sphere

        real function sphere_fn(p, c)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: c

            sphere_fn = sqrt(p%x**2+p%y**2+p%z**2) - c

        end function sphere_fn


        function torus_init(oradius, iradius, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(torus) :: out
            
            real,            intent(IN) :: oradius, iradius, mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%oradius = oradius
            out%iradius = iradius
            out%layer = layer
            out%transform = t

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function torus_init

        real function eval_torus(this, pos)

            implicit none

            class(torus) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_torus = torus_fn(p, this%oradius, this%iradius)

        end function eval_torus

        real function torus_fn(p, or, ir)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: or, ir


            type(vector) :: q

            q = vector(length(vector(p%x, 0., p%z)) - or, p%y, 0.)
            torus_fn = length(q) - ir

        end function torus_fn


        function triprisim_init(h1, h2, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(triprisim) :: out
            
            real,            intent(IN) :: h1, h2, mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%h1 = h1
            out%h2 = h2
            out%layer = layer
            out%transform = t

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function triprisim_init

        real function eval_triprisim(this, pos)

            implicit none

            class(triprisim) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_triprisim = triprisim_fn(p, this%h1, this%h2)

        end function eval_triprisim

        real function triprisim_fn(p, h1, h2)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: h1, h2


            type(vector) :: q

            q = abs(p)
            triprisim_fn = max(q%z - h2, max(q%x*.866025 + p%y*.5, -p%y) - h1*.5) 

        end function triprisim_fn

        function cone_init(a, b, ra, rb, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(cone) :: out
            
            type(vector),    intent(IN) :: a, b
            real,            intent(IN) :: ra, rb, mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1d-9)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function cone_init

        real function eval_cone(this, pos)

            implicit none

            class(cone) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_cone = cone_fn(p, this%a, this%b, this%ra, this%rb)

        end function eval_cone

        real function cone_fn(p, a, b, ra, rb)

            use utils, only : clamp

            implicit none

            type(vector), intent(IN) :: p, a, b
            real,         intent(IN) :: ra, rb


            type(vector) :: q
            real :: rba, baba, papa, paba, x, cax, cay, k, f, cbx, cby, s

            rba = rb - ra
            baba = (b-a) .dot. (b-a)
            papa = (p-a) .dot. (p-a)
            paba =  ((p-a) .dot. (b-a))/ baba
            x = sqrt(papa - baba*paba**2)
            if(paba < 0.5)then
                cax = max(0., x - ra)
            else
                cax = max(0., x - rb)
            end if
            cay = abs(paba - 0.5) - .5
            k = rba**2 + baba
            f = clamp((rba * (x - ra) + paba*baba) / k, 0., 1.)
            cbx = x - ra - f*rba
            cby = paba - f
            if(cbx < 0. .and. cay <0.)then
                s = -1.
            else
                s = 1.
            end if 
            cone_fn = s * sqrt(min(cax**2 + baba*cay**2, cbx**2 + baba*cby**2)) 

        end function cone_fn


        function translate(o) result(out)

            implicit none

            type(vector), intent(IN) :: o

            real :: out(4, 4)

            out(:, 1) = [1., 0., 0., o%x] 
            out(:, 2) = [0., 1., 0., o%y] 
            out(:, 3) = [0., 0., 1., o%z] 
            out(:, 4) = [0., 0., 0., 1.] 

        end function translate

        real function union(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            union = min(d1, d2)
        end function union


        real function SmoothUnion(d1, d2)

            use utils, only : mix, clamp

            implicit none

            real, intent(IN) :: d1, d2
            real :: k=0.01, h

            h = clamp(0.5 +.5*(d2-d1)/k, 0., 1.)
            SmoothUnion = mix(d2, d1, h) - k*h*(1.-h)

        end function SmoothUnion

        real function subtraction(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            subtraction = max(-d1, d2)

        end function subtraction

        real function intersection(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            intersection = max(d1, d2)

        end function intersection

        type(elongate) function elongate_init(prim, size) result(out)

            implicit none

            type(vector), intent(IN) :: size
            class(sdf), target :: prim

            out%size = size
            out%prim => prim

            out%mus = prim%mus
            out%mua = prim%mua
            out%hgg = prim%hgg
            out%g2 = prim%g2
            out%n = prim%n
            out%kappa = prim%kappa
            out%albedo = prim%kappa
            out%layer = prim%layer
            out%transform = identity()

            end function elongate_init

        real function eval_elongate(this, pos)

            implicit none

            class(elongate) :: this
            type(vector), intent(IN) :: pos

            eval_elongate = elongate_fn(pos, this%size, this%prim)

        end function eval_elongate

        real function elongate_fn(p, size, prim)

            implicit none

            class(sdf) :: prim

            type(vector), intent(IN) :: size
            type(vector), intent(IN) :: p

            real :: w
            type(vector) :: q

            q = abs(p) - size
            w = min(max(q%x, max(q%y, q%z)), 0.)

            elongate_fn = prim%evaluate(max(q, 0.)) + w

        end function elongate_fn

        type(bend) function bend_init(prim, k) result(out)

            implicit none

            real, intent(IN) :: k
            class(sdf), target :: prim

            out%k = k
            out%prim => prim

            out%mus = prim%mus
            out%mua = prim%mua
            out%hgg = prim%hgg
            out%g2 = prim%g2
            out%n = prim%n
            out%kappa = prim%kappa
            out%albedo = prim%kappa
            out%layer = prim%layer
            out%transform = identity()

            end function bend_init

        real function eval_bend(this, pos)

            implicit none

            class(bend) :: this
            type(vector), intent(IN) :: pos

            eval_bend = bend_fn(pos, this%k, this%prim)

        end function eval_bend

        real function bend_fn(p, k, prim)

            implicit none

            class(sdf) :: prim

            real, intent(IN)         :: k
            type(vector), intent(IN) :: p

            real :: c, s, x2, y2, z2

            c = cos(k * p%x)
            s = sin(k * p%x)
            x2 = c * p%x - s * p%y
            y2 = s * p%x + c * p%y
            z2 = p%z

            bend_fn = prim%evaluate(vector(x2, y2, z2))

        end function bend_fn

        type(displacement) function displacement_init(prim, func) result(out)

            implicit none

            class(sdf), target :: prim
            procedure(primitive) :: func

            out%func => func
            out%prim => prim

            out%mus = prim%mus
            out%mua = prim%mua
            out%hgg = prim%hgg
            out%g2 = prim%g2
            out%n = prim%n
            out%kappa = prim%kappa
            out%albedo = prim%kappa
            out%layer = prim%layer
            out%transform = identity()

            end function displacement_init

        real function eval_disp(this, pos)

            implicit none

            class(displacement) :: this
            type(vector), intent(IN) :: pos

            eval_disp = displacement_fn(pos, this%prim, this%func)


        end function eval_disp

        real function displacement_fn(p, prim, disp)

            implicit none

            class(sdf) :: prim
            procedure(primitive) :: disp
            type(vector), intent(IN) :: p

            real :: d1, d2

            d1 = prim%evaluate(p)
            d2 = disp(p)

            displacement_fn = d1 + d2

        end function displacement_fn

        type(twist) function twist_init(prim, k) result(out)

            implicit none

            class(sdf), target :: prim
            real :: k

            out%k = k
            out%prim => prim
            out%mus = prim%mus
            out%mua = prim%mua
            out%hgg = prim%hgg
            out%g2 = prim%g2
            out%n = prim%n
            out%kappa = prim%kappa
            out%albedo = prim%kappa
            out%layer = prim%layer
            out%transform = identity()

            end function twist_init

        real function eval_twist(this, pos)

            implicit none

            class(twist) :: this
            type(vector), intent(IN) :: pos

            eval_twist = twist_fn(pos, this%k, this%prim)

        end function eval_twist

        real function twist_fn(p, k, prim)

            implicit none

            class(sdf) :: prim
            type(vector) :: p
            real :: k

            real :: c, s, x2, y2, z2, x, y, z

            c = cos(k * p%z)
            s = sin(k * p%z)
            x2 = c*p%x - s*p%y
            y2 = s*p%x + c*p%y
            z2 = p%z

            twist_fn = prim%evaluate(vector(x2, y2, z2))

        end function twist_fn


        function rotate_x(angle) result(r)
        ! rotation funcs from https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/
            use utils, only : deg2rad

            implicit none
            
            real, intent(IN) :: angle
            real :: r(4, 4), c, s, a

            a = deg2rad(angle)
            c = cos(a)
            s = sin(a)

            r(:, 1) = [1., 0., 0., 0.]
            r(:, 2) = [0., c,  s,  0.]
            r(:, 3) = [0.,-s,  c,  0.]
            r(:, 4) = [0., 0., 0., 1.]

        end function rotate_x

        function rotate_y(angle) result(r)
            
            use utils, only : deg2rad

            implicit none
            
            real, intent(IN) :: angle
            real :: r(4, 4), c, s, a

            a = deg2rad(angle)
            ! c = cos(a)
            s = sin(a)
            c = sqrt(1. - s**2)

            r(:, 1) = [c,  0., -s,  0.]
            r(:, 2) = [0., 1.,  0., 0.]
            r(:, 3) = [s,  0.,  c,  0.]
            r(:, 4) = [0., 0.,  0., 1.]

        end function rotate_y

        function rotate_z(angle) result(r)
            
            use utils, only : deg2rad

            implicit none
            
            real, intent(IN) :: angle
            real :: r(4, 4), c, s, a

            a = deg2rad(angle)
            c = cos(a)
            s = sin(a)

            r(:, 1) = [c, -s,  0., 0.]
            r(:, 2) = [s,  c,  0., 0.]
            r(:, 3) = [0., 0., 1., 0.]
            r(:, 4) = [0., 0., 0., 1.]

        end function rotate_z

        function identity() result(r)
            
            implicit none
            
            real :: r(4, 4)

            r(:, 1) = [1., 0., 0., 0.]
            r(:, 2) = [0., 1., 0., 0.]
            r(:, 3) = [0., 0., 1., 0.]
            r(:, 4) = [0., 0., 0., 1.]

        end function identity


        subroutine render(cnt, extent, samples, fname)

            implicit none
            
            type(container),        intent(IN) :: cnt(:)
            integer,                intent(IN) :: samples
            type(vector),           intent(IN) :: extent
            character(*), optional, intent(IN) :: fname

            type(vector)      :: pos, wid
            integer           :: i, j, k, u, ns
            real              :: x, y, z, ds(size(cnt))
            real, allocatable :: image(:, :, :)
            
            character(len=:), allocatable  :: filename

            if(present(fname))then
                filename = fname
            else
                filename = "../model.dat"
            end if

            ns = int(samples / 2)
            allocate(image(samples, samples, samples))
            wid = extent/real(ns)
!$omp parallel default(none) shared(cnt, ns, wid, image, samples)&
!$omp private(i, x, y, z, pos, j, k, u, ds)
!$omp do
            do i = 1, samples
                x = (i-ns) *wid%x
                do j = 1, samples
                    y = (j-ns) *wid%y
                    do k = 1, samples
                        z = (k-ns) * wid%z
                        pos = vector(x, y, z)
                        ds = 0.
                        do u = 1, size(ds)
                            ds(u) = cnt(u)%p%evaluate(pos)
                        end do

                        image(i, j, k) = minval(ds)
                    end do
                end do
            end do
!$OMP end  do
!$OMP end parallel
            open(newunit=u,file=filename, access="stream", form="unformatted")
            write(u)image
            close(u)
        end subroutine render
end module sdfs