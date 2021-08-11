module sdfs

    use vector_class

    implicit none

    type :: sdf
        real :: mus, mua, kappa, albedo, hgg, g2, n
        real :: transform(4, 4)
        type(vector) :: centre
        integer :: layer
        contains
        procedure :: evaluate => evaluate_fn
    end type sdf

    type :: container
        class(sdf), pointer :: p => null()
    end type container

    type, extends(sdf) :: sphere
        real         :: radius
        contains
        procedure :: evaluate => eval_sphere
    end type sphere

    type, extends(sdf) :: cylinder
        real         :: radius
        type(vector) :: a, b
        contains
        procedure :: evaluate => eval_cylinder
    end type cylinder

    type, extends(sdf) :: box
        type(vector) :: lengths
        contains
        procedure :: evaluate => eval_box
    end type box

    type, extends(sdf) :: torus
        real :: oradius, iradius
        contains
        procedure :: evaluate => eval_torus
    end type torus

    type, extends(sdf) :: model
        type(container), allocatable  :: array(:)
        procedure(op),nopass, pointer :: func
        contains
        procedure :: evaluate => eval_model
    end type model

    interface cylinder
        module procedure cylinder_init
    end interface cylinder


    interface sphere
        module procedure sphere_init
    end interface sphere

    interface torus
        module procedure torus_init
    end interface torus

    interface box
        module procedure box_init_vec
        module procedure box_init_scal
    end interface box

    interface model
        module procedure model_init
    end interface model

    abstract interface
        real function op(d1, d2)
            implicit none
            real, intent(IN) :: d1, d2
        end function op
    end interface


    private
    public :: sdf, model, cylinder, sphere, box, container, model_init, render, torus
    public :: union, intersection, subtraction, calcNormal, onion, SmoothUnion, op
    public :: rotate_x, rotate_y, rotate_z, identity, translate

    contains

        type(vector) function calcNormal(p, obj)

            implicit none

            type(vector), intent(IN) :: p
            class(sdf) :: obj

            real :: h
            type(vector) :: xyy, yyx, yxy, xxx

            h = 1e-8
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

        real function evaluate_fn(this, pos)
            
            implicit none
            class(sdf) :: this
            type(vector), intent(IN) :: pos
        end function evaluate_fn

        function cylinder_init(a, b, radius, mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(cylinder) :: out
            
            real,         intent(IN) :: radius, mus, mua, hgg, n
            integer,      intent(IN) :: layer
            type(vector), intent(IN) :: a, b

            real,         optional :: transform(4, 4)

            real         :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%a = a
            out%b = b
            out%radius = radius
            out%centre = vector(0.,0.,0.)
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

        function box_init(lengths, mus, mua, hgg, n, layer, c, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            type(vector), intent(IN) :: lengths
            real,         intent(IN) :: mus, mua, hgg, n
            integer,      intent(IN) :: layer

            real,         optional, intent(in) :: transform(4, 4)
            type(vector), optional, intent(IN) :: c

            real :: t(4, 4)
            type(vector) :: centre

            if(present(c))then
                centre = c
            else
                centre = vector(0., 0., 0.)
            end if

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if
            out%lengths = .5*lengths! as only half lengths
            out%layer = layer
            out%centre = centre
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

        function box_init_vec(lengths, mus, mua, hgg, n, layer, c, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            type(vector), intent(IN) :: lengths
            real,         intent(IN) :: mus, mua, hgg, n
            integer,      intent(IN) :: layer

            real,         optional, intent(in) :: transform(4, 4)
            type(vector), optional, intent(IN) :: c

            out = box_init(lengths, mus, mua, hgg, n, layer, c, transform)

        end function box_init_vec


        function box_init_scal(length, mus, mua, hgg, n, layer, c, transform) result(out)
        
            implicit none
        
            type(box) :: out
            
            real,         intent(IN) :: length, mus, mua, hgg, n
            integer,      intent(IN) :: layer

            real,         optional, intent(in) :: transform(4, 4)
            type(vector), optional, intent(IN) :: c

            type(vector) :: lengths

            lengths = vector(length, length, length)

            out = box_init(lengths, mus, mua, hgg, n, layer, c, transform)

        end function box_init_scal


        function sphere_init(radius, mus, mua, hgg, n, layer, c, transform) result(out)
        
            implicit none
        
            type(sphere) :: out
            
            real, intent(IN) :: radius, mus, mua, hgg, n
            integer, intent(IN) :: layer
            type(vector), optional, intent(IN) :: c
            real,  optional, intent(IN) :: transform

            type(vector) :: centre
            real         :: t(4, 4)

            if(present(c))then
                centre = c
            else
                centre = vector(0., 0., 0.)
            end if

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%radius = radius
            out%layer = layer
            out%centre = centre

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


        function torus_init(oradius, iradius, mus, mua, hgg, n, layer, c, transform) result(out)
                             ! .25,    0.,  100.,0.,  1.,1,     c=vector(0., 0., -0.75)
            implicit none
        
            type(torus) :: out
            
            real, intent(IN) :: oradius, iradius, mus, mua, hgg, n
            integer, intent(IN) :: layer
            type(vector), optional, intent(IN) :: c
            real,  optional, intent(IN) :: transform

            type(vector) :: centre
            real         :: t(4, 4)

            if(present(c))then
                centre = c
            else
                centre = vector(0., 0., 0.)
            end if

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%oradius = oradius
            out%iradius = iradius
            out%layer = layer
            out%centre = centre

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

        real function eval_cylinder(this, pos)

            implicit none

            class(cylinder) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_cylinder = cylinder_fn(p, this%a, this%b, this%radius)

        end function eval_cylinder


        real function eval_sphere(this, pos)

            implicit none

            class(sphere) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_sphere = sphere_fn(p, this%radius)

        end function eval_sphere


        real function eval_box(this, pos)

            implicit none

            class(box) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_box = box_fn(p, this%lengths)

        end function eval_box


        real function eval_torus(this, pos)

            implicit none

            class(torus) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_torus = torus_fn(p, this%oradius, this%iradius)

        end function eval_torus

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

            !lies along y plane
            ! d = abs(vector(sqrt(p%x**2 + p%z**2), p%y, 0.)) - vector(h, r, 0.)
            ! cylinder_fn = min(max(d%x, d%y), 0.) + length(max(d, 0.))

        end function cylinder_fn


        real function sphere_fn(p, c)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: c

            sphere_fn = sqrt(p%x**2+p%y**2+p%z**2) - c

        end function sphere_fn


        real function box_fn(p, b)

            implicit none

            type(vector), intent(IN) :: p, b

            type(vector) :: q

            q = abs(p) - b
            box_fn = length(max(q, 0.)) + min(max(q%x, max(q%y, q%z)), 0.)

        end function box_fn

        real function torus_fn(p, or, ir)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: or, ir


            type(vector) :: q

            q = vector(length(vector(p%x, 0., p%z)) - or, p%y, 0.)
            torus_fn = length(q) - ir

        end function torus_fn

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
            real :: k=0.1, h

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

        real function extrude(p, h, prim)

            implicit none

            type(vector), intent(IN) :: p
            real, intent(IN) :: h
            procedure(primitive) :: prim

            real :: d
            type(vector) :: w

            d = prim(vector(p%x, p%y, 0.))
            w = vector(d, abs(p%z) - h, 0.)
            extrude = min(max(w%x, w%y), 0.) + length(max(w, 0.))

        end function extrude


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

        subroutine render(cnt, extent, samples)

            implicit none
            
            type(container), intent(IN) :: cnt(:)
            integer,         intent(IN) :: samples
            type(vector),    intent(IN) :: extent

            type(vector)      :: pos, wid
            integer           :: i, j, k, u, ns
            real              :: x, y, z, ds(size(cnt))
            real, allocatable :: image(:, :, :)

            ns = int(samples / 2)
            allocate(image(-ns:ns, -ns:ns, -ns:ns))
            wid = extent/100.
!$omp parallel default(none) shared(cnt, ns, wid, image)&
!$omp private(i, x, y, z, pos, j, k, u, ds)
!$omp do
            do i = -ns, ns
                x = i *wid%x
                do j = -ns, ns
                    y = j *wid%y
                    do k = -ns, ns
                        z = k * wid%z
                        pos = vector(x, y, z)
                        ds = 0.
                        do u = 1, size(ds)
                            ds(u) = cnt(u)%p%evaluate(pos)
                        end do
                        image(i, j, k) = minval(abs(ds))
                    end do
                end do
            end do
!$OMP end  do
!$OMP end parallel
            open(newunit=u,file="../model.dat", access="stream", form="unformatted")
            write(u)image
            close(u)
        end subroutine render
end module sdfs