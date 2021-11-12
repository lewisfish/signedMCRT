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

    type, extends(sdf) :: capsule
        type(vector) :: a, b
        real         :: r
        contains
        procedure :: evaluate => eval_capsule
    end type capsule

    interface capsule
        module procedure capsule_init
    end interface capsule

    type, extends(sdf) :: plane
        type(vector) :: a
        contains
        procedure :: evaluate => eval_plane
    end type plane

    interface plane
        module procedure plane_init
    end interface plane

    type, extends(sdf) :: moon
        real :: d, ra, rb
        contains
        procedure :: evaluate => eval_moon
    end type moon

    interface moon
        module procedure moon_init
    end interface moon

    type, extends(sdf) :: neural
        contains
        procedure :: evaluate => eval_neural
    end type neural

    interface neural
        module procedure neural_init
    end interface neural


    type, extends(sdf) :: model
        type(container), allocatable  :: array(:)
        procedure(op), nopass, pointer :: func
        real :: k
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

    type, extends(sdf) :: repeat
        real :: c
        type(vector) :: la, lb
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_repeat
    end type repeat

    interface repeat
        module procedure repeat_init
    end interface repeat

    type, extends(sdf) :: extrude
        real :: h
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_extrude
    end type extrude

    interface extrude
        module procedure extrude_init
    end interface extrude



    abstract interface
        real function op(d1, d2, k)
            implicit none
            real, intent(IN) :: d1, d2, k
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
    public :: capsule, plane, moon, neural!, text
    ! meta
    public :: model, container
    ! boolean ops
    public :: union, intersection, subtraction, SmoothUnion, op
    ! move ops
    public :: rotate_x, rotate_y, rotate_z, identity, translate
    ! deform ops
    public :: displacement, bend, twist, elongate, repeat, extrude
    ! utility funcs
    public :: calcNormal, model_init, render

    contains

        function neural_init(mus, mua, hgg, n, layer, transform) result(out)
        
            implicit none
        
            type(neural) :: out
            
            real,           intent(IN) :: mus, mua, hgg, n
            integer,        intent(IN) :: layer
            real, optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

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

        end function neural_init



real function eval_neural(this, pos) result(out)

    use vec4_class
    use mat_class

    implicit none

    class(neural) :: this
    type(vector), intent(IN) :: pos

    type(vec4) :: f0_0, f0_1, f0_2, f0_3, f1_0, f1_1, f1_2, f1_3, f2_0, f2_1, f2_2, f2_3
 
f0_0=sin(pos%y*vec4(-2.845,2.169,-3.165,-3.379)+pos%z*vec4(-1.938,2.632,-1.872,3.019)+pos%x*vec4(1.042,4.283,-.778,-3.063)+vec4(8.246,7.046,6.947,4.930))
f0_1=sin(pos%y*vec4(-.870,-.224,-3.854,-2.134)+pos%z*vec4(-1.474,-2.259,4.302,2.238)+pos%x*vec4(1.576,-.826,1.014,2.511)+vec4(.639,-.425,6.931,7.475))
f0_2=sin(pos%y*vec4(-2.161,1.539,-3.024,1.739)+pos%z*vec4(-3.388,.743,-3.374,-.536)+pos%x*vec4(3.843,-3.060,3.127,1.533)+vec4(1.916,-2.519,-.403,-1.918))
f0_3=sin(pos%y*vec4(.076,1.427,-2.403,-.442)+pos%z*vec4(3.012,1.584,-3.421,1.513)+pos%x*vec4(.061,-3.592,-.935,.125)+vec4(7.654,-.949,-7.044,4.319))
f1_0=sin(mat([.230,.244,-.196,-.230,-.303,.187,.122,-.410,.388,-.741,.731,.031,-.999,-.044,.107,.171])*f0_0+&
   mat([.141,-.542,-.726,-.280,-.445,-.038,-.397,.348,-.086,.413,.544,.298,-.075,.122,-.315,.293])*f0_1+&
   mat([.329,-.106,-.690,.012,-.073,-.862,.381,.628,.378,-.012,-.137,-.071,.580,.060,-1.094,.652])*f0_2+&
   mat([.616,.363,.211,-.366,-.396,.388,-.090,.217,.479,-.339,-.039,-.132,.176,.223,.312,.082])*f0_3+&
   vec4(-2.491,3.220,2.193,2.989))/1.0+f0_0
f1_1=sin(mat([.137,.498,.156,.391,-.059,.084,-.182,-.027,.650,.302,.676,-.202,-.488,.297,-.136,-.366])*f0_0+&
   mat([.211,.210,-.557,.021,.808,-.041,-.726,.028,.003,-.012,-.280,-.539,-.624,.261,.262,-.123])*f0_1+&
   mat([-.220,-.208,.211,-.740,-.009,-.103,.741,-.820,-.312,.074,.432,-.258,.390,-.249,.375,.169])*f0_2+&
   mat([-.085,-.105,.309,.291,-.476,-.764,.048,-.148,-.318,.182,-.170,-.220,-1.041,-.116,-1.077,.565])*f0_3+&
   vec4(1.529,2.058,-.449,-1.700))/1.0+f0_1
f1_2=sin(mat([.385,.698,-.630,-.482,-.099,-.552,.173,-.148,-.711,-.721,.342,.079,-.470,-.177,.011,-.117])*f0_0+&
   mat([.834,-.758,.162,.761,-.348,.367,.375,.721,.731,.078,.024,.329,.340,-.429,-.006,.084])*f0_1+&
   mat([.277,.015,-.449,-.100,-.223,-.897,.623,-.004,.049,-.371,.040,-.200,.367,.221,-.063,.665])*f0_2+&
   mat([1.595,.428,.417,.268,-.241,.145,.141,-.185,.430,-.096,.521,.109,-.239,.770,-.100,-.128])*f0_3+&
   vec4(-1.674,2.261,-2.165,-3.244))/1.0+f0_2
f1_3=sin(mat([.228,-.315,.459,.085,-.060,-.070,-.217,-.605,-.217,-.237,.400,.358,.198,-.722,-.062,-.805])*f0_0+&
   mat([.108,.313,.793,-.378,-.168,-.350,.504,.664,.251,-.058,.087,-.405,.068,.077,.148,-.741])*f0_1+&
   mat([.324,-.130,.018,.100,-.658,.087,-.077,-.082,.101,-.103,-.132,.159,-.117,-.563,-.606,-.278])*f0_2+&
   mat([.181,.584,.390,-.498,.096,-.461,-.324,.573,-.261,.849,-.012,.017,-.289,.578,.161,-1.158])*f0_3+&
   vec4(-.938,-.744,2.675,1.662))/1.0+f0_3
f2_0=sin(mat([-.165,-.648,-.274,-.443,.307,-.478,-.204,-.103,-.017,1.027,.680,.228,-.121,.419,.137,-.253])*f1_0+&
   mat([-.252,.249,.450,-.249,-1.085,-.292,.609,.047,-.999,-.082,-.005,-.836,-.356,.411,-.402,.699])*f1_1+&
   mat([-.103,-.218,-.000,-.639,-.057,.447,-.502,-.505,-.298,.029,-.363,-.607,.116,-.942,2.064,.363])*f1_2+&
   mat([.394,.030,.070,1.402,-.552,-.391,.255,-.733,-.102,.263,-.291,.190,.489,.614,.683,.061])*f1_3+&
   vec4(-.896,.069,3.048,2.889))/1.4+f1_0
f2_1=sin(mat([.553,-.627,-.284,-1.032,-.524,.294,-.078,-.128,.100,.112,.089,.318,.560,.407,-.317,-1.634])*f1_0+&
   mat([.024,.420,.599,.063,-.317,.728,.197,-1.037,.425,-.035,-.240,-1.530,-.246,.333,.209,-1.455])*f1_1+&
   mat([.514,-.095,.139,-.551,-.225,-.067,.568,-.958,.038,.161,-.057,-.871,.468,-.633,.440,1.798])*f1_2+&
   mat([.185,-.798,.206,.941,1.032,-.016,.629,-.416,.026,-.777,.206,-.360,-.456,.463,1.117,.883])*f1_3+&
   vec4(2.594,2.602,-4.067,.754))/1.4+f1_1
f2_2=sin(mat([.601,.096,-.071,.363,-1.159,-.053,.049,-.140,.056,.701,-.047,-.031,1.257,-.421,-.360,-1.596])*f1_0+&
   mat([.357,.525,-.130,-.001,-.643,-.206,.237,-.837,.813,-.308,.407,-.157,1.845,.737,-.374,.390])*f1_1+&
   mat([-.569,.109,.476,-.377,.124,.171,-.100,-.606,-.560,.414,.176,.483,-2.082,-.375,.928,.125])*f1_2+&
   mat([-.599,.549,-.414,.642,-.525,.025,.385,-.797,1.270,.342,-.068,.965,-1.026,-.177,.890,.156])*f1_3+&
   vec4(-4.226,1.935,3.830,-2.055))/1.4+f1_2
f2_3=sin(mat([.463,1.053,.008,-.354,-.124,-.566,-.734,.673,.970,-.193,.377,.283,-.481,-1.345,.024,.687])*f1_0+&
   mat([-.358,-.069,.358,.373,.745,-.056,.054,1.070,-.459,-.404,-.654,.665,-.400,.127,.422,.202])*f1_1+&
   mat([.255,-.513,-.741,-.345,-.404,.188,.118,.486,.783,.610,-.790,-.099,-.555,-3.111,-1.241,-.215])*f1_2+&
   mat([.776,.131,.035,.434,.355,-2.543,-.158,-.140,.454,1.154,.156,.543,.555,-.726,.369,-1.032])*f1_3+&
   vec4(1.108,.703,-.637,.098))/1.4+f1_3
out= (f2_0.dot.vec4(-.041,-.048,-.054,.034))+&
    (f2_1.dot.vec4(.066,.103,.066,-.022))+&
    (f2_2.dot.vec4(.011,.059,-.092,.064))+&
    (f2_3.dot.vec4(-.059,-.031,.063,-.052))+&
    (0.079)
end function eval_neural







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

        function model_init(array, func, kopt) result(out)
        !TODO make sure optical properties are same in all inputs            
            implicit none

            type(model) :: out

            procedure(op) :: func
            type(container), intent(IN) :: array(:)
            real, optional,  intent(IN) :: kopt

            out%array = array
            out%func => func
            if(present(kopt))then
                out%k = kopt
            else
                out%k = 0.
            end if

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
                eval_model = this%func(eval_model, this%array(i)%p%evaluate(pos), this%k)
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
        !h1 is height
        !h2 is length
        !
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

        function capsule_init(a, b, r, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(capsule) :: out
            
            type(vector),    intent(IN) :: a, b
            real,            intent(IN) :: r, mus, mua, hgg, n
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
            out%r = r
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

        end function capsule_init

        real function eval_capsule(this, pos)

            implicit none

            class(capsule) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_capsule = capsule_fn(p, this%a, this%b, this%r)

        end function eval_capsule

        real function capsule_fn(p, a, b, r)

            use utils, only : clamp

            implicit none

            type(vector), intent(IN) :: p, a, b
            real,         intent(IN) :: r

            type(vector) :: pa, ba
            real :: h

            pa = p - a
            ba = b - a
            h = clamp((pa .dot. ba) / (ba .dot. ba), 0., 1.)
            capsule_fn = length(pa - ba*h) - r

        end function capsule_fn


        function plane_init(a, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(plane) :: out
            
            type(vector),    intent(IN) :: a
            real,            intent(IN) :: mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%a = a
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

        end function plane_init

        real function eval_plane(this, pos)

            implicit none

            class(plane) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_plane = plane_fn(p, this%a)

        end function eval_plane

        real function plane_fn(p, n)

            use utils, only : clamp

            implicit none

            type(vector), intent(IN) :: p, n

            !n must be normalised
            plane_fn = (p .dot. n)

        end function plane_fn


        function moon_init(d, ra, rb, mus, mua, hgg, n, layer, transform) result(out)

            implicit none
        
            type(moon) :: out
            
            real,            intent(IN) :: d, ra, rb, mus, mua, hgg, n
            integer,         intent(IN) :: layer
            real,  optional, intent(IN) :: transform(4, 4)

            real :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%d = d
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

        end function moon_init

        real function eval_moon(this, pos)

            implicit none

            class(moon) :: this
            type(vector), intent(IN) :: pos

            type(vector) :: p

            p = pos .dot. this%transform
            eval_moon = moon_fn(p, this%d, this%ra, this%rb)

        end function eval_moon

        real function moon_fn(p, d, ra, rb)

            use utils, only : clamp

            implicit none

            type(vector), intent(IN) :: p
            real, intent(IN) :: d, ra, rb

            type(vector) :: pos
            real :: a, b, ra2, rb2, d2

            pos = vector(p%x, abs(p%y), 0.)
            ra2 = ra*ra
            rb2 = rb*rb
            d2 = d*d
            a = (ra2 - rb2 + d2) / (2.*d)
            b = sqrt(max(ra2 - a**2, 0.))
            if(d*(pos%x*b - pos%y*a) > d2*max(b - pos%y, 0.))then
                moon_fn = length(pos - vector(a, b, 0.))
            else
                moon_fn = max(-length(pos) - ra, length(pos - vector(d, 0., 0.)) - rb)
            end if

        end function moon_fn



        function translate(o) result(out)

            implicit none

            type(vector), intent(IN) :: o

            real :: out(4, 4)

            out(:, 1) = [1., 0., 0., o%x] 
            out(:, 2) = [0., 1., 0., o%y] 
            out(:, 3) = [0., 0., 1., o%z] 
            out(:, 4) = [0., 0., 0., 1.] 

        end function translate

        real function union(d1, d2, k)

            implicit none

            real, intent(IN) :: d1, d2, k

            union = min(d1, d2)
        end function union


        real function SmoothUnion(d1, d2, k)

            use utils, only : mix, clamp

            implicit none

            real, intent(IN) :: d1, d2, k
            real :: h

            h = max(k - abs(d1 - d2), 0.) / k
            SmoothUnion = min(d1, d2) - h*h*h*k*(1./6.)
            ! h = clamp(0.5 +.5*(d2-d1)/k, 0., 1.)
            ! SmoothUnion = mix(d2, d1, h) - k*h*(1.-h)

        end function SmoothUnion

        real function subtraction(d1, d2, k)

            implicit none

            real, intent(IN) :: d1, d2, k

            subtraction = max(-d1, d2)

        end function subtraction

        real function intersection(d1, d2, k)

            implicit none

            real, intent(IN) :: d1, d2, k

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
            type(vector), intent(IN) :: p
            real,         intent(IN) :: k

            real :: c, s, x2, y2, z2

            c = cos(k * p%z)
            s = sin(k * p%z)
            x2 = c*p%x - s*p%y
            y2 = s*p%x + c*p%y
            z2 = p%z

            twist_fn = prim%evaluate(vector(x2, y2, z2))

        end function twist_fn


        type(repeat) function repeat_init(prim, c, la, lb) result(out)

            implicit none

            class(sdf), target :: prim
            type(vector), intent(IN) :: la, lb
            real,         intent(IN) :: c

            out%c = c
            out%la = la
            out%lb = lb
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

            end function repeat_init

        real function eval_repeat(this, pos)

            implicit none

            class(repeat) :: this
            type(vector), intent(IN) :: pos

            eval_repeat = repeat_fn(pos, this%c, this%la, this%lb, this%prim)

        end function eval_repeat

        real function repeat_fn(p, c, la, lb, prim)

            use vector_class

            implicit none

            class(sdf) :: prim
            type(vector), intent(IN) :: p, la, lb
            real,         intent(IN) :: c

            type(vector) :: q

            q = p - c*clamp_vec(nint(p/c), la, lb)
            repeat_fn = prim%evaluate(q)

        end function repeat_fn

        type(extrude) function extrude_init(prim, h) result(out)

            implicit none

            class(sdf), target :: prim
            real, intent(IN)   :: h

            out%h = h
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

        end function extrude_init

        real function eval_extrude(this, pos)

            implicit none

            class(extrude) :: this
            type(vector), intent(IN) :: pos

            eval_extrude = extrude_fn(pos, this%h, this%prim)

        end function eval_extrude


        real function extrude_fn(p, h, prim)

            implicit none

            class(sdf) :: prim
            type(vector), intent(IN) :: p
            real,         intent(IN) ::  h

            type(vector) :: w
            real :: d

            d = prim%evaluate(p)
            w = vector(d, p%z - h, 0.)
            extrude_fn = min(max(w%x, w%y), 0.) + length(max(w, 0.))

        end function extrude_fn



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
            c = cos(a)
            s = sin(a)

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

            use utils, only : pbar

            implicit none
            
            type(container),        intent(IN) :: cnt(:)
            integer,                intent(IN) :: samples
            type(vector),           intent(IN) :: extent
            character(*), optional, intent(IN) :: fname

            type(vector)      :: pos, wid
            integer           :: i, j, k, u, ns, id
            real              :: x, y, z, ds(size(cnt)-1)
            real, allocatable :: image(:, :, :)
            type(pbar)        :: bar

            character(len=:), allocatable  :: filename

            if(present(fname))then
                filename = fname
            else
                filename = "model.dat"
            end if
            ns = int(samples / 2)
            allocate(image(samples, samples, samples))
            wid = extent/real(ns)
            bar = pbar(samples)
!$omp parallel default(none) shared(cnt, ns, wid, image, samples, bar)&
!$omp private(i, x, y, z, pos, j, k, u, ds, id)
!$omp do
            do i = 1, samples
                call bar%progress()
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

                        if(all(ds > 0.))then
                            id=0.
                        else
                            if(maxval(ds) < 0.)then
                                id = cnt(maxloc(ds,dim=1))%p%layer
                            else
                                id = cnt(minloc(ds,dim=1))%p%layer
                            end if
                        end if
                        if(minval(ds) > 0)then
                            id = 0
                        else
                            id =  minval(ds)!minloc(ds, dim=1)
                        end if
                        image(i, j, k) = minval(ds)!id
                    end do
                end do
            end do
!$OMP end  do
!$OMP end parallel
            open(newunit=u,file=filename, access="stream", form="unformatted", status="replace")
            write(u)image
            close(u)
        end subroutine render
end module sdfs