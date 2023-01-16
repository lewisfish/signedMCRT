module sdfs
!   Module provides signed distance functions (SDFs) for various shapes 
!   and some operations to adjust, move, rotate etc them
!   All SDF functions are adapted from Inigo Quilex exhaustive list at:
!   https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
!
    use vector_class
    use constants, only : wp

    implicit none

    type, abstract :: sdf
    ! base class.
    ! provides deferred evaluation func to get distance to shape
    ! and gives the transform and optical prperties.
    ! Layer is an important bookeeping integer.
        real(kind=wp) :: mus, mua, kappa, albedo, hgg, g2, n
        real(kind=wp) :: transform(4, 4)
        integer :: layer
        contains
        procedure(evalInterface), deferred :: evaluate
    end type sdf

    abstract interface
        function evalInterface(this, pos) result(res)
            use vector_class
            use constants, only : wp
            import sdf
            implicit none
            class(sdf) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res
        end function evalInterface
    end interface



    type, extends(sdf) :: sphere
        real(kind=wp) :: radius
        contains
        procedure :: evaluate => eval_sphere
    end type sphere

    interface sphere
        module procedure sphere_init
    end interface sphere


    type, extends(sdf) :: cylinder
        real(kind=wp) :: radius
        type(vector)  :: a, b
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


    type, extends(sdf) :: bevel_box
        type(vector)  :: size
        real(kind=wp) :: box_r
        contains
        procedure :: evaluate => eval_bevel_box
    end type bevel_box

    interface bevel_box
        module procedure bevel_box_init
    end interface bevel_box


    type, extends(sdf) :: torus
        real(kind=wp) :: oradius, iradius
        contains
        procedure :: evaluate => eval_torus
    end type torus

    interface torus
        module procedure torus_init
    end interface torus


    type, extends(sdf) :: triprism
        real(kind=wp) :: h1, h2
        contains
        procedure :: evaluate => eval_triprism
    end type triprism

    interface triprism
        module procedure triprism_init
    end interface triprism


    type, extends(sdf) :: cone
        type(vector)  :: a, b
        real(kind=wp) :: ra, rb
        contains
        procedure :: evaluate => eval_cone
    end type cone

    interface cone
        module procedure cone_init
    end interface cone

    type, extends(sdf) :: capsule
        type(vector)  :: a, b
        real(kind=wp) :: r
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

    type, extends(sdf) :: segment
        type(vector) :: a, b
        contains
        procedure :: evaluate => eval_segment
    end type segment

    interface segment
        module procedure segment_init
    end interface segment

    type, extends(sdf) :: egg
        real(kind=wp) :: r1, r2, h
        contains
        procedure :: evaluate => eval_egg
    end type egg

    interface egg
        module procedure egg_init
    end interface egg

    type, extends(sdf) :: moon
        real(kind=wp) :: d, ra, rb
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

    type :: container
        class(sdf), pointer :: p => null()
    end type container

    type, extends(sdf) :: model
        type(container), allocatable   :: array(:)
        procedure(op), nopass, pointer :: func
        real(kind=wp) :: k
        contains
        procedure :: evaluate => eval_model
    end type model

    interface model
        module procedure model_init
    end interface model





    type, extends(sdf) :: twist
        real(kind=wp) :: k
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
        real(kind=wp) :: k
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
        real(kind=wp) :: c
        type(vector) :: la, lb
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_repeat
    end type repeat

    interface repeat
        module procedure repeat_init
    end interface repeat

    type, extends(sdf) :: extrude
        real(kind=wp) :: h
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_extrude
    end type extrude

    interface extrude
        module procedure extrude_init
    end interface extrude

    type, extends(sdf) :: revolution
        real(kind=wp) :: o
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_revolution
    end type revolution

    interface revolution
        module procedure revolution_init
    end interface revolution

    type, extends(sdf) :: onion
        real(kind=wp) :: thickness
        class(sdf), pointer :: prim
        contains
        procedure :: evaluate => eval_onion
    end type onion

    interface onion
        module procedure onion_init
    end interface onion

    interface render
        module procedure render_scaler
        module procedure render_vec
    end interface render



    abstract interface
        function op(d1, d2, k) result(res)
            use constants, only : wp
            implicit none
            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res
        end function op
    end interface

    abstract interface
        function primitive(pos) result(res)
            use vector_class, only : vector
            use constants,    only : wp
            implicit none
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res
        end function primitive
    end interface

    private
    ! shapes
    public :: sdf, cylinder, sphere, box, torus, cone, triprism
    public :: capsule, plane, moon, segment, egg, neural
    ! meta
    public :: model, container
    ! boolean ops
    public :: union, intersection, subtraction, SmoothUnion, op
    ! move ops
    public :: rotate_x, rotate_y, rotate_z, identity, translate, rotationAlign
    ! deform ops
    public :: displacement, bend, twist, elongate, repeat, extrude, revolution, onion
    ! utility funcs
    public :: calcNormal, model_init, render, skewSymm, rotmat

    contains


        function segment_init(a, b, mus, mua, hgg, n, layer, transform) result(out)

            type(segment) :: out

            type(vector),            intent(IN) :: a, b
            real(kind=wp),           intent(IN) :: mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1.
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function segment_init

        function eval_segment(this, pos) result(res)

            class(segment) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = segment_fn(p, this%a, this%b) - 0.1_wp

        end function eval_segment

        function segment_fn(p, a, b) result(res)
            !p = pos
            !a = pt1
            !b = pt2
            !draws segment along the axis between 2 points a and b

            use utils, only : clamp

            type(vector), intent(IN) :: p, a, b
            real(kind=wp) :: res

            type(vector)  :: pa, ba
            real(kind=wp) :: h
           
            pa = p - a
            ba = b - a
            h = clamp((pa .dot. ba) / (ba .dot. ba), 0.0_wp, 1.0_wp)

            res = length(pa - ba*h)

        end function segment_fn


        function neural_init(mus, mua, hgg, n, layer, transform) result(out)
                
            type(neural) :: out
            
            real(kind=wp),           intent(IN) :: mus, mua, hgg, n
            integer,                 intent(IN) :: layer
            real(kind=wp), optional, intent(IN) :: transform(4, 4)

            real(kind=wp) :: t(4, 4)

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
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function neural_init



function eval_neural(this, pos) result(out)
! this function is generated using code. do not edit
    use vec4_class
    use mat_class

    class(neural) :: this
    type(vector), intent(IN) :: pos
    real(kind=wp) :: out

    type(vec4) :: f0_0, f0_1, f0_2, f0_3, f1_0, f1_1, f1_2, f1_3, f2_0, f2_1, f2_2, f2_3
 
f0_0=sin(pos%y*vec4(-2.845,2.169,-3.165,-3.379)+pos%z*vec4(-1.938,2.632,-1.872,3.019)+&
         pos%x*vec4(1.042,4.283,-.778,-3.063)+vec4(8.246,7.046,6.947,4.930))
f0_1=sin(pos%y*vec4(-.870,-.224,-3.854,-2.134)+pos%z*vec4(-1.474,-2.259,4.302,2.238)+&
         pos%x*vec4(1.576,-.826,1.014,2.511)+vec4(.639,-.425,6.931,7.475))
f0_2=sin(pos%y*vec4(-2.161,1.539,-3.024,1.739)+pos%z*vec4(-3.388,.743,-3.374,-.536)+&
         pos%x*vec4(3.843,-3.060,3.127,1.533)+vec4(1.916,-2.519,-.403,-1.918))
f0_3=sin(pos%y*vec4(.076,1.427,-2.403,-.442)+pos%z*vec4(3.012,1.584,-3.421,1.513)+&
         pos%x*vec4(.061,-3.592,-.935,.125)+vec4(7.654,-.949,-7.044,4.319))
f1_0=sin(mat(&
       [.230_wp,.244_wp,-.196_wp,-.230_wp,-.303_wp,.187_wp,.122_wp,-.410_wp,.388_wp,-.741_wp,.731_wp,.031_wp,&
       -.999_wp,-.044_wp,.107_wp,.171_wp])*f0_0+&
   mat([.141_wp,-.542_wp,-.726_wp,-.280_wp,-.445_wp,-.038_wp,-.397_wp,.348_wp,-.086_wp,.413_wp,.544_wp,.298_wp,&
        -.075_wp,.122_wp,-.315_wp,.293_wp])*f0_1+&
   mat([.329_wp,-.106_wp,-.690_wp,.012_wp,-.073_wp,-.862_wp,.381_wp,.628_wp,.378_wp,-.012_wp,-.137_wp,-.071_wp,&
        .580_wp,.060_wp,-1.094_wp,.652_wp])*f0_2+&
   mat([.616_wp,.363_wp,.211_wp,-.366_wp,-.396_wp,.388_wp,-.090_wp,.217_wp,.479_wp,-.339_wp,-.039_wp,-.132_wp,&
        .176_wp,.223_wp,.312_wp,.082_wp])*f0_3+&
   vec4(-2.491,3.220,2.193,2.989))/1.0_wp+f0_0
f1_1=sin(mat([.137_wp,.498_wp,.156_wp,.391_wp,-.059_wp,.084_wp,-.182_wp,-.027_wp,.650_wp,.302_wp,.676_wp,&
        -.202_wp,-.488_wp,.297_wp,-.136_wp,-.366_wp])*f0_0+&
   mat([.211_wp,.210_wp,-.557_wp,.021_wp,.808_wp,-.041_wp,-.726_wp,.028_wp,.003_wp,-.012_wp,-.280_wp,&
        -.539_wp,-.624_wp,.261_wp,.262_wp,-.123_wp])*f0_1+&
   mat([-.220_wp,-.208_wp,.211_wp,-.740_wp,-.009_wp,-.103_wp,.741_wp,-.820_wp,-.312_wp,.074_wp,.432_wp,&
        -.258_wp,.390_wp,-.249_wp,.375_wp,.169_wp])*f0_2+&
   mat([-.085_wp,-.105_wp,.309_wp,.291_wp,-.476_wp,-.764_wp,.048_wp,-.148_wp,-.318_wp,.182_wp,-.170_wp,&
        -.220_wp,-1.041_wp,-.116_wp,-1.077_wp,.565_wp])*f0_3+&
   vec4(1.529_wp,2.058_wp,-.449_wp,-1.700_wp))/1.0_wp+f0_1
f1_2=sin(mat([.385_wp,.698_wp,-.630_wp,-.482_wp,-.099_wp,-.552_wp,.173_wp,-.148_wp,-.711_wp,-.721_wp,&
        .342_wp,.079_wp,-.470_wp,-.177_wp,.011_wp,-.117_wp])*f0_0+&
   mat([.834_wp,-.758_wp,.162_wp,.761_wp,-.348_wp,.367_wp,.375_wp,.721_wp,.731_wp,.078_wp,.024_wp,&
        .329_wp,.340_wp,-.429_wp,-.006_wp,.084_wp])*f0_1+&
   mat([.277_wp,.015_wp,-.449_wp,-.100_wp,-.223_wp,-.897_wp,.623_wp,-.004_wp,.049_wp,-.371_wp,.040_wp,&
        -.200_wp,.367_wp,.221_wp,-.063_wp,.665_wp])*f0_2+&
   mat([1.595_wp,.428_wp,.417_wp,.268_wp,-.241_wp,.145_wp,.141_wp,-.185_wp,.430_wp,-.096_wp,.521_wp,&
        .109_wp,-.239_wp,.770_wp,-.100_wp,-.128_wp])*f0_3+&
   vec4(-1.674_wp,2.261_wp,-2.165_wp,-3.244_wp))/1.0_wp+f0_2
f1_3=sin(mat([.228_wp,-.315_wp,.459_wp,.085_wp,-.060_wp,-.070_wp,-.217_wp,-.605_wp,-.217_wp,-.237_wp,&
        .400_wp,.358_wp,.198_wp,-.722_wp,-.062_wp,-.805_wp])*f0_0+&
   mat([.108_wp,.313_wp,.793_wp,-.378_wp,-.168_wp,-.350_wp,.504_wp,.664_wp,.251_wp,-.058_wp,.087_wp,&
        -.405_wp,.068_wp,.077_wp,.148_wp,-.741_wp])*f0_1+&
   mat([.324_wp,-.130_wp,.018_wp,.100_wp,-.658_wp,.087_wp,-.077_wp,-.082_wp,.101_wp,-.103_wp,-.132_wp,&
        .159_wp,-.117_wp,-.563_wp,-.606_wp,-.278_wp])*f0_2+&
   mat([.181_wp,.584_wp,.390_wp,-.498_wp,.096_wp,-.461_wp,-.324_wp,.573_wp,-.261_wp,.849_wp,-.012_wp,&
        .017_wp,-.289_wp,.578_wp,.161_wp,-1.158_wp])*f0_3+&
   vec4(-.938_wp,-.744_wp,2.675_wp,1.662_wp))/1.0_wp+f0_3
f2_0=sin(mat([-.165_wp,-.648_wp,-.274_wp,-.443_wp,.307_wp,-.478_wp,-.204_wp,-.103_wp,-.017_wp,1.027_wp,&
        .680_wp,.228_wp,-.121_wp,.419_wp,.137_wp,-.253_wp])*f1_0+&
   mat([-.252_wp,.249_wp,.450_wp,-.249_wp,-1.085_wp,-.292_wp,.609_wp,.047_wp,-.999_wp,-.082_wp,-.005_wp,&
        -.836_wp,-.356_wp,.411_wp,-.402_wp,.699_wp])*f1_1+&
   mat([-.103_wp,-.218_wp,-.000_wp,-.639_wp,-.057_wp,.447_wp,-.502_wp,-.505_wp,-.298_wp,.029_wp,-.363_wp,&
        -.607_wp,.116_wp,-.942_wp,2.064_wp,.363_wp])*f1_2+&
   mat([.394_wp,.030_wp,.070_wp,1.402_wp,-.552_wp,-.391_wp,.255_wp,-.733_wp,-.102_wp,.263_wp,-.291_wp,&
        .190_wp,.489_wp,.614_wp,.683_wp,.061_wp])*f1_3+&
   vec4(-.896_wp,.069_wp,3.048_wp,2.889_wp))/1.4_wp+f1_0
f2_1=sin(mat([.553_wp,-.627_wp,-.284_wp,-1.032_wp,-.524_wp,.294_wp,-.078_wp,-.128_wp,.100_wp,.112_wp,.089_wp,&
        .318_wp,.560_wp,.407_wp,-.317_wp,-1.634_wp])*f1_0+&
   mat([.024_wp,.420_wp,.599_wp,.063_wp,-.317_wp,.728_wp,.197_wp,-1.037_wp,.425_wp,-.035_wp,-.240_wp,-1.530_wp,&
        -.246_wp,.333_wp,.209_wp,-1.455_wp])*f1_1+&
   mat([.514_wp,-.095_wp,.139_wp,-.551_wp,-.225_wp,-.067_wp,.568_wp,-.958_wp,.038_wp,.161_wp,-.057_wp,-.871_wp,&
        .468_wp,-.633_wp,.440_wp,1.798_wp])*f1_2+&
   mat([.185_wp,-.798_wp,.206_wp,.941_wp,1.032_wp,-.016_wp,.629_wp,-.416_wp,.026_wp,-.777_wp,.206_wp,-.360_wp,&
        -.456_wp,.463_wp,1.117_wp,.883_wp])*f1_3+&
   vec4(2.594_wp,2.602_wp,-4.067_wp,.754_wp))/1.4_wp+f1_1
f2_2=sin(mat([.601_wp,.096_wp,-.071_wp,.363_wp,-1.159_wp,-.053_wp,.049_wp,-.140_wp,.056_wp,.701_wp,-.047_wp,&
        -.031_wp,1.257_wp,-.421_wp,-.360_wp,-1.596_wp])*f1_0+&
   mat([.357_wp,.525_wp,-.130_wp,-.001_wp,-.643_wp,-.206_wp,.237_wp,-.837_wp,.813_wp,-.308_wp,.407_wp,-.157_wp,&
        1.845_wp,.737_wp,-.374_wp,.390_wp])*f1_1+&
   mat([-.569_wp,.109_wp,.476_wp,-.377_wp,.124_wp,.171_wp,-.100_wp,-.606_wp,-.560_wp,.414_wp,.176_wp,.483_wp,&
        -2.082_wp,-.375_wp,.928_wp,.125_wp])*f1_2+&
   mat([-.599_wp,.549_wp,-.414_wp,.642_wp,-.525_wp,.025_wp,.385_wp,-.797_wp,1.270_wp,.342_wp,-.068_wp,.965_wp,&
        -1.026_wp,-.177_wp,.890_wp,.156_wp])*f1_3+&
   vec4(-4.226_wp,1.935_wp,3.830_wp,-2.055_wp))/1.4_wp+f1_2
f2_3=sin(mat([.463_wp,1.053_wp,.008_wp,-.354_wp,-.124_wp,-.566_wp,-.734_wp,.673_wp,.970_wp,-.193_wp,.377_wp,&
        .283_wp,-.481_wp,-1.345_wp,.024_wp,.687_wp])*f1_0+&
   mat([-.358_wp,-.069_wp,.358_wp,.373_wp,.745_wp,-.056_wp,.054_wp,1.070_wp,-.459_wp,-.404_wp,-.654_wp,.665_wp,&
        -.400_wp,.127_wp,.422_wp,.202_wp])*f1_1+&
   mat([.255_wp,-.513_wp,-.741_wp,-.345_wp,-.404_wp,.188_wp,.118_wp,.486_wp,.783_wp,.610_wp,-.790_wp,-.099_wp,&
        -.555_wp,-3.111_wp,-1.241_wp,-.215_wp])*f1_2+&
   mat([.776_wp,.131_wp,.035_wp,.434_wp,.355_wp,-2.543_wp,-.158_wp,-.140_wp,.454_wp,1.154_wp,.156_wp,.543_wp,&
        .555_wp,-.726_wp,.369_wp,-1.032_wp])*f1_3+&
   vec4(1.108_wp,.703_wp,-.637_wp,.098_wp))/1.4_wp+f1_3
out= (f2_0.dot.vec4(-.041,-.048,-.054,.034))+&
    (f2_1.dot.vec4(.066,.103,.066,-.022))+&
    (f2_2.dot.vec4(.011,.059,-.092,.064))+&
    (f2_3.dot.vec4(-.059,-.031,.063,-.052))+&
    (0.079)
end function eval_neural


        type(vector) function calcNormal(p, obj)

            type(vector), intent(IN) :: p
            class(sdf) :: obj

            real(kind=wp) :: h
            type(vector) :: xyy, yyx, yxy, xxx

            h = 1e-5_wp
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

        function model_init(array, func, kopt) result(out)

            type(model) :: out

            procedure(op) :: func
            type(container),          intent(IN) :: array(:)
            real(kind=wp), optional,  intent(IN) :: kopt

            out%array = array
            out%func => func
            if(present(kopt))then
                out%k = kopt
            else
                out%k = 0._wp
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

        function eval_model(this, pos) result(res)

            class(model) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            integer :: i

            res = this%array(1)%p%evaluate(pos)
            do i = 2, size(this%array)
                res = this%func(res, this%array(i)%p%evaluate(pos), this%k)
            end do

        end function eval_model


        function cylinder_init(a, b, radius, mus, mua, hgg, n, layer, transform) result(out)
                
            type(cylinder) :: out
            
            real(kind=wp),           intent(IN) :: radius, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function cylinder_init

        function eval_cylinder(this, pos) result(res)

            class(cylinder) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = cylinder_fn(p, this%a, this%b, this%radius)

        end function eval_cylinder

        function cylinder_fn(p, a, b, r) result(res)
            !p = pos
            !a = pt1
            !b = pt2
            !r = radius
            !draws cylinder along the axis between 2 points a and b

            type(vector),  intent(IN) :: p, a, b
            real(kind=wp), intent(IN) :: r
            real(kind=wp) :: res

            type(vector)  :: ba, pa
            real(kind=wp) :: x, y, x2, y2, d, baba, paba

            ba = b - a
            pa = p - a
            baba = ba .dot. ba
            paba = pa .dot. ba
            x = length(pa * baba - ba*paba) - r*baba
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

        end function cylinder_fn

        function box_init(lengths, mus, mua, hgg, n, layer, transform) result(out)
                
            type(box) :: out
            
            type(vector),            intent(IN) :: lengths
            real(kind=wp),           intent(IN) :: mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function box_init

        function box_init_vec(widths, mus, mua, hgg, n, layer, transform) result(out)
                
            type(box) :: out
            
            type(vector),            intent(IN) :: widths
            real(kind=wp),           intent(IN) :: mus, mua, hgg, n
            integer,                 intent(IN) :: layer
            real(kind=wp), optional, intent(in) :: transform(4, 4)

            out = box_init(widths, mus, mua, hgg, n, layer, transform)

        end function box_init_vec

        function box_init_scal(width, mus, mua, hgg, n, layer, transform) result(res)
                
            type(box) :: res
            
            real(kind=wp),           intent(IN) :: width, mus, mua, hgg, n
            integer,                 intent(IN) :: layer
            real(kind=wp), optional, intent(in) :: transform(4, 4)

            type(vector) :: widths

            widths = vector(width, width, width)

            res = box_init(widths, mus, mua, hgg, n, layer, transform)

        end function box_init_scal

        function eval_box(this, pos) result(res)

            class(box) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = box_fn(p, this%lengths)

        end function eval_box

        function box_fn(p, b) result(res)

            type(vector), intent(IN) :: p, b
            real(kind=wp) :: res

            type(vector) :: q

            q = abs(p) - b
            res = length(max(q, 0._wp)) + min(max(q%x, max(q%y, q%z)), 0._wp)

        end function box_fn


        function bevel_box_init(size, box_r, mus, mua, hgg, n, layer, transform) result(out)
                
            type(bevel_box) :: out
            
            type(vector),            intent(IN) :: size
            real(kind=wp),           intent(IN) :: box_r, mus, mua, hgg, n
            integer,                 intent(IN) :: layer
            real(kind=wp), optional, intent(in) :: transform(4, 4)

            real(kind=wp) :: t(4, 4)

            if(present(transform))then
                t = transform
            else
                t = identity()
            end if

            out%size = size
            out%box_r = box_r
            out%layer = layer
            out%transform = t

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function bevel_box_init

        function eval_bevel_box(this, pos) result(res)

            class(bevel_box) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = bevel_box_fn(p, this%size, this%box_r)

        end function eval_bevel_box

        function bevel_box_fn(p, size, box_r) result(res)

            type(vector),  intent(IN) :: p, size
            real(kind=wp), intent(IN) :: box_r
            real(kind=wp) :: res

            type(vector)  :: box_edge, dd
            real(kind=wp) :: maxdd, ddd

            box_edge = size - box_r * 0.5_wp
            dd = abs(p) - box_edge

            maxdd = max(max(dd%x, dd%y), dd%z)
            maxdd = min(maxdd, 0.0_wp)

            dd = max(dd, 0.0_wp)
            ddd = (length(dd) - box_r)
            ddd  = ddd + maxdd

            res = ddd

        end function bevel_box_fn



        function sphere_init(radius, mus, mua, hgg, n, layer, transform) result(out)
                
            type(sphere) :: out
            
            real(kind=wp),           intent(IN) :: radius, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function sphere_init

        function eval_sphere(this, pos) result(res)

            class(sphere) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = sphere_fn(p, this%radius)

        end function eval_sphere

        function sphere_fn(p, c) result(res)

            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: c
            real(kind=wp) :: res

            res = sqrt(p%x**2+p%y**2+p%z**2) - c

        end function sphere_fn


        function torus_init(oradius, iradius, mus, mua, hgg, n, layer, transform) result(out)
        
            type(torus) :: out
            
            real(kind=wp),           intent(IN) :: oradius, iradius, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function torus_init

        function eval_torus(this, pos) result(res)

            class(torus) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = torus_fn(p, this%oradius, this%iradius)

        end function eval_torus

        function torus_fn(p, or, ir) result(res)

            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: or, ir
            real(kind=wp) :: res

            type(vector) :: q

            q = vector(length(vector(p%x, 0._wp, p%z)) - or, p%y, 0._wp)
            res = length(q) - ir

        end function torus_fn


        function triprism_init(h1, h2, mus, mua, hgg, n, layer, transform) result(out)
        !h1 is height
        !h2 is length
        !        
            type(triprism) :: out
            
            real(kind=wp),           intent(IN) :: h1, h2, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function triprism_init

        function eval_triprism(this, pos) result(res)

            class(triprism) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = triprism_fn(p, this%h1, this%h2)

        end function eval_triprism

        function triprism_fn(p, h1, h2) result(res)

            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: h1, h2
            real(kind=wp) :: res

            type(vector) :: q

            q = abs(p)
            res = max(q%z - h2, max(q%x*.866025_wp + p%y*.5_wp, -p%y) - h1*.5_wp) 

        end function triprism_fn

        function cone_init(a, b, ra, rb, mus, mua, hgg, n, layer, transform) result(out)
        
            type(cone) :: out
            
            type(vector),            intent(IN) :: a, b
            real(kind=wp),           intent(IN) :: ra, rb, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function cone_init

        function eval_cone(this, pos) result(res)

            class(cone) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = cone_fn(p, this%a, this%b, this%ra, this%rb)

        end function eval_cone

        function cone_fn(p, a, b, ra, rb) result(res)

            use utils, only : clamp

            type(vector),  intent(IN) :: p, a, b
            real(kind=wp), intent(IN) :: ra, rb
            real(kind=wp) :: res

            real(kind=wp) :: rba, baba, papa, paba, x, cax, cay, k, f, cbx, cby, s

            rba = rb - ra
            baba = (b-a) .dot. (b-a)
            papa = (p-a) .dot. (p-a)
            paba =  ((p-a) .dot. (b-a))/ baba
            x = sqrt(papa - baba*paba**2)
            if(paba < 0.5_wp)then
                cax = max(0._wp, x - ra)
            else
                cax = max(0._wp, x - rb)
            end if
            cay = abs(paba - 0.5_wp) - .5_wp
            k = rba**2 + baba
            f = clamp((rba * (x - ra) + paba*baba) / k, 0._wp, 1._wp)
            cbx = x - ra - f*rba
            cby = paba - f
            if(cbx < 0._wp .and. cay < 0._wp)then
                s = -1._wp
            else
                s = 1._wp
            end if 
            res = s * sqrt(min(cax**2 + baba*cay**2, cbx**2 + baba*cby**2)) 

        end function cone_fn

        function capsule_init(a, b, r, mus, mua, hgg, n, layer, transform) result(out)
        
            type(capsule) :: out
            
            type(vector),            intent(IN) :: a, b
            real(kind=wp),           intent(IN) :: r, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function capsule_init

        function eval_capsule(this, pos) result(res)

            class(capsule) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = capsule_fn(p, this%a, this%b, this%r)

        end function eval_capsule

        function capsule_fn(p, a, b, r) result(res)

            use utils, only : clamp

            type(vector),  intent(IN) :: p, a, b
            real(kind=wp), intent(IN) :: r
            real(kind=wp) :: res

            type(vector) :: pa, ba
            real(kind=wp) :: h

            pa = p - a
            ba = b - a
            h = clamp((pa .dot. ba) / (ba .dot. ba), 0._wp, 1._wp)
            res = length(pa - ba*h) - r

        end function capsule_fn


        function plane_init(a, mus, mua, hgg, n, layer, transform) result(out)
        
            type(plane) :: out
            
            type(vector),             intent(IN) :: a
            real(kind=wp),            intent(IN) :: mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function plane_init

        function eval_plane(this, pos) result(res)

            class(plane) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = plane_fn(p, this%a)

        end function eval_plane

        function plane_fn(p, n) result(res)

            type(vector), intent(IN) :: p, n
            real(kind=wp) :: res

            !n must be normalised
            res = (p .dot. n)

        end function plane_fn


        function moon_init(d, ra, rb, mus, mua, hgg, n, layer, transform) result(out)
        
            type(moon) :: out
            
            real(kind=wp),            intent(IN) :: d, ra, rb, mus, mua, hgg, n
            integer,                  intent(IN) :: layer
            real(kind=wp),  optional, intent(IN) :: transform(4, 4)

            real(kind=wp) :: t(4, 4)

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
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function moon_init

        function eval_moon(this, pos) result(res)

            class(moon) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = moon_fn(p, this%d, this%ra, this%rb)

        end function eval_moon

        function moon_fn(p, d, ra, rb) result(res)

            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: d, ra, rb
            real(kind=wp) :: res

            type(vector)  :: pos
            real(kind=wp) :: a, b, ra2, rb2, d2

            pos = vector(p%x, abs(p%y), 0._wp)
            ra2 = ra*ra
            rb2 = rb*rb
            d2 = d*d
            a = (ra2 - rb2 + d2) / (2._wp*d)
            b = sqrt(max(ra2 - a**2, 0._wp))
            if(d*(pos%x*b - pos%y*a) > d2*max(b - pos%y, 0._wp))then
                res = length(pos - vector(a, b, 0._wp))
            else
                res = max(-length(pos) - ra, length(pos - vector(d, 0._wp, 0._wp)) - rb)
            end if

        end function moon_fn

        function egg_init(r1, r2, h, mus, mua, hgg, n, layer, transform) result(out)
        ! makes a Moss egg. https://www.shadertoy.com/view/WsjfRt
        ! R1 controls "fatness" of the egg. Actually controls the base circle radius.
        ! R2 contorls the pointiness of the egg. Actually controls radius of top circle.
        ! h controls the height of the egg. Actually controls y position of top circle.
            type(egg) :: out
            
            real(kind=wp),            intent(IN) :: r1, r2, h, mus, mua, hgg, n
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

            out%mus = mus
            out%mua = mua
            out%kappa = mus + mua
            if(out%mua < 1e-9_wp)then
                out%albedo = 1._wp
            else
                out%albedo = mus / out%kappa
            end if
            out%hgg = hgg
            out%g2 = hgg**2
            out%n = n

        end function egg_init


        function eval_egg(this, pos) result(res)

            class(egg) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            type(vector) :: p

            p = pos .dot. this%transform
            res = egg_fn(p, this%r1, this%r2, this%h)

        end function eval_egg

        function egg_fn(p, r1, r2, h) result(res)
        !https://www.shadertoy.com/view/WsjfRt

            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: r1, r2, h
            real(kind=wp) :: res

            real(kind=wp) :: r, l, h_in
            type(vector) :: p_in

            p_in = p

            p_in%x = abs(p%x)
            r = r1 - r2
            h_in = h + r
            l = (h_in**2 - r**2) / (2._wp * r)

            if(p_in%y <= 0._wp)then
                res = length(p_in) - r1
            else
                if((p_in%y - h_in) * l > p_in%x*h_in)then
                    res = length(p_in - vector(0._wp, h_in, 0._wp)) - ((r1+l) - length(vector(h_in,l, 0._wp)))
                else
                    res = length(p_in + vector(l, 0._wp, 0._wp)) - (r1+l)
                end if
            end if
        end function egg_fn


        function translate(o) result(out)

            type(vector), intent(IN) :: o

            real(kind=wp) :: out(4, 4)

            out(:, 1) = [1._wp, 0._wp, 0._wp, o%x] 
            out(:, 2) = [0._wp, 1._wp, 0._wp, o%y] 
            out(:, 3) = [0._wp, 0._wp, 1._wp, o%z] 
            out(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp] 

        end function translate

        function union(d1, d2, k) result(res)

            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res

            res = min(d1, d2)
        end function union


        function SmoothUnion(d1, d2, k) result(res)

            ! use utils, only : mix, clamp

            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res

            real(kind=wp) :: h

            h = max(k - abs(d1 - d2), 0._wp) / k
            res = min(d1, d2) - h*h*h*k*(1._wp/6._wp)
            ! h = clamp(0.5 +.5*(d2-d1)/k, 0., 1.)
            ! SmoothUnion = mix(d2, d1, h) - k*h*(1.-h)

        end function SmoothUnion

        function subtraction(d1, d2, k) result(res)

            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res

            res = max(-d1, d2)

        end function subtraction

        function intersection(d1, d2, k) result(res)

            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res

            res = max(d1, d2)

        end function intersection

        type(elongate) function elongate_init(prim, size) result(out)

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

        function eval_elongate(this, pos) result(res)

            class(elongate) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = elongate_fn(pos, this%size, this%prim)

        end function eval_elongate

        function elongate_fn(p, size, prim) result(res)

            class(sdf) :: prim

            type(vector), intent(IN) :: size
            type(vector), intent(IN) :: p
            real(kind=wp) :: res

            real(kind=wp) :: w
            type(vector) :: q

            q = abs(p) - size
            w = min(max(q%x, max(q%y, q%z)), 0._wp)

            res = prim%evaluate(max(q, 0._wp)) + w

        end function elongate_fn

        type(bend) function bend_init(prim, k) result(out)

            real(kind=wp), intent(IN) :: k
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

        function eval_bend(this, pos) result(res)

            class(bend) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = bend_fn(pos, this%k, this%prim)

        end function eval_bend

        function bend_fn(p, k, prim) result(res)

            class(sdf) :: prim

            real(kind=wp), intent(IN) :: k
            type(vector),  intent(IN) :: p
            real(kind=wp) :: res

            real(kind=wp) :: c, s, x2, y2, z2

            c = cos(k * p%x)
            s = sin(k * p%x)
            x2 = c * p%x - s * p%y
            y2 = s * p%x + c * p%y
            z2 = p%z

            res = prim%evaluate(vector(x2, y2, z2))

        end function bend_fn

        type(displacement) function displacement_init(prim, func) result(out)

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

        function eval_disp(this, pos) result(res)

            class(displacement) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = displacement_fn(pos, this%prim, this%func)

        end function eval_disp

        function displacement_fn(p, prim, disp) result(res)

            class(sdf) :: prim
            procedure(primitive) :: disp
            type(vector), intent(IN) :: p
            real(kind=wp) :: res

            real(kind=wp) :: d1, d2

            d1 = prim%evaluate(p)
            d2 = disp(p)

            res = d1 + d2

        end function displacement_fn

        type(twist) function twist_init(prim, k) result(out)

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

        function eval_twist(this, pos) result(res)

            class(twist) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = twist_fn(pos, this%k, this%prim)

        end function eval_twist

        function twist_fn(p, k, prim) result(res)

            class(sdf) :: prim
            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: k
            real(kind=wp) :: res

            real(kind=wp) :: c, s, x2, y2, z2

            c = cos(k * p%z)
            s = sin(k * p%z)
            x2 = c*p%x - s*p%y
            y2 = s*p%x + c*p%y
            z2 = p%z

            res = prim%evaluate(vector(x2, y2, z2))

        end function twist_fn


        type(repeat) function repeat_init(prim, c, la, lb) result(out)

            class(sdf), target :: prim
            type(vector),  intent(IN) :: la, lb
            real(kind=wp), intent(IN) :: c

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

        function eval_repeat(this, pos) result(res)

            class(repeat) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = repeat_fn(pos, this%c, this%la, this%lb, this%prim)

        end function eval_repeat

        function repeat_fn(p, c, la, lb, prim) result(res)

            ! use utils, only : clamp

            class(sdf) :: prim
            type(vector),  intent(IN) :: p, la, lb
            real(kind=wp), intent(IN) :: c
            real(kind=wp) :: res

            type(vector) :: q

            error stop "Not implmented as no vector dependacny in utils yet!"
            ! q = p - c*clamp(nint(p/c), la, lb)
            res = prim%evaluate(q)

        end function repeat_fn

        type(extrude) function extrude_init(prim, h) result(out)

            class(sdf), target :: prim
            real(kind=wp), intent(IN)   :: h

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

        function eval_extrude(this, pos) result(res)

            class(extrude) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = extrude_fn(pos, this%h, this%prim)

        end function eval_extrude


        function extrude_fn(p, h, prim) result(res)

            class(sdf) :: prim
            type(vector),  intent(IN) :: p
            real(kind=wp), intent(IN) :: h
            real(kind=wp) :: res

            type(vector)  :: w
            real(kind=wp) :: d

            d = prim%evaluate(p)
            w = vector(d, abs(p%z) - h, 0._wp)
            res = min(max(w%x, w%y), 0._wp) + length(max(w, 0._wp))

        end function extrude_fn

        type(revolution) function revolution_init(prim, o) result(out)

            class(sdf), target :: prim
            real(kind=wp), intent(IN) :: o

            out%o = o
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

            end function revolution_init

        function eval_revolution(this, pos) result(res)

            class(revolution) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = revolution_fn(pos, this%o, this%prim)

        end function eval_revolution

        real(kind=wp) function revolution_fn(p, o, prim) result(res)

            class(sdf) :: prim
            type(vector), intent(IN) :: p
            real(kind=wp), intent(IN) :: o

            type(vector) :: pxz, q

            pxz = vector(p%x, p%z, 0._wp)

            q = vector(length(pxz) - o, p%y, 0._wp)
            res = prim%evaluate(q)

        end function revolution_fn


        type(onion) function onion_init(prim, thickness) result(out)

            class(sdf), target :: prim
            real(kind=wp), intent(IN) :: thickness

            out%thickness = thickness
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

        end function onion_init

        function eval_onion(this, pos) result(res)

            class(onion) :: this
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res

            res = onion_fn(pos, this%thickness, this%prim)

        end function eval_onion

        real(kind=wp) function onion_fn(p, thickness, prim) result(res)

            class(sdf) :: prim
            type(vector), intent(IN) :: p
            real(kind=wp), intent(IN) :: thickness

            res = abs(prim%evaluate(p)) - thickness

        end function onion_fn

        function rotate_x(angle) result(r)
        ! rotation funcs from https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/
            use utils, only : deg2rad
            
            real(kind=wp), intent(IN) :: angle
            
            real(kind=wp) :: r(4, 4), c, s, a

            a = deg2rad(angle)
            c = cos(a)
            s = sin(a)

            r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
            r(:, 2) = [0._wp, c,  s,  0._wp]
            r(:, 3) = [0._wp,-s,  c,  0._wp]
            r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

        end function rotate_x

        function rotate_y(angle) result(r)
            
            use utils, only : deg2rad
            
            real(kind=wp), intent(IN) :: angle

            real(kind=wp) :: r(4, 4), c, s, a

            a = deg2rad(angle)
            c = cos(a)
            s = sin(a)

            r(:, 1) = [c,  0._wp, -s,  0._wp]
            r(:, 2) = [0._wp, 1._wp,  0._wp, 0._wp]
            r(:, 3) = [s,  0._wp,  c,  0._wp]
            r(:, 4) = [0._wp, 0._wp,  0._wp, 1._wp]

        end function rotate_y

        function rotate_z(angle) result(r)
            
            use utils, only : deg2rad
            
            real(kind=wp), intent(IN) :: angle

            real(kind=wp) :: r(4, 4), c, s, a

            a = deg2rad(angle)
            c = cos(a)
            s = sin(a)

            r(:, 1) = [c, -s,  0._wp, 0._wp]
            r(:, 2) = [s,  c,  0._wp, 0._wp]
            r(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
            r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

        end function rotate_z

        function rotmat(axis, angle)
        ! http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/

            type(vector),  intent(in) :: axis
            real(kind=wp), intent(in) :: angle

            type(vector) :: axist

            real(kind=wp) :: rotmat(4, 4), s, c, oc

            axist = axis%magnitude()
            s = sqrt(1. - angle**2)
            c = angle
            oc = 1._wp - c

        rotmat(:, 1) = [oc * axist%x * axist%x + c, oc * axist%x * axist%y - axist%z * s,&
                        oc * axist%z * axist%x + axist%y * s, 0.0_wp]
        rotmat(:, 2) = [oc * axist%x * axist%y + axist%z * s, oc * axist%y * axist%y + c,&
                        oc * axist%y * axist%z - axist%x * s, 0.0_wp]
        rotmat(:, 3) = [oc * axist%z * axist%x - axist%y * s, oc * axist%y * axist%z + axist%x * s,&
                        oc * axist%z * axist%z + c, 0.0_wp]
        rotmat(:, 4) = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]

        end function rotmat

        ! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        ! https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        function rotationAlign(a, b) result(res)
            
            type(vector), intent(in) :: a, b
            
            type(vector)  :: v
            real(kind=wp) :: c, k, res(4, 4), v_x(4, 4), v_x2(4, 4)
            
            
            v = a .cross. b
            c = a .dot. b
            k = 1._wp / (1._wp + c)
            
            !skew-symmetric matrix
            v_x(:, 1) = [0._wp     , -1._wp*v%z, v%y       , 0._wp]
            v_x(:, 2) = [v%z       , 0._wp     , -1._wp*v%x, 0._wp]
            v_x(:, 3) = [-1._wp*v%y, v%x       , 0._wp     , 0._wp]
            v_x(:, 4) = [0._wp, 0._wp, 0._wp, 0._wp]
            
            v_x2 = matmul(v_x, v_x)
            res = identity() + v_x + v_x2*k

        end function rotationAlign

        function identity() result(r)
                        
            real(kind=wp) :: r(4, 4)

            r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
            r(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
            r(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
            r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

        end function identity


        function skewSymm(a) result(out)

            type(vector), intent(in) :: a
            real(kind=wp) :: out(4,4)

            out(:, 1) = [0._wp, -a%z, a%y, 0._wp]
            out(:, 2) = [a%z, 0._wp, -a%x, 0._wp]
            out(:, 3) = [-a%y, a%x, 0._wp, 0._wp]
            out(:, 4) = [0._wp, 0._wp, 0._wp, 0._wp]

        end function skewSymm

        subroutine render_vec(cnt, extent, samples, fname)
                                    
            type(container),        intent(IN) :: cnt(:)
            integer,                intent(IN) :: samples(3)
            type(vector),           intent(IN) :: extent
            character(*), optional, intent(IN) :: fname

            character(len=:), allocatable  :: filename

            if(present(fname))then
                filename = fname
            else
                filename = "model.dat"
            end if

            call render_sub(cnt, extent, samples, filename)

        end subroutine render_vec

        subroutine render_scaler(cnt, extent, sample, fname)

            type(container),        intent(IN) :: cnt(:)
            integer,                intent(IN) :: sample
            type(vector),           intent(IN) :: extent
            character(*), optional, intent(IN) :: fname

            character(len=:), allocatable  :: filename
            integer :: samples(3)

            if(present(fname))then
                filename = fname
            else
                filename = "model.dat"
            end if

            samples = [sample, sample, sample]

            call render_sub(cnt, extent, samples, filename)

        end subroutine render_scaler

        subroutine render_sub(cnt, extent, samples, fname)

            use utils,     only : pbar
            use constants, only : fileplace
            use writer_mod
                        
            type(container),        intent(IN) :: cnt(:)
            integer,                intent(IN) :: samples(3)
            type(vector),           intent(IN) :: extent
            character(*),           intent(IN) :: fname

            type(vector)               :: pos, wid
            integer                    :: i, j, k, u, id
            real(kind=wp)              :: x, y, z, ds(size(cnt)), ns(3), minvalue
            real(kind=wp), allocatable :: image(:, :, :)
            type(pbar)                 :: bar

            ns = nint(samples / 2._wp)
            allocate(image(samples(1), samples(2), samples(3)))
            wid = vector(extent%x/ns(1), extent%y/ns(2), extent%z/ns(3))
            bar = pbar(samples(1))
!$omp parallel default(none) shared(cnt, ns, wid, image, samples, bar)&
!$omp private(i, x, y, z, pos, j, k, u, ds, id, minvalue)
!$omp do
            do i = 1, samples(1)
                x = (i-ns(1)) *wid%x
                do j = 1, samples(2)
                    y = (j-ns(2)) *wid%y
                    do k = 1, samples(3)
                        z = (k-ns(3)) * wid%z
                        pos = vector(x, y, z)
                        ds = 0._wp
                        do u = 1, size(ds)
                            ds(u) = cnt(u)%p%evaluate(pos)
                        end do
                        if(all(ds > 0._wp))then
                            id=0
                        else
                            if(all(ds < 0._wp))then
                                id = cnt(maxloc(ds,dim=1))%p%layer
                            else
                                ds = -1.*ds
                                minvalue = 1000000._wp
                                do u = 1, size(ds)
                                    if(ds(u) > 0. .and. ds(u) < minvalue)then
                                        minvalue = ds(u)
                                        id = u
                                    end if
                                end do
                                ! id = cnt(minloc(ds),dim=1))%p%layer
                            end if
                        end if
                        ! if(id == 0._wp)then
                        !     image(i,j,k)=-99._wp
                        ! else
                            image(i, j, k) = real(id)!cnt(id)%p%mua
                        ! end if
                    end do
                end do
                call bar%progress()
            end do
!$OMP end  do
!$OMP end parallel
            call write_fluence(image, trim(fileplace)//fname, overwrite=.true.)
        end subroutine render_sub
end module sdfs