module subs

    use constants, only : wp

    implicit none

    public  :: setup_simulation
    private :: directory, alloc_array, zarray

    contains

        subroutine setup_simulation(sdfarray, dict)
        ! Read in parameters
        ! Setup up various simulation parameters and routines
        !
            use vector_class
            use sdfs,          only : container
            use fhash,         only : fhash_tbl_t
            use sim_state_mod, only : settings => state


            implicit none

            type(fhash_tbl_t),  optional, intent(IN)  :: dict
            type(container), allocatable, intent(OUT) :: sdfarray(:)

            !set directory paths
            call directory

            !allocate and set arrays to 0
            call alloc_array(settings%grid%nxg, settings%grid%nyg, settings%grid%nzg)
            call zarray()

            ! setup geometry using SDFs
            select case(settings%experiment)
                case("logo")
                    sdfarray = setup_logo()
                case("neural")
                    sdfarray = setup_neural_sdf()
                case("omg")
                    sdfarray = setup_omg_sdf()
                case("scat_test")
                    sdfarray = setup_scat_test(dict)
                case("skin")
                    sdfarray = setup_skin_model()
                case("aptran")
                    sdfarray = setup_sphere()
                case("exp")
                    sdfarray = setup_exp(dict)
                case("jacques")
                    sdfarray = setup_jacques()
                case("vessels")
                    sdfarray = get_vessels()
                case("lens")
                    sdfarray = lens_test_setup()
                case("blobby")
                    sdfarray = blobby()
                case("sphere_scene")
                    sdfarray = setup_sphere_scene(dict)
                case default
                    error stop "no such routine"
            end select

        end subroutine setup_simulation


        function setup_sphere_scene(dict) result(array)

            use sdfs,         only : container, sphere, box, translate
            use fhash,        only : fhash_tbl_t, key=>fhash_key
            use random,       only : ranu
            use vector_class, only : vector, invert

            implicit none

            type(fhash_tbl_t), intent(in) :: dict
            type(container), allocatable :: array(:)
            
            type(sphere), target, save, allocatable :: sphs(:)
            type(box),    target, save :: bbox
            integer :: num_spheres, i
            real(kind=wp) :: t(4,4), mus, mua, hgg, n, radius
            type(vector) :: pos

            call dict%get(key("num_spheres"), num_spheres)
            allocate(sphs(num_spheres+1))

            mus = 1e-17_wp
            mua = 1e-17_wp
            hgg = 0.0_wp
            n   = 1.0_wp
            bbox = box(2._wp, mus, mua, hgg, n, num_spheres+1)

            do i = 1, num_spheres
                radius = ranu(0.001_wp, 0.25_wp)
                pos = vector(ranu(-1._wp+radius, 1._wp-radius), ranu(-1._wp+radius, 1._wp-radius), ranu(-1._wp+radius, 1._wp-radius))
                t = invert(translate(pos))
                mus = ranu(1._wp, 50._wp)
                mua = ranu(0.01_wp, 1._wp)
                hgg = 0.9_wp
                n = 1.37_wp
                sphs(i) = sphere(radius, mus, mua, hgg, n, i, transform=t)
            end do

            allocate(array(num_spheres+1))
            do i = 1, num_spheres
                allocate(array(i)%p, source=sphs(i))
                array(i)%p => sphs(i)
            end do
            
            allocate(array(num_spheres+1)%p, source=bbox)
            array(num_spheres+1)%p => bbox

        end function setup_sphere_scene

        function blobby() result(array)
            
            use sdfs, only : box, capsule, sphere, translate, rotate_x, rotate_y, rotate_z, container, model, smoothunion, model_init
            use vector_class, only : vector, invert

            implicit none
            
            type(container), allocatable :: array(:)
            type(container), target, save, allocatable :: cnta(:)
            type(sphere), target, save :: sph(7)
            type(capsule), target, save :: cap(3)
            type(box), target, save :: boxy
            type(model), target, save :: m
            type(vector) :: a, b
            real(kind=wp) :: t(4,4)
            real(kind=wp) :: mus, mua, hgg,n 
            integer :: layer, i
            
            mua = 0.1_wp
            mus = 10._wp
            n = 1.3_wp
            layer = 1

            sph(1) = sphere(1.5_wp, mus, mua, hgg, n, layer)

            a = vector(0._wp, 0._wp, 3._wp)
            t = invert(translate(a))
            sph(2) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)
            b = vector(0._wp, 0._wp, -3._wp)
            t = invert(translate(b))
            sph(3) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)

            a = vector(0._wp, 3._wp, 0._wp)
            t = invert(translate(a))
            sph(4) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)
            b = vector(0._wp, -3._wp, 0._wp)
            t = invert(translate(b))
            sph(5) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)

            a = vector(3._wp, 0._wp, 0._wp)
            t = invert(translate(a))
            sph(6) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)
            b = vector(-3._wp, 0._wp, 0._wp)
            t = invert(translate(b))
            sph(7) = sphere(0.75_wp, mus, mua, hgg, n, layer, transform=t)


            a = vector(0._wp, 0._wp, -3._wp)
            b = vector(0._wp, 0._wp, 3._wp)
            cap(1) = capsule(a, b, 0.5_wp, mus, mua, hgg, n, layer)
            a = vector(-3._wp, 0._wp, 0._wp)
            b = vector(3._wp, 0._wp, 0._wp)
            cap(2) = capsule(a, b, 0.5_wp, mus, mua, hgg, n, layer)
            a = vector(0._wp, -3._wp, 0._wp)
            b = vector(0._wp, 3._wp, 0._wp)
            cap(3) = capsule(a, b, 0.5_wp, mus, mua, hgg, n, layer)

            boxy = box(10._wp, 0.0_wp, 1.e-17_wp, 0.0_wp, n=1.0_wp, layer=2)

            allocate(array(2), cnta(10))
            do i = 1, 7
                allocate(cnta(i)%p, source=sph(i))
                cnta(i)%p => sph(i)
            end do
            
            do i = 1, 3
                allocate(cnta(i+7)%p, source=cap(i))
                cnta(i+7)%p => cap(i)
            end do
            m = model_init(cnta, SmoothUnion, 0.5_wp)

            allocate(array(1)%p, source=m)
            array(1)%p => m
            allocate(array(2)%p, source=boxy)
            array(2)%p => boxy
        end function blobby


        function lens_test_setup() result(array)
        !create lens geometry from tran and Jacques interpolated normals approach paper
            use sdfs, only : sphere, container, box, model, intersection, model_init, translate
            use vector_class, only : vector, invert

            implicit none


            type(container),target, save :: cnta(2)
            type(container) :: array(2)
            type(box),     target, save :: bbox
            type(sphere), target, save :: sph(2)
            type(model), target, save :: m
            
            type(vector)  :: a, b
            real(kind=wp) :: hgg, mus, mua, n1, n2, t(4,4)
            integer       :: layer, i

            mus = 0._wp
            mua = 1.e-17_wp
            hgg = 0.0_wp
            n1 = 1.33_wp
            n2 = 1.52_wp
            layer = 1

            bbox = box(vector(2._wp, 2._wp, 2._wp), mus, mua, hgg, n1, 2) 
            !sphere1
            a = vector(0._wp, 0._wp, 0.05_wp)
            t = 0._wp
            t = invert(translate(a))
            sph(1) = sphere(.2_wp, mus, mua, hgg, n2, 1, transform=t)
            !sphere2
            b = vector(0._wp, 0._wp, -0.15_wp)
            t = 0._wp
            t = invert(translate(b))
            sph(2) = sphere(.25_wp, mus, mua, hgg, n2, 1, transform=t)

            do i = 1, size(cnta)
                allocate(cnta(i)%p, source=sph(i))
                cnta(i)%p => sph(i)
            end do

            m = model_init(cnta, intersection)
            allocate(array(1)%p, source=m)
            array(1)%p => m

            allocate(array(2)%p, source=bbox)
            array(2)%p => bbox

        end function lens_test_setup


        function setup_logo() result(array)
        ! setup uni crest geometry
        !
        !
            use vector_class
            use sdfs, only : container, box, segment, extrude, model, union, model_init

            implicit none


            type(container), allocatable :: cnta(:), array(:)
            type(box),     target, save :: bbox
            type(segment), allocatable, target, save :: seg(:)
            type(extrude), allocatable, target, save :: ex(:)
            type(model), target, save :: m
            
            type(vector)  :: a, b
            real(kind=wp) :: hgg, mus, mua, n
            integer       :: layer, i
            logical       :: fexists

            allocate(array(2), cnta(725), seg(725), ex(725))

            mus = 10._wp
            mua = .1_wp
            hgg = 0.9_wp
            n = 1.5_wp
            layer = 1

            bbox = box(vector(10._wp, 10._wp, 2.001_wp), 0._wp, 0._wp, 0._wp, 1._wp, 2) 
            inquire(file="res/svg.f90", exist=fexists)
            if(.not. fexists)error stop "need to generate svg.f90 and place in res/"
            include "../res/svg.f90"

            do i = 1, size(cnta)
                allocate(cnta(i)%p, source=ex(i))
                cnta(i)%p => ex(i)
            end do

            m = model_init(cnta, union)
            allocate(array(1)%p, source=m)
            array(1)%p => m

            allocate(array(2)%p, source=bbox)
            array(2)%p => bbox

        end function setup_logo


        function setup_neural_sdf() result(array)
        ! Setup up neural SDF geometry
        !
        !
            use vector_class
            use sdfs, only : container, box, neural

            implicit none

            type(container)         :: array(2)
            type(box), target, save :: bbox
            type(neural), target, save :: neu
            
            real(kind=wp) :: hgg

            hgg = 0.9_wp

            bbox = box(vector(10._wp, 10._wp, 2.001_wp), 0._wp, 0._wp, 0._wp, 1._wp, 2) 
            !420
            neu = neural(82._wp/(1._wp-hgg), 1.8_wp, hgg, 1.38_wp, 1)

            allocate(array(1)%p, source=neu)
            allocate(array(2)%p, source=bbox)

            array(1)%p => neu
            array(2)%p => bbox

        end function setup_neural_sdf


        function setup_jacques() result(array)
        ! setup the classic jacques test geometry

            use sdfs, only : container, box
            use vector_class

            implicit none

            type(container)         :: array(2)
            type(box), target, save :: medium, bbox
            real(kind=wp)           :: hgg

            hgg = 0.9_wp

            bbox = box(vector(10._wp, 10._wp, 2.001_wp), 0._wp, 0._wp, 0._wp, 1._wp, 2) 
            !420
            medium = box(vector(10._wp, 10._wp, 2._wp), 82._wp/(1._wp-hgg), 1.8_wp, hgg, 1.38_wp, 1)
            !630
            ! medium = box(vector(10., 10., 2.), 21./(1.-hgg), .23, hgg, 1.38, 1)


            allocate(array(1)%p, source=medium)
            allocate(array(2)%p, source=bbox)

            array(1)%p => medium
            array(2)%p => bbox

        end function setup_jacques

        function setup_sphere() result(array)
        ! setup the sphere test case from tran and jacques paper.
        !
        !
            use sdfs,         only : sphere, box, container, translate
            use vector_class, only : vector, invert

            implicit none
            
            type(container), allocatable :: array(:)
            type(sphere),   target, save :: sph
            type(box),      target, save :: bbox(2)

            real(kind=wp) :: mus, mua, n, hgg, t(4, 4)
            type(vector)  :: a
            
            mus = 0._wp; mua = 1.e-17_wp; hgg = 0._wp; n = 1._wp;
            bbox(1) = box(2._wp, mus, mua, hgg, n, 2)
            bbox(2) = box(2.01_wp, mus, 10000000._wp, hgg, n, 3)

            mus = 0._wp; mua = 1.e-17_wp; hgg = 0._wp; n = 1.33_wp;
            a = vector(.0_wp, 0._wp, 0._wp)
            t = invert(translate(a))
            sph = sphere(0.5_wp, mus, mua, hgg, n, 1, transform=t)

            allocate(array(3))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox(1))
            allocate(array(3)%p, source=bbox(2))

            array(1)%p => sph
            array(2)%p => bbox(1)
            array(3)%p => bbox(2)

        end function setup_sphere


        function interior_test() result(array)

            use sdfs, only : container, model, union, sphere, box, translate, model_init
            
            use vector_class

            implicit none

            type(container) :: array(1), cnta(3)
            type(sphere), target, save :: sph
            type(box),    target, save :: boxes(2), bbox
            type(model),  target, save :: m

            real(kind=wp) :: t(4, 4)

            ! packet = photon("point")

            bbox = box(2._wp, 0._wp, 0._wp, 0._wp, 1._wp, 2)

            t = invert(translate(vector(.25_wp-.2_wp, -.1_wp, 0._wp)))
            boxes(1) = box(vector(.5_wp, 1._wp, 1._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(0._wp-.2_wp, .6_wp, 0._wp)))
            boxes(2) = box(vector(1._wp, .5_wp, 1._wp), -.1_wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(.5_wp-.2_wp, -.1_wp, 0._wp)))
            sph = sphere(0.5_wp, 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            allocate(cnta(1)%p, source=sph)
            allocate(cnta(2)%p, source=boxes(1))
            allocate(cnta(3)%p, source=boxes(2))

            allocate(array(1)%p, source=m)

            cnta(1)%p => sph
            cnta(2)%p => boxes(1)
            cnta(3)%p => boxes(2)

            m = model_init(cnta, union)

            array(1)%p => m
            ! array(2)%p => bbox

        end function interior_test


        function exterior_test() result(array)

            use sdfs, only : container, moon, box, translate
            
            use vector_class

            implicit none

            type(container) :: array(1)
            ! type(sphere), target, save :: sph
            type(box),    target, save :: boxes(6), bbox
            ! type(model),  target, save :: m
            type(moon), target, save :: luna
            ! type(extrude), target, save :: ex
            real(kind=wp) :: t(4, 4)
            ! integer :: i

            ! packet = photon("point")
! d, ra, rb, mus, mua, hgg, n, layer, transform
            t = invert(translate(vector(-0.4_wp-1.0_wp, .3_wp, 0._wp)))

            luna = moon(1.2_wp+cos(3.9_wp), 1._wp, .8_wp, 0._wp, 0._wp, 0._wp, 0._wp, 1)
            ! ex = extrude(luna, 5.)

            bbox = box(2._wp, 0._wp, 0._wp, 0._wp, 1._wp, 2)
            t = invert(translate(vector(0.0_wp, 1.0_wp, 0._wp)))
            boxes(1) = box(2._wp*vector(2._wp, 0.2_wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(1.2_wp, 1.0_wp, 0._wp)))
            boxes(2) = box(2._wp*vector(0.8_wp, 1._wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(1.4_wp, -0.3_wp, 0._wp)))
            boxes(3) = box(2._wp*vector(0.6_wp, 0.9_wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(0.0_wp, -1.0_wp, 0._wp)))
            boxes(4) = box(2._wp*vector(1.0_wp, 0.2_wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)

            t = invert(translate(vector(-1.2_wp, -0.8_wp, 0._wp)))
            boxes(5) = box(2._wp*vector(0.8_wp, 0.6_wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)
            
            t = invert(translate(vector(-1.5_wp, 0.3_wp, 0._wp)))
            boxes(6) = box(2._wp*vector(0.6_wp, 0.7_wp, 0._wp), 0._wp, 0._wp, 0._wp, 1._wp, 1, transform=t)


            ! allocate(array(7))
            ! do i = 1, 6
            !     allocate(array(i)%p, source=boxes(i))
            !     array(i)%p => boxes(i)
            ! end do
            allocate(array(1)%p, source=luna)
            array(1)%p => luna

            ! allocate(array(7)%p, source=bbox)
            ! array(7)%p => bbox

        end function exterior_test


        function setup_exp(dict) result(array)
        ! Setup experimental geometry from Georgies paper. i.e a glass bottle with contents
        !
        !
            use sdfs,  only : container, box, cylinder, rotate_y, subtraction, translate
            use utils, only : deg2rad

            use vector_class, only : vector, invert
            use fhash,        only : fhash_tbl_t, key=>fhash_key


            implicit none

            type(fhash_tbl_t), intent(IN)  :: dict

            type(container), allocatable :: array(:)
            ! type(cylinder), target, save :: cyl(2)
            type(box),      target, save :: bbox, slab

            type(vector)  :: a, b
            real(kind=wp) :: n, t(4, 4), optprop(5)

            ! packet = photon("annulus")

            call dict%get(key("musb"), optprop(1))
            call dict%get(key("muab"), optprop(2))
            call dict%get(key("musc"), optprop(3))
            call dict%get(key("muac"), optprop(4))
            call dict%get(key("hgg"), optprop(5))
            n = 1._wp

            a = vector(-10._wp, 0._wp, 0._wp)
            b = vector(10._wp, 0._wp, 0._wp)
            !bottle
            ! cyl(2) = cylinder(a, b, 1.75, optprop(1), optprop(2), optprop(5), 1.5, 2)
            ! contents
            ! cyl(1) = cylinder(a, b, 1.55, optprop(3), optprop(4), optprop(5), 1.3, 1)

            t = invert(translate(vector(0._wp, 0._wp, -5._wp+1.75_wp)))
            slab = box(vector(10._wp, 10._wp, 10._wp), optprop(3), optprop(4), optprop(5), 1.3_wp, 1, transform=t)
            bbox = box(4._wp, 0.0_wp, 0.0_wp, 0.0_wp, n, 2)

            allocate(array(2))
            ! allocate(array(1)%p, source=cyl(1))
            ! allocate(array(2)%p, source=cyl(2))
            allocate(array(1)%p, source=slab)
            allocate(array(2)%p, source=bbox)

            ! array(1)%p => cyl(1)
            ! array(2)%p => cyl(2)
            ! array(3)%p => bbox
            array(1)%p => slab
            array(2)%p => bbox
        end function setup_exp

        function setup_scat_test(dict) result(array)

            use sdfs,  only : container, sphere, box
            use fhash, only : fhash_tbl_t, key=>fhash_key

            use vector_class

            implicit none

            type(fhash_tbl_t), intent(IN)  :: dict
            type(container), allocatable :: array(:)

            type(sphere), target, save :: sph
            type(box),    target, save :: bbox

            real(kind=wp) :: mus, mua, hgg, n, tau

            call dict%get(key("tau"), tau)

            n = 1._wp
            hgg = 0._wp
            mua = 1e-17_wp
            mus = tau

            sph = sphere(1._wp, mus, 1.5_wp, hgg, n, 1)
            bbox = box(2._wp, 0._wp, mua, hgg, n, 2)

            allocate(array(2))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox)

            array(1)%p => sph
            array(2)%p => bbox

        end function setup_scat_test

        function setup_omg_sdf() result(array)
            
            use sdfs,      only : container, cylinder, torus, model, box, union, rotate_y, model_init, translate
            
            use vector_class, only : vector, invert

            implicit none

            type(container), allocatable :: array(:)

            type(container) :: cnta(10)
            type(cylinder), target, save :: m(4), g(5)!save needed as these
            type(torus),    target, save :: sph       !are deallocated on sub exit
            type(box),      target, save :: boxy      !should be safe as this sub only called once
            type(model),    target, save :: omg_sdf
            
            type(vector)  :: a, b
            real(kind=wp) :: t(4, 4), mus, mua, hgg, n
            integer       :: j, layer

            mus = 10._wp
            mua = 0.16_wp
            hgg = 0.0_wp
            n = 2.65_wp
            layer = 1_wp

            ! x
            ! |
            ! |
            ! |
            ! |
            ! |_____z

            !O letter
            a = vector(0._wp, 0._wp, -0.7_wp)
            t = invert(translate(a))
            sph = torus(.2_wp, 0.05_wp, mus, mua, hgg, n, layer, transform=t)

            !M letter
            a = vector(-.25_wp, 0._wp, -.25_wp)
            b = vector(-.25_wp, 0._wp, .25_wp)
            t = invert(rotate_y(90._wp))
            m(1) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer, transform=t)
            
            a = vector(-.25_wp, 0._wp, -.25_wp)
            b = vector(.25_wp, 0._wp, .0_wp)
            m(2) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)
            
            a = vector(.25_wp, 0._wp, .0_wp)
            b = vector(-.25_wp, 0._wp, .25_wp)
            m(3) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            a = vector(-.25_wp, 0._wp, .25_wp)
            b = vector(.25_wp, 0._wp, .25_wp)
            m(4) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            !G letter
            a = vector(-.25_wp, 0._wp, .5_wp)
            b = vector(.25_wp, 0._wp, .5_wp)
            g(1) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            a = vector(-.25_wp, 0._wp, .5_wp)
            b = vector(-.25_wp, 0._wp, .75_wp)
            g(2) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            a = vector(.25_wp, 0._wp, .5_wp)
            b = vector(.25_wp, 0._wp, .75_wp)
            g(3) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            a = vector(.25_wp, 0._wp, .75_wp)
            b = vector(0._wp, 0._wp, .75_wp)
            g(4) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            a = vector(0._wp, 0._wp, .625_wp)
            b = vector(0._wp, 0._wp, .75_wp)
            g(5) = cylinder(a, b, .05_wp, mus, mua, hgg, n, layer)

            !bbox
            boxy = box(2._wp, 0._wp, 0._wp, 0._wp, 1.0_wp, 2)

            allocate(array(2))
            allocate(array(1)%p, source=omg_sdf)
            allocate(array(2)%p, source=boxy)

            do j = 1, size(m)
                allocate(cnta(j)%p, source=m(j))
            end do

            do j = 1, size(g)
                allocate(cnta(j+size(m))%p, source=g(j))
            end do

            cnta(1)%p => m(1)
            cnta(2)%p => m(2)
            cnta(3)%p => m(3)
            cnta(4)%p => m(4)
            cnta(5)%p => sph
            cnta(6)%p => g(1)
            cnta(7)%p => g(2)
            cnta(8)%p => g(3)
            cnta(9)%p => g(4)
            cnta(10)%p => g(5)

            omg_sdf = model_init(cnta, union, 0.05_wp)

            array(1)%p => omg_sdf ! model
            array(2)%p => boxy    ! bbox

        end function setup_omg_sdf


        function get_vessels() result(array)

            use sdfs, only : container, model, capsule, model_init, union, box, union
            use vector_class, only : vector, invert

            implicit none

            type(container), allocatable :: array(:), cnta(:)
            type(model), target, save :: vessels
            type(capsule), allocatable, target, save :: cyls(:)
            type(box), target, save :: bbox

            real(kind=wp), allocatable :: nodes(:, :), radii(:)
            integer, allocatable :: edges(:, :)
            integer :: io, edge_cnt, tmp1, tmp2, u, node_cnt, i, counter
            real(kind=wp) :: x, y, z, radius, res, maxx, maxy, maxz
            real(kind=wp) :: musv, muav, gv, nv
            real(kind=wp) :: musd, muad, gd, nd
            type(vector) :: a, b

            ! packet = photon("uniform")
            ! mua, mus, g, n]
            !MCmatlab: an open-source, user-friendly, MATLAB-integrated three-dimensional Monte Carlo light transport solver with heat diffusion and tissue damage
            muav = 231._wp
            musv = 94._wp
            gv = 0.9_wp
            nv = 1.37_wp

            muad = 0.458_wp
            musd = 357._wp
            gd = 0.9_wp
            nd = 1.37_wp

            !get number of edges
            open(newunit=u, file="res/edges.dat", iostat=io)
            edge_cnt = 0
            do
                read(u,*,iostat=io)tmp1, tmp2
                if(io /= 0)exit
                edge_cnt = edge_cnt + 1
            end do
            close(u)

            !get number of nodes and radii
            open(newunit=u, file="res/nodes.dat", iostat=io)
            node_cnt = 0
            do
                read(u,*,iostat=io)x, y, z
                if(io /= 0)exit
                node_cnt = node_cnt + 1
            end do
            allocate(edges(edge_cnt, 2), nodes(node_cnt, 3), radii(node_cnt))
            ! print*,node_cnt,edge_cnt
            ! stop
            !read in edges
            open(newunit=u, file="res/edges.dat", iostat=io)
            do i = 1, edge_cnt
                read(u,*,iostat=io)edges(i, :)
                if(io /= 0)exit
            end do
            close(u)

            !read in nodes
            open(newunit=u, file="res/nodes.dat", iostat=io)
            do i = 1, edge_cnt
                read(u,*,iostat=io)nodes(i, :)
                if(io /= 0)exit
            end do
            close(u)

            !read in radii
            open(newunit=u, file="res/radii.dat", iostat=io)
            do i = 1, node_cnt
                read(u,*,iostat=io)radii(i)
                if(io /= 0)exit
            end do
            close(u)

            res = 0.001_wp!0.01mm
            maxx = maxval(abs(nodes(:, 1)))
            maxy = maxval(abs(nodes(:, 2)))
            maxz = maxval(abs(nodes(:, 3)))

            nodes(:, 1) = (nodes(:, 1) / maxx) - 0.5_wp
            nodes(:, 2) = (nodes(:, 2) / maxy) - 0.5_wp
            nodes(:, 3) = (nodes(:, 3) / maxz) - 0.5_wp
            nodes(:, 1) = nodes(:, 1) * maxx * res
            nodes(:, 2) = nodes(:, 2) * maxy * res
            nodes(:, 3) = nodes(:, 3) * maxz * res

            !allocate SDFs and container
            allocate(cyls(edge_cnt), cnta(edge_cnt))
            do i = 1, edge_cnt
                allocate(cnta(i)%p, source=cyls(1))
            end do

            counter = 1
            do i = 1, edge_cnt
                a = vector(nodes(edges(i, 1), 1), nodes(edges(i, 1), 2), nodes(edges(i, 1), 3))
                b = vector(nodes(edges(i, 2), 1), nodes(edges(i, 2), 2), nodes(edges(i, 2), 3))
                radius = radii(edges(i, 1)) * res

                                                     ! mus, mua, hgg, n, layer
                cyls(counter) = capsule(a, b, radius, musv, muav, gv, nv, 1)
                cnta(counter)%p => cyls(counter)
                counter = counter + 1
            end do

            allocate(array(2))
            allocate(array(1)%p, source=vessels)
            vessels = model_init(cnta, union, 0.005_wp)

            bbox = box(vector(.32_wp, .18_wp, .26_wp), musd, muad, gd, nd, 2)
            allocate(array(2)%p, source=bbox)

            array(1)%p => vessels
            array(2)%p => bbox

        end function get_vessels


        function setup_skin_model() result(array)

            use sdfs, only : container, model, capsule, box, translate, model_init, smoothunion

            use vector_class, only : vector, invert

            implicit none

            type(container), allocatable :: array(:), cnta(:)

            type(capsule), target, save :: cyls(15)
            type(box),      target, save :: skin(3)
            type(model),    target, save :: vessels

            integer       :: i
            real(kind=wp) :: mus, mua, hgg, n, t(4, 4)
            real(kind=wp) :: mus_epi, mus_derm, mua_epi, mua_derm, n_epi, n_derm, hgg_epi, hgg_derm
            real(kind=wp) :: mus_b, mua_b, hgg_b, n_b
            type(vector)  :: a, b, c

            n = 1.e0_wp
            hgg = 1.e0_wp
            mua = .00036_wp
            mus = 10._wp

            mus_epi = 376._wp
            mua_epi = 16.6_wp
            hgg_epi = 0.9_wp
            n_epi = 1._wp

            mus_derm = 357._wp
            mua_derm = 0.459_wp
            hgg_derm = 0.9_wp
            n_derm = 1._wp

            mus_b = 94._wp
            mua_b = 231000._wp
            hgg_b = 0.9_wp
            n_b = 1._wp

            ! total 0.1 cm 
            c = vector(0._wp, 0._wp, 0.045_wp)
            t = invert(translate(c))
            skin(3) = box(vector(.1_wp, .1_wp, .01_wp), mus, mua, hgg, n, 4, transform=t)!water

            c = vector(0._wp, 0._wp, .037_wp)
            t = invert(translate(c))
            skin(2) = box(vector(.1_wp, .1_wp, .006_wp), mus_epi, mua_epi, hgg_epi, n_epi, 3, transform=t)!epidermis

            c = vector(0._wp, 0._wp, -.008_wp)
            t = invert(translate(c))
            skin(1) = box(vector(.1_wp, .1_wp, 0.084_wp), mus_derm, mua_derm, hgg_derm, n_derm, 2, transform=t)!dermis
            ! skin(1) = box(vector(.1, .1, .1), mus, mua, hgg, n, 2)!bbox for debug
            

            !   x
            !   |
            !   |
            !   |
            !   |
            !   |_______y

            a = vector(0._wp, -0.05_wp, 0._wp)
            b = vector(0._wp, -0.03_wp, 0._wp)
            cyls(1) = capsule(a, b, .005_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.005_wp, -0.025_wp, 0._wp)
            b = vector(0.03_wp, -0.02_wp, 0._wp)
            cyls(2) = capsule(a, b, .005_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.03_wp, -0.02_wp, 0._wp)
            b = vector(0._wp, 0.01_wp, -0.02_wp)
            cyls(3) = capsule(a, b, .005_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0._wp, 0.01_wp, -0.02_wp)
            b = vector(0._wp, 0.03_wp, 0._wp)
            cyls(4) = capsule(a, b, .005_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0._wp, 0.03_wp, 0._wp)
            b = vector(0.03_wp, 0.05_wp, 0._wp)
            cyls(5) = capsule(a, b, .005_wp, mus_b, mua_b, hgg_b, n_b, 1)

            !branch 1
            a = vector(-0.0025_wp, -0.04_wp, 0._wp)
            b = vector(-0.03_wp, -0.03_wp, 0._wp)
            cyls(6) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03_wp, -0.028_wp, 0._wp)
            b = vector(-0.03_wp, -0.02_wp, 0._wp)
            cyls(7) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03_wp, -0.018_wp, 0._wp)
            b = vector(0.0025_wp, 0.011_wp, -0.02_wp)
            cyls(8) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)


            !branch 2
            a = vector(0.015_wp, -0.025_wp, 0._wp)
            b = vector(0.048_wp, -0.02_wp, 0._wp)
            cyls(9) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.048_wp, -0.018_wp, 0._wp)
            b = vector(0.04_wp, 0.01_wp, 0._wp)
            cyls(10) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.04_wp, 0.012_wp, 0._wp)
            b = vector(0.03_wp, 0.025_wp, 0._wp)
            cyls(11) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.028_wp, 0.025_wp, 0._wp)
            b = vector(0.00_wp, 0.025_wp, 0._wp)
            cyls(12) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            !branch 3
            a = vector(0.015_wp, -0.025_wp, 0._wp)
            b = vector(0.015_wp, 0.00_wp, 0.02_wp)
            cyls(13) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.015_wp, 0.002_wp, 0.02_wp)
            b = vector(0.02_wp, 0.025_wp, 0.01_wp)
            cyls(14) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.02_wp, 0.027_wp, 0.01_wp)
            b = vector(0.015_wp, 0.04_wp, 0._wp)
            cyls(15) = capsule(a, b, .001_wp, mus_b, mua_b, hgg_b, n_b, 1)


            allocate(array(1))

            allocate(array(1)%p, source=vessels)
            ! allocate(array(2)%p, source=skin(1))
            ! allocate(array(3)%p, source=skin(2))
            ! allocate(array(4)%p, source=skin(3))

            allocate(cnta(size(cyls)))
            do i = 1, size(cnta)
                allocate(cnta(i)%p, source=cyls(i))
                cnta(i)%p => cyls(i)
            end do

            vessels = model_init(cnta, smoothunion, 0.005_wp)
            array(1)%p => vessels
            array(2)%p => skin(1)
            array(3)%p => skin(2)
            array(4)%p => skin(3)

        end function setup_skin_model

        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            character(len=:), allocatable :: mkdirCMD
            logical :: dirExists

            !get current working directory
            call get_environment_variable('PWD', cwd)
  
            ! get 'home' dir from cwd
            homedir = trim(cwd)
            ! get data dir
            fileplace = trim(homedir)//'/data/'
            !check if data directory exists. if not create it
#ifdef __GFORTRAN__
            inquire(file=trim(fileplace)//"/.", exist=dirExists)
#elif __INTEL_COMPILER
            inquire(directory=trim(fileplace), exist=dirExists)
#endif
            if(.not. dirExists)then
                mkdirCMD = "mkdir -p "//trim(fileplace)
                call execute_command_line(mkdirCMD)
                mkdirCMD = "mkdir -p "//trim(fileplace)//"jmean/"
                call execute_command_line(mkdirCMD)
            end if

            ! get res dir
            resdir = trim(homedir)//'/res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none

            jmean = 0._wp
            jmeanGLOBAL = 0._wp

        end subroutine zarray


        subroutine alloc_array(nxg, nyg, nzg)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iarray

            implicit none

            integer, intent(IN) :: nxg, nyg, nzg

            allocate(jmean(nxg, nyg, nzg), jmeanGLOBAL(nxg, nyg, nzg))

        end subroutine alloc_array
end module subs