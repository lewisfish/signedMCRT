module subs

    use dict_mod

    implicit none

    type :: settings_t
        integer :: nphotons, iseed
        character(len=:), allocatable :: experiment, outfile, renderfile
        logical :: render_bool       
    end type settings_t

    public  :: setup_simulation, print_time, get_time, Sellmeier
    private :: directory, alloc_array, zarray

    contains

        subroutine setup_simulation(grid, packet, sdfarray, settings, dict)
        ! Read in parameters
            use constants, only : resdir
            use gridMod,   only : cart_grid
            use sdfs,      only : container
            
            use photonMod
            use ch_opt
            use vector_class

            implicit none

            type(dict_t),    optional,    intent(IN)  :: dict
            type(container), allocatable, intent(OUT) :: sdfarray(:)
            type(cart_grid),              intent(OUT) :: grid
            type(photon),                 intent(OUT) :: packet
            type(settings_t),             intent(OUT) :: settings

            real    :: xmax, ymax, zmax, albedo(3), hgg(3), g2(3), kappa(3)
            integer :: nxg, nyg, nzg, u
            character(len=256) :: experiment, outfile, renderfile

            !set directory paths
            call directory

            open(newunit=u,file=trim(resdir)//'input.params',status='old')
                read(u,*) settings%nphotons
                
                read(u,"(A)") experiment
                settings%experiment = trim(experiment)
                
                read(u,"(A)") outfile
                settings%outfile = trim(outfile)
                
                read(u,"(A)") renderfile
                settings%renderfile = trim(renderfile)
                
                read(u,*) settings%render_bool
                read(u,*) settings%iseed
                read(u,*) xmax
                read(u,*) ymax
                read(u,*) zmax
                read(u,*) nxg
                read(u,*) nyg
                read(u,*) nzg
            close(u)

            !allocate and set arrays to 0
            call alloc_array(nxg, nyg, nzg)
            call zarray

            ! Set up grid
            grid = cart_grid(nxg, nyg, nzg, xmax, ymax, zmax)
            call init_opt1(kappa, albedo, hgg, g2)

            select case(settings%experiment)
                case("logo")
                    sdfarray = setup_logo(packet)
                case("neural")
                    sdfarray = setup_neural_sdf(packet)
                case("omg")
                    sdfarray = setup_omg_sdf(packet)
                case("scat_test")
                    sdfarray = setup_scat_test(packet, dict)
                case("fresnel_test")
                    sdfarray = setup_fresnel_test(packet)
                case("skin")
                    sdfarray = setup_skin_model(packet)
                case("interior")
                    sdfarray = interior_test(packet)
                case("exterior")
                    sdfarray = exterior_test(packet)
                case("sphere")
                    sdfarray = setup_sphere(packet)
                case("exp")
                    sdfarray = setup_exp(packet, dict)
                case("jacques")
                    sdfarray = setup_jacques(packet)
                case("vessels")
                    sdfarray = get_vessels(packet)
                case("lens")
                    sdfarray = lens_test_setup(packet)
                case default
                    error stop "no such routine"
            end select

        end subroutine setup_simulation


        function lens_test_setup(packet) result(array)

            use sdfs, only : sphere, container, box, model, intersection, model_init, translate
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container),target, save :: cnta(2)
            type(container) :: array(2)
            type(box),     target, save :: bbox
            type(sphere), target, save :: sph(2)
            type(model), target, save :: m
            
            type(vector) :: a, b
            real :: hgg, mus, mua, n1, n2, t(4,4)
            integer :: layer, i

            packet = photon("uniform")

            mus = 0.d0
            mua = 1.d-17
            hgg = 0.0
            n1 = 1.33
            n2 = 1.52
            layer = 1

            bbox = box(vector(2., 2., 2.), mus, mua, hgg, n1, 2) 
            !sphere1
            a = vector(0., 0., 0.05)
            t = 0.d0
            t = invert(translate(a))
            sph(1) = sphere(.2, mus, mua, hgg, n2, 1, transform=t)
            !sphere2
            b = vector(0., 0., -0.15)
            t = 0.d0
            t = invert(translate(b))
            sph(2) = sphere(.25, mus, mua, hgg, n2, 1, transform=t)

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


        function setup_logo(packet) result(array)

            use sdfs, only : container, box, segment, extrude, model, union, model_init
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container), allocatable :: cnta(:), array(:)
            type(box),     target, save :: bbox
            type(segment), allocatable, target, save :: seg(:)
            type(extrude), allocatable, target, save :: ex(:)
            type(model), target, save :: m
            
            real :: hgg, mus, mua, n
            integer :: layer, i

            packet = photon("uniform")

            allocate(array(2), cnta(725), seg(725), ex(725))

            mus = 10.
            mua = .1
            hgg = 0.9
            n = 1.5
            layer = 1

            bbox = box(vector(10., 10., 2.001), 0., 0., 0., 1., 2) 
            !420

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


        function setup_neural_sdf(packet) result(array)

            use sdfs, only : container, box, neural
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container)         :: array(2)
            type(box), target, save :: bbox
            type(neural), target, save :: neu
            
            real :: hgg

            packet = photon("point")

            hgg = 0.9

            bbox = box(vector(10., 10., 2.001), 0., 0., 0., 1., 2) 
            !420
            neu = neural(82./(1.-hgg), 1.8, hgg, 1.38, 1)

            allocate(array(1)%p, source=neu)
            allocate(array(2)%p, source=bbox)

            array(1)%p => neu
            array(2)%p => bbox

        end function setup_neural_sdf


        function setup_jacques(packet) result(array)

            use sdfs, only : container, box
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container)         :: array(2)
            type(box), target, save :: medium, bbox
            real                 :: hgg

            packet = photon("uniform")

            hgg = 0.9

            bbox = box(vector(10., 10., 2.001), 0., 0., 0., 1., 2) 
            !420
            medium = box(vector(10., 10., 2.), 82./(1.-hgg), 1.8, hgg, 1.38, 1)
            !630
            ! medium = box(vector(10., 10., 2.), 21./(1.-hgg), .23, hgg, 1.38, 1)


            allocate(array(1)%p, source=medium)
            allocate(array(2)%p, source=bbox)

            array(1)%p => medium
            array(2)%p => bbox

        end function setup_jacques

        function setup_sphere(packet) result(array)

            use sdfs, only : sphere, box, container, translate
            use vector_class
            use photonMod

            implicit none
            
            type(photon), intent(OUT) :: packet

            type(container), allocatable :: array(:)
            type(sphere),   target, save :: sph
            type(box),      target, save :: bbox(2)

            real :: mus, mua, n, hgg, t(4, 4)
            type(vector) :: a

            packet = photon("uniform")
            
            mus = 0.; mua = 1.d-17; hgg = 0.; n = 1.;
            bbox(1) = box(2., mus, mua, hgg, n, 2)
            bbox(2) = box(2.01, mus, 10000000.d0, hgg, n, 2)

            mus = 0.; mua = 1.d-17; hgg = .0; n = 1.46;
            a = vector(.5, 0., -1.)
            t = invert(translate(a))
            sph = sphere(1., mus, mua, hgg, n, 1, transform=t)

            allocate(array(3))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox(1))
            allocate(array(3)%p, source=bbox(2))

            array(1)%p => sph
            array(2)%p => bbox(1)
            array(3)%p => bbox(2)

        end function setup_sphere


        function interior_test(packet) result(array)

            use sdfs, only : container, model, union, sphere, box, translate, model_init
            
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container) :: array(1), cnta(3)
            type(sphere), target, save :: sph
            type(box),    target, save :: boxes(2), bbox
            type(model),  target, save :: m

            real :: t(4, 4)

            packet = photon("point")

            bbox = box(2., 0., 0., 0., 1., 2)

            t = invert(translate(vector(.25-.2, -.1, 0.)))
            boxes(1) = box(vector(.5, 1., 1.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(0.-.2, .6, 0.)))
            boxes(2) = box(vector(1., .5, 1.), -.1, 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(.5-.2, -.1, 0.)))
            sph = sphere(0.5, 0., 0., 0., 1., 1, transform=t)

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


        function exterior_test(packet) result(array)

            use sdfs, only : container, moon, box, translate
            
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container) :: array(1)
            ! type(sphere), target, save :: sph
            type(box),    target, save :: boxes(6), bbox
            ! type(model),  target, save :: m
            type(moon), target, save :: luna
            ! type(extrude), target, save :: ex
            real :: t(4, 4)
            ! integer :: i

            packet = photon("point")
! d, ra, rb, mus, mua, hgg, n, layer, transform
            t = invert(translate(vector(-0.4-1.0, .3, 0.)))

            luna = moon(1.2+cos(3.9), 1., .8, 0., 0., 0., 0., 1)
            ! ex = extrude(luna, 5.)

            bbox = box(2., 0., 0., 0., 1., 2)
            t = invert(translate(vector(0.0, 1.0, 0.)))
            boxes(1) = box(2.*vector(2., 0.2, 0.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(1.2, 1.0, 0.)))
            boxes(2) = box(2.*vector(0.8, 1., 0.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(1.4, -0.3, 0.)))
            boxes(3) = box(2.*vector(0.6, 0.9, 0.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(0.0, -1.0, 0.)))
            boxes(4) = box(2.*vector(1.0, 0.2, 0.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(-1.2, -0.8, 0.)))
            boxes(5) = box(2.*vector(0.8, 0.6, 0.), 0., 0., 0., 1., 1, transform=t)
            
            t = invert(translate(vector(-1.5, 0.3, 0.)))
            boxes(6) = box(2.*vector(0.6, 0.7, 0.), 0., 0., 0., 1., 1, transform=t)


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


        function setup_exp(packet, dict) result(array)
            
            use sdfs,  only : container, box, cylinder, rotate_y, subtraction, translate
            use utils, only : deg2rad

            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet
            type(dict_t), intent(IN)  :: dict

            type(container), allocatable :: array(:)
            ! type(cylinder), target, save :: cyl(2)
            type(box),      target, save :: bbox, slab

            type(vector) :: a, b
            real         :: n, t(4, 4), optprop(5)

            packet = photon("annulus")

            optprop(1) = dict%get_value_real("musb")
            optprop(2) = dict%get_value_real("muab")
            optprop(3) = dict%get_value_real("musc")
            optprop(4) = dict%get_value_real("muac")
            optprop(5) = dict%get_value_real("hgg")

            n = 1.d0

            a = vector(-10., 0., 0.)
            b = vector(10., 0., 0.)
            !bottle
            ! cyl(2) = cylinder(a, b, 1.75, optprop(1), optprop(2), optprop(5), 1.5, 2)
            ! contents
            ! cyl(1) = cylinder(a, b, 1.55, optprop(3), optprop(4), optprop(5), 1.3, 1)

            t = invert(translate(vector(0., 0., -5.+1.75)))
            slab = box(vector(10., 10., 10.), optprop(3), optprop(4), optprop(5), 1.3, 1, transform=t)
            bbox = box(4., 0.0, 0.0, 0.0, n, 2)

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


        function setup_fresnel_test(packet) result(array)

            use sdfs, only : container, box
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container), allocatable :: array(:)
            type(box), target, save :: bbox, fibre

            integer :: layer
            real    :: mus, mua, hgg, n

            packet = photon("uniform")

            layer = 1
            n = 1.d0
            hgg = 0.d0
            mua = .1
            mus = 0.d0

            fibre = box(1., mus, mua, hgg, 1.5, layer)
            bbox  = box(2., mus, mua, hgg, n, 2)

            allocate(array(2))
            allocate(array(1)%p, source=fibre)
            allocate(array(2)%p, source=bbox)

            array(1)%p => fibre
            array(2)%p => bbox

        end function setup_fresnel_test

        function setup_scat_test(packet, dict) result(array)

            use sdfs, only : container, sphere, box
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet
            type(container), allocatable :: array(:)
            type(dict_t), intent(IN) :: dict

            type(sphere), target, save :: sph
            type(box), target, save :: bbox

            real :: mus, mua, hgg, n, tau

            tau = dict%get_value_real("tau")
            packet = photon("point")

            n = 1.d0
            hgg = 0.d0
            mua = 1d-17
            mus = tau

            sph = sphere(1., mus, mua, hgg, n, 1)
            bbox = box(2., 0.d0, mua, hgg, n, 2)

            allocate(array(2))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox)

            array(1)%p => sph
            array(2)%p => bbox

        end function setup_scat_test

        function setup_omg_sdf(packet) result(array)
            
            use sdfs,      only : container, cylinder, torus, model, box, smoothunion, rotate_y, model_init, translate
            
            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container), allocatable :: array(:)

            type(container) :: cnta(10)
            type(cylinder), target, save :: m(4), g(5)!save needed as these
            type(torus),    target, save :: sph       !are deallocated on sub exit
            type(box),      target, save :: boxy      !should be safe as this sub only called once
            type(model),    target, save :: omg_sdf
            
            type(vector) :: a, b
            real         :: t(4, 4), mus, mua, hgg, n
            integer      :: j, layer

            packet = photon("uniform")

            mus = 10.
            mua = 0.16
            hgg = 0.0d0
            n = 2.65
            layer = 1

            ! x
            ! |
            ! |
            ! |
            ! |
            ! |_____z

            !O letter
            a = vector(0., 0., -0.7)
            t = invert(translate(a))
            sph = torus(.2, 0.05, mus, mua, hgg, n, layer, transform=t)

            !M letter
            a = vector(-.25, 0., -.25)
            b = vector(-.25, 0., .25)
            t = invert(rotate_y(90.))
            m(1) = cylinder(a, b, .05, mus, mua, hgg, n, layer, transform=t)
            
            a = vector(-.25, 0., -.25)
            b = vector(.25, 0., .0)
            m(2) = cylinder(a, b, .05, mus, mua, hgg, n, layer)
            
            a = vector(.25, 0., .0)
            b = vector(-.25, 0., .25)
            m(3) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            a = vector(-.25, 0., .25)
            b = vector(.25, 0., .25)
            m(4) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            !G letter
            a = vector(-.25, 0., .5)
            b = vector(.25, 0., .5)
            g(1) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            a = vector(-.25, 0., .5)
            b = vector(-.25, 0., .75)
            g(2) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            a = vector(.25, 0., .5)
            b = vector(.25, 0., .75)
            g(3) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            a = vector(.25, 0., .75)
            b = vector(0., 0., .75)
            g(4) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            a = vector(0., 0., .625)
            b = vector(0., 0., .75)
            g(5) = cylinder(a, b, .05, mus, mua, hgg, n, layer)

            !bbox
            boxy = box(2., 0.d0, 0., 0.d0, 1.0, 2)

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

            omg_sdf = model_init(cnta, smoothunion, 0.05)

            array(1)%p => omg_sdf ! model
            array(2)%p => boxy    ! bbox

        end function setup_omg_sdf


        function get_vessels(packet) result(array)

            use sdfs, only : container, model, capsule, model_init, union, box, union
            use photonMod
            use vector_class

            implicit none

            type(photon) :: packet
            type(container), allocatable :: array(:), cnta(:)
            type(model), target, save :: vessels
            type(capsule), allocatable, target, save :: cyls(:)
            type(box), target, save :: bbox

            real, allocatable :: nodes(:, :), radii(:)
            integer, allocatable :: edges(:, :)
            integer :: io, edge_cnt, tmp1, tmp2, u, node_cnt, i, counter
            real :: x, y, z, radius, res, maxx, maxy, maxz
            real :: musv, muav, gv, nv
            real :: musd, muad, gd, nd
            type(vector) :: a, b

            packet = photon("uniform")
            ! mua, mus, g, n]
            !MCmatlab: an open-source, user-friendly, MATLAB-integrated three-dimensional Monte Carlo light transport solver with heat diffusion and tissue damage
            muav = 231.
            musv = 94.
            gv = 0.9
            nv = 1.37

            muad = 0.458
            musd = 357.
            gd = 0.9
            nd = 1.37

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

            res = 0.001!0.01mm
            maxx = maxval(abs(nodes(:, 1)))
            maxy = maxval(abs(nodes(:, 2)))
            maxz = maxval(abs(nodes(:, 3)))

            nodes(:, 1) = (nodes(:, 1) / maxx) - 0.5
            nodes(:, 2) = (nodes(:, 2) / maxy) - 0.5
            nodes(:, 3) = (nodes(:, 3) / maxz) - 0.5
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
            vessels = model_init(cnta, union, 0.005)

            bbox = box(vector(.32, .18,.26), musd, muad, gd, nd, 2)
            allocate(array(2)%p, source=bbox)

            array(1)%p => vessels
            array(2)%p => bbox

        end function get_vessels


        function setup_skin_model(packet) result(array)

            use sdfs, only : container, model, capsule, box, translate, model_init, smoothunion

            use vector_class
            use photonMod

            implicit none

            type(photon), intent(OUT) :: packet

            type(container), allocatable :: array(:), cnta(:)

            type(capsule), target, save :: cyls(15)
            type(box),      target, save :: skin(3)
            type(model),    target, save :: vessels

            integer :: i
            real    :: mus, mua, hgg, n, t(4, 4)

            real :: mus_epi, mus_derm, mua_epi, mua_derm, n_epi, n_derm, hgg_epi, hgg_derm
            real :: mus_b, mua_b, hgg_b, n_b
            type(vector) :: a, b, c

            packet = photon("uniform")

            n = 1.d0
            hgg = 1.d0
            mua = .00036
            mus = 10.

            mus_epi = 376.
            mua_epi = 16.6
            hgg_epi = 0.9
            n_epi = 1.

            mus_derm = 357.
            mua_derm = 0.459
            hgg_derm = 0.9
            n_derm = 1.

            mus_b = 94.
            mua_b = 231000.
            hgg_b = 0.9
            n_b = 1.

            ! total 0.1 cm 
            c = vector(0., 0., 0.045)
            t = invert(translate(c))
            skin(3) = box(vector(.1, .1, .01), mus, mua, hgg, n, 4, transform=t)!water

            c = vector(0., 0., .037)
            t = invert(translate(c))
            skin(2) = box(vector(.1, .1, .006), mus_epi, mua_epi, hgg_epi, n_epi, 3, transform=t)!epidermis

            c = vector(0., 0., -.008)
            t = invert(translate(c))
            skin(1) = box(vector(.1, .1, 0.084), mus_derm, mua_derm, hgg_derm, n_derm, 2, transform=t)!dermis
            ! skin(1) = box(vector(.1, .1, .1), mus, mua, hgg, n, 2)!bbox for debug
            

            !   x
            !   |
            !   |
            !   |
            !   |
            !   |_______y

            a = vector(0., -0.05, 0.)
            b = vector(0., -0.03, 0.)
            cyls(1) = capsule(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.005, -0.025, 0.)
            b = vector(0.03, -0.02, 0.)
            cyls(2) = capsule(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.03, -0.02, 0.)
            b = vector(0., 0.01, -0.02)
            cyls(3) = capsule(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0., 0.01, -0.02)
            b = vector(0., 0.03, 0.)
            cyls(4) = capsule(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0., 0.03, 0.)
            b = vector(0.03, 0.05, 0.)
            cyls(5) = capsule(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            !branch 1
            a = vector(-0.0025, -0.04, 0.)
            b = vector(-0.03, -0.03, 0.)
            cyls(6) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03, -0.028, 0.)
            b = vector(-0.03, -0.02, 0.)
            cyls(7) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03, -0.018, 0.)
            b = vector(0.0025, 0.011, -0.02)
            cyls(8) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)


            !branch 2
            a = vector(0.015, -0.025, 0.)
            b = vector(0.048, -0.02, 0.)
            cyls(9) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.048, -0.018, 0.)
            b = vector(0.04, 0.01, 0.)
            cyls(10) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.04, 0.012, 0.)
            b = vector(0.03, 0.025, 0.)
            cyls(11) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.028, 0.025, 0.)
            b = vector(0.00, 0.025, 0.)
            cyls(12) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            !branch 3
            a = vector(0.015, -0.025, 0.)
            b = vector(0.015, 0.00, 0.02)
            cyls(13) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.015, 0.002, 0.02)
            b = vector(0.02, 0.025, 0.01)
            cyls(14) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.02, 0.027, 0.01)
            b = vector(0.015, 0.04, 0.)
            cyls(15) = capsule(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)


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

            vessels = model_init(cnta, smoothunion, 0.005)
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

            !get current working directory

            call get_environment_variable('PWD', cwd)
  
            ! get 'home' dir from cwd
            homedir = trim(cwd)
            ! get data dir
            fileplace = trim(homedir)//'/data/'
            ! get res dir
            resdir = trim(homedir)//'/res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none

            jmean = 0.
            jmeanGLOBAL = 0.

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


        real function get_time()

#ifdef _OPENMP
            use omp_lib
#endif
            implicit none

#ifdef _OPENMP
                get_time = omp_get_wtime()
#else
                call cpu_time(get_time)
#endif

        end function get_time


        subroutine print_time(time, id)

            implicit none

            real,    intent(IN) :: time
            integer, intent(IN) :: id

            if(id == 0)then
                ! if(time >= 60. .and. id == 0)then
                !    print*, floor((time)/60.) + mod(time, 60.)/100.
                ! else
                   print*, 'time taken ~',time,'s'
                ! end if
            end if

        end subroutine print_time


        real function Sellmeier(wave)
        ! Sellmeier equation for fused quatrz
        ! I. H. Malitson. Interspecimen comparison of the refractive index of fused silica
            implicit none

            real, intent(IN) :: wave
            real :: wave2, a, b ,c

            wave2 = wave**2

            a = (0.6961663d0*wave2)/(wave2 - 0.0684043d0**2)
            b = (0.4079426*wave2) / (wave2 - .1162414**2)
            c = (0.8974794*wave2) / (wave2 - 9.896161**2)

            Sellmeier = sqrt(1.d0 + (a + b + c))

        end function Sellmeier


    subroutine makeImage(image, dir, pos, diameter)
        use vector_class
        implicit none

        integer,      intent(INOUT) :: image(-200:200,-200:200)
        real,         intent(IN)    :: diameter
        type(vector), intent(IN)    :: pos, dir

        real :: binwid, angle, na, bottom, top
        type(vector) :: n, d
        integer :: xp, yp

        n = vector(0., 0., -1.)
        n = n%magnitude()
        d = dir%magnitude()
        d = (-1.)*d

        top = n .dot. d
        bottom = sqrt(d .dot. d) * sqrt(n .dot. n)
        angle = acos(top / bottom)
        na = asin(0.22)

        if(angle > na)then
            return
        end if
        binwid = diameter / 401.
        xp = floor(pos%x / binwid)
        yp = floor(pos%y / binwid)
        if(abs(xp) > 200 .or. abs(yp) > 200)then
            return
        end if
!$omp atomic
        image(xp,yp) = image(xp,yp) + 1

    end subroutine makeImage

end MODULE subs