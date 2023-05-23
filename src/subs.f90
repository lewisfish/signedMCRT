module subs

    use constants, only : wp
    use tomlf

    implicit none

    private
    public  :: setup_simulation, dealloc_array, directory

    contains

        subroutine setup_simulation(sdfarray, dict)
        ! Read in parameters
        ! Setup up various simulation parameters and routines
        !
            use vector_class
            use sdfNew,        only : sdf
            use sim_state_mod, only : settings => state

            type(toml_table),  optional,  intent(INOUT) :: dict
            type(sdf), allocatable, intent(OUT)   :: sdfarray(:)

            !allocate and set arrays to 0
            call alloc_array(settings%grid%nxg, settings%grid%nyg, settings%grid%nzg)
            call zarray()

            ! setup geometry using SDFs
            select case(settings%experiment)
                case("logo")
                    sdfarray = setup_logo()
                case("omg")
                    sdfarray = setup_omg_sdf()
                case("scat_test")
                    sdfarray = setup_scat_test(dict)
                case("scat_test2")
                    sdfarray = setup_scat_test2(dict)
                case("aptran")
                    sdfarray = setup_sphere()
                case("vessels")
                    sdfarray = get_vessels()
                case("sphere_scene")
                    sdfarray = setup_sphere_scene(dict)
                case("test_egg")
                    sdfarray = setup_egg()
                case default
                    error stop "no such routine"
            end select

        end subroutine setup_simulation

        function setup_egg() result(array)

            use sdfNew, only : sdf, onion, sphere, box, revolution, egg
            use vector_class
            use opticalProperties

            type(sdf), allocatable :: array(:)
            type(box) :: bbox
            type(revolution), save :: albumen, rev1
            type(onion) :: shell
            type(sphere) :: yolk
            type(opticalProp_t) :: opt(4)
            type(egg), save :: egg_shell, egg_albumen

            real(kind=wp) :: r1, r2, h
            
            r1 = 3._wp
            r2 = 3._wp * sqrt(2._wp - sqrt(2._wp))
            h = r2
            
            !width = 42mm
            !height = 62mm

            !shell
            opt(1) = mono(100._wp, 10._wp, 0.0_wp, 1.37_wp)
            egg_shell = egg(r1, r2, h, opt(1), 2)
            rev1 = revolution(egg_shell, .2_wp)
            shell = onion(rev1, .2_wp)

            !albumen
            opt(2) = mono(1._wp, 0._wp, 0.0_wp, 1.37_wp)
            egg_albumen = egg(r1-.2_wp, r2, h, opt(2), 3)
            albumen = revolution(egg_albumen, .2_wp)

            !yolk
            opt(3) = mono(10._wp, 1._wp, 0.9_wp, 1.37_wp)
            yolk = sphere(1.5_wp, opt(3), 1)

            !bounding box
            opt(4) = mono(0._wp, 0._wp, 0.0_wp, 1._wp)
            bbox = box(vector(20.001_wp, 20.001_wp, 20.001_wp), opt(4), 4)
            
            allocate(array(4))
            
            array(1) = yolk
            array(2) = albumen
            array(3) = shell
            array(4) = bbox

        end function setup_egg

        function setup_sphere_scene(dict) result(array)

            use sdfNew,       only : sdf, sphere, box
            use sdfHelpers,   only : translate
            use random,       only : ranu
            use vector_class, only : vector
            use mat_class,    only : invert
            use opticalProperties, only : opticalProp_t, mono

            type(toml_table), intent(inout) :: dict
            type(sdf), allocatable :: array(:)
            
            integer :: num_spheres, i
            real(kind=wp) :: t(4,4), mus, mua, hgg, n, radius
            type(vector) :: pos
            type(opticalProp_t) :: opt(2)

            call get_value(dict, "num_spheres", num_spheres)
            allocate(array(num_spheres+1))

            mus = 1e-17_wp
            mua = 1e-17_wp
            hgg = 0.0_wp
            n   = 1.0_wp

            opt(2) = mono(mus, mua, hgg, n)

            array(num_spheres+1) = box(vector(2._wp, 2._wp, 2._wp), opt(2), num_spheres+1)
            
            mus = 0.0_wp!ranu(1._wp, 50._wp)
            mua = 0.0_wp!ranu(0.01_wp, 1._wp)
            hgg = 0.9_wp
            n = 1.37_wp
            opt(1) = mono(mus, mua, hgg, n)
            do i = 1, num_spheres
                radius = ranu(0.001_wp, 0.25_wp)
                pos = vector(ranu(-1._wp+radius, 1._wp-radius), ranu(-1._wp+radius, 1._wp-radius),&
                             ranu(-1._wp+radius, 1._wp-radius))
                t = invert(translate(pos))
    
                array(i) = sphere(radius, opt(1), i, transform=t)
            end do

        end function setup_sphere_scene


        function setup_logo() result(array)
        ! setup uni crest geometry
        !
        !
            use vector_class
            use sdfNew, only : sdf, box, segment, extrude
            use opticalProperties

            type(sdf), allocatable :: array(:)
            type(segment), allocatable, save :: seg(:)

            type(opticalProp_t) :: opt(2)

            type(vector)  :: a, b
            real(kind=wp) :: hgg, mus, mua, n
            integer       :: layer
            logical       :: fexists

            allocate(array(726), seg(725))

            mus = 10._wp
            mua = .1_wp
            hgg = 0.9_wp
            n = 1.5_wp
            layer = 1

            opt(1) = mono(0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp)
            opt(2) = mono(mus, mua, hgg, n)

            inquire(file="res/svg.f90", exist=fexists)
            if(.not. fexists)error stop "need to generate svg.f90 and place in res/"
            error stop "need to uncomment inlcude line!"
            ! include "../res/svg.f90"
            array(726) = box(vector(10._wp, 10._wp, 2.001_wp), opt(1), 2) 

        end function setup_logo


        function setup_sphere() result(array)
        ! setup the sphere test case from tran and jacques paper.
        !
        !
            use sdfNew,       only : sdf, box, sphere
            use sdfHelpers,   only : translate
            use vector_class, only : vector
            use mat_class,    only : invert
            use opticalProperties, only : mono, opticalProp_t
            
            type(sdf), allocatable :: array(:)
            type(opticalProp_t) :: opt(3)

            real(kind=wp) :: mus, mua, n, hgg, t(4, 4)
            type(vector)  :: a
            
            allocate(array(3))
            mus = 0._wp; mua = 1.e-17_wp; hgg = 0._wp; n = 1._wp;
            opt(1) = mono(mus, mua, hgg, n)
            array(2) = box(vector(2._wp, 2._wp, 2._wp), opt(1), 2)
            opt(2) = mono(mus, 10000000._wp, hgg, n)
            array(3) = box(vector(2.01_wp, 2.01_wp, 2.01_wp), opt(2), 3)

            mus = 0._wp; mua = 1.e-17_wp; hgg = 0._wp; n = 1.33_wp;
            opt(3) = mono(mus, mua, hgg, n)
            a = vector(.0_wp, 0._wp, 0._wp)
            t = invert(translate(a))
            array(1) = sphere(0.5_wp, opt(3), 1, transform=t)

        end function setup_sphere


        function setup_exp(dict) result(array)
        ! Setup experimental geometry from Georgies paper. i.e a glass bottle with contents
        !
        !
            use sdfNew,       only : sdf, box, cylinder!, subtraction
            use sdfHelpers,   only : rotate_y, translate
            use utils,        only : deg2rad
            use vector_class, only : vector
            use mat_class,    only : invert
            use opticalProperties, only : mono, opticalProp_t

            type(toml_table), intent(inout)  :: dict

            type(sdf), allocatable :: array(:)
            type(opticalProp_t) :: opt(3)

            type(vector)  :: a, b
            real(kind=wp) :: n, t(4, 4), optprop(5)

            call get_value(dict, "musb", optprop(1))
            call get_value(dict, "muab", optprop(2))
            call get_value(dict, "musc", optprop(3))
            call get_value(dict, "muac", optprop(4))
            call get_value(dict, "hgg", optprop(5))
            n = 1._wp

            opt(1) = mono(optprop(1), optprop(2), optprop(5), 1.5_wp)
            opt(2) = mono(optprop(3), optprop(4), optprop(5), 1.3_wp)

            a = vector(-10._wp, 0._wp, 0._wp)
            b = vector(10._wp, 0._wp, 0._wp)
            !bottle
            array(2) = cylinder(a, b, 1.75_wp, opt(1), 2)
            ! contents
            array(1) = cylinder(a, b, 1.55_wp, opt(2), 1)

            ! t = invert(translate(vector(0._wp, 0._wp, -5._wp+1.75_wp)))
            ! slab = box(vector(10._wp, 10._wp, 10._wp), optprop(3), optprop(4), optprop(5), 1.3_wp, 1, transform=t)
            opt(3) = mono(0.0_wp, 0.0_wp, 0.0_wp, n)
            array(3) = box(vector(4._wp, 4._wp, 4._wp), opt(3), 2)

        end function setup_exp

        function setup_scat_test(dict) result(array)

            use sdfNew, only : sdf, sphere, box

            use vector_class
            use opticalProperties

            type(toml_table), intent(inout) :: dict
            type(sdf), allocatable :: array(:)

            type(opticalProp_t) :: opt(2)
            real(kind=wp) :: mus, mua, hgg, n, tau

            call get_value(dict, "tau", tau)
            allocate(array(2))
            n = 1._wp
            hgg = 0.0_wp
            mua = 0.00_wp
            mus = tau

            opt(1) = mono(mus, mua, hgg, n)
            array(1) = sphere(1._wp, opt(1), 1)

            opt(2) = mono(0.0_wp, mua, hgg, n)
            array(2) = box(vector(2._wp, 2._wp, 2._wp), opt(2), 2)

        end function setup_scat_test

        function setup_scat_test2(dict) result(array)

            use sdfNew,    only : sdf, box

            use vector_class
            use opticalProperties

            type(toml_table), intent(inout) :: dict
            type(sdf), allocatable :: array(:)

            type(opticalProp_t) :: opt
            real(kind=wp) :: mus, mua, hgg, n, tau

            allocate(array(1))
            call get_value(dict, "tau", tau)
            call get_value(dict, "hgg", hgg)

            n = 1._wp
            hgg = hgg
            mua = 1e-17_wp
            mus = tau

            opt = mono(mus, mua, hgg, n)
            array(1) = box(vector(200._wp, 200._wp, 200._wp), opt, 2)

        end function setup_scat_test2

        function setup_omg_sdf() result(array)
            
            use sdfNew, only : sdf, cylinder, torus, box
            use sdfHelpers, only : translate, rotate_y
            use vector_class, only : vector
            use mat_class,    only : invert
            use opticalProperties

            type(sdf), allocatable :: array(:)
            
            type(opticalProp_t) :: opt(2)
            type(vector)  :: a, b
            real(kind=wp) :: t(4, 4), mus, mua, hgg, n
            integer       :: layer

            allocate(array(11))

            mus = 10._wp
            mua = 0.16_wp
            hgg = 0.0_wp
            n = 2.65_wp
            layer = 1

            opt(1) = mono(mus, mua, hgg, n)

            opt(2) = mono(0._wp, 0._wp, 0._wp, 1.0_wp)
            ! x
            ! |
            ! |
            ! |
            ! |
            ! |_____z

            !O letter
            a = vector(0._wp, 0._wp, -0.7_wp)
            t = invert(translate(a))
            array(1) = torus(.2_wp, 0.05_wp, opt(1), layer, transform=t)

            !M letter
            a = vector(-.25_wp, 0._wp, -.25_wp)
            b = vector(-.25_wp, 0._wp, .25_wp)
            t = invert(rotate_y(90._wp))
            array(2) = cylinder(a, b, .05_wp, opt(1), layer, transform=t)
            
            a = vector(-.25_wp, 0._wp, -.25_wp)
            b = vector(.25_wp, 0._wp, .0_wp)
            array(3) = cylinder(a, b, .05_wp, opt(1), layer)
            
            a = vector(.25_wp, 0._wp, .0_wp)
            b = vector(-.25_wp, 0._wp, .25_wp)
            array(4) = cylinder(a, b, .05_wp, opt(1), layer)

            a = vector(-.25_wp, 0._wp, .25_wp)
            b = vector(.25_wp, 0._wp, .25_wp)
            array(5) = cylinder(a, b, .05_wp, opt(1), layer)

            !G letter
            a = vector(-.25_wp, 0._wp, .5_wp)
            b = vector(.25_wp, 0._wp, .5_wp)
            array(6) = cylinder(a, b, .05_wp, opt(1), layer)

            a = vector(-.25_wp, 0._wp, .5_wp)
            b = vector(-.25_wp, 0._wp, .75_wp)
            array(7) = cylinder(a, b, .05_wp, opt(1), layer)

            a = vector(.25_wp, 0._wp, .5_wp)
            b = vector(.25_wp, 0._wp, .75_wp)
            array(8) = cylinder(a, b, .05_wp, opt(1), layer)

            a = vector(.25_wp, 0._wp, .75_wp)
            b = vector(0._wp, 0._wp, .75_wp)
            array(9) = cylinder(a, b, .05_wp, opt(1), layer)

            a = vector(0._wp, 0._wp, .625_wp)
            b = vector(0._wp, 0._wp, .75_wp)
            array(10) = cylinder(a, b, .05_wp, opt(1), layer)

            !bbox
            array(11) = box(vector(2._wp, 2._wp, 2._wp), opt(2), 2)

            ! allocate(array(2))
            ! allocate(array(1)%p, source=omg_sdf)
            ! allocate(array(2)%p, source=boxy)

            ! do j = 1, size(m)
            !     allocate(cnta(j)%p, source=m(j))
            ! end do

            ! do j = 1, size(g)
            !     allocate(cnta(j+size(m))%p, source=g(j))
            ! end do

            ! cnta(1)%p => m(1)
            ! cnta(2)%p => m(2)
            ! cnta(3)%p => m(3)
            ! cnta(4)%p => m(4)
            ! cnta(5)%p => sph
            ! cnta(6)%p => g(1)
            ! cnta(7)%p => g(2)
            ! cnta(8)%p => g(3)
            ! cnta(9)%p => g(4)
            ! cnta(10)%p => g(5)

            ! omg_sdf = model_init(cnta, union, 0.05_wp)

            ! array(1)%p => omg_sdf ! model
            ! array(2)%p => boxy    ! bbox

        end function setup_omg_sdf


        function get_vessels() result(array)

            ! use sdfs,         only : container, model, capsule, model_init, union, box, union
            use sdfNew, only : sdf, capsule, box
            use vector_class, only : vector

            use opticalProperties

            type(sdf), allocatable :: array(:)

            real(kind=wp), allocatable :: nodes(:, :), radii(:)
            integer, allocatable :: edges(:, :)
            integer :: io, edge_cnt, tmp1, tmp2, u, node_cnt, i
            real(kind=wp) :: x, y, z, radius, res, maxx, maxy, maxz
            real(kind=wp) :: musv, muav, gv, nv
            real(kind=wp) :: musd, muad, gd, nd
            type(vector) :: a, b

            type(opticalProp_t) :: opt(2)

            !MCmatlab: an open-source, user-friendly, MATLAB-integrated three-dimensional Monte Carlo light transport solver with heat diffusion and tissue damage
            muav = 231._wp
            musv = 94._wp
            gv = 0.9_wp
            nv = 1.37_wp

            muad = 0.458_wp
            musd = 357._wp
            gd = 0.9_wp
            nd = 1.37_wp

            opt(1) = mono(musv, muav, gv, nv)
            opt(2) = mono(musd, muad, gd, nd)

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

            allocate(array(edge_cnt+1))

            do i = 1, edge_cnt
                a = vector(nodes(edges(i, 1), 1), nodes(edges(i, 1), 2), nodes(edges(i, 1), 3))
                b = vector(nodes(edges(i, 2), 1), nodes(edges(i, 2), 2), nodes(edges(i, 2), 3))
                radius = radii(edges(i, 1)) * res
                array(i) = capsule(a, b, radius, opt(1), 1)
            end do

            array(i) = box(vector(.32_wp, .18_wp, .26_wp), opt(2), 2)

        end function get_vessels

        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : homedir, fileplace, resdir

            character(len=256) :: cwd
            logical :: dataExists, jmeanExists, depositExists, detectorsExists, phasorExists

            !get current working directory
            call get_environment_variable('PWD', cwd)
  
            ! get 'home' dir from cwd
            homedir = trim(cwd)
            ! get data dir
            fileplace = trim(homedir)//'/data/'
            !check if data directory and subdirectories exists. if not create it
#ifdef __GFORTRAN__
            inquire(file=trim(fileplace)//"/.", exist=dataExists)
            inquire(file=trim(fileplace)//"/jmean/.", exist=jmeanExists)
            inquire(file=trim(fileplace)//"/deposit/.", exist=depositExists)
            inquire(file=trim(fileplace)//"/detectors/.", exist=detectorsExists)
            inquire(file=trim(fileplace)//"/phasor/.", exist=phasorExists)
#elif __INTEL_COMPILER
            inquire(directory=trim(fileplace), exist=dataExists)
            inquire(directory=trim(fileplace)//"/jmean", exist=jmeanExists)
            inquire(directory=trim(fileplace)//"/deposit", exist=depositExists)
            inquire(directory=trim(fileplace)//"/detectors", exist=detectorsExists)
            inquire(directory=trim(fileplace)//"/phasor", exist=phasorExists)
#else 
    error stop "Compiler not supported!"
#endif
            if(.not. dataExists)then
                call create_directory("", dataExists, "", .false.)
                call create_directory("jmean/", jmeanExists, "data/", .false.)
                call create_directory("deposit/", depositExists, "data/", .false.)
                call create_directory("detectors/", detectorsExists, "data/", .false.)
                call create_directory("phasor/", phasorExists, "data/", .false.)
            else
                call create_directory("jmean/", jmeanExists, "data/", .true.)
                call create_directory("deposit/", depositExists, "data/", .true.)
                call create_directory("detectors/", detectorsExists, "data/", .true.)
                call create_directory("phasor/", phasorExists, "data/", .true.)
            end if

            ! get res dir
            resdir = trim(homedir)//'/res/'

        end subroutine directory


        subroutine create_directory(name, flag, appendname, newline)

            use constants, only : fileplace

            character(*),      intent(in) :: name, appendname
            logical,           intent(in) :: flag
            logical, optional, intent(in) :: newline

            character(len=:), allocatable :: mkdirCMD

            if(.not. flag)then
                mkdirCMD = "mkdir -p "//trim(fileplace)//name
                call execute_command_line(mkdirCMD)
                ! output correct message for base data dir
                if(len(name) == 0)then
                    mkdirCMD = "Created "//appendname//"data/"                    
                else
                    mkdirCMD = "Created "//appendname//name
                end if
                if(newline)mkdirCMD = mkdirCMD//new_line("a")
                print*,mkdirCMD
            end if

        end subroutine create_directory

        subroutine zarray

            use iarray

            !sets all arrays to zer

            phasor = 0._wp
            phasorGLOBAL = 0._wp
            jmean = 0._wp
            jmeanGLOBAL = 0._wp
            absorb = 0.0_wp
            absorbGLOBAL = 0.0_wp

        end subroutine zarray


        subroutine alloc_array(nxg, nyg, nzg)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iarray

            integer, intent(IN) :: nxg, nyg, nzg

            allocate(phasor(nxg, nyg, nzg), phasorGLOBAL(nxg, nyg, nzg))
            allocate(jmean(nxg, nyg, nzg), jmeanGLOBAL(nxg, nyg, nzg))
            allocate(absorb(nxg, nyg, nzg), absorbGLOBAL(nxg, nyg, nzg))

        end subroutine alloc_array

        subroutine dealloc_array()

            use iarray

            deallocate(jmean)
            deallocate(jmeanGLOBAL)
            deallocate(absorb)
            deallocate(absorbGLOBAL)
            deallocate(phasor)
            deallocate(phasorGLOBAL)
        end subroutine dealloc_array
end module subs