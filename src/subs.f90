MODULE subs

implicit none

public  :: setup_simulation, print_time, get_time, Sellmeier
private :: directory, alloc_array, zarray


    contains

        subroutine setup_simulation(nphotons, grid, sdfarray, choice, tau)
        ! Read in parameters
            use constants, only : resdir
            use ch_opt
            use gridMod,   only : cart_grid
            use sdfs, only : container

            implicit none

            integer,                      intent(OUT) :: nphotons
            type(cart_grid),              intent(OUT) :: grid
            character(*),                 intent(IN)  :: choice
            type(container), allocatable, intent(OUT) :: sdfarray(:)
            real,            optional,    intent(IN)  :: tau

            real    :: xmax, ymax, zmax, albedo(3), hgg(3), g2(3),kappa(3)
            integer :: nxg, nyg, nzg, u

            !set directory paths
            call directory

            open(newunit=u,file=trim(resdir)//'input.params',status='old')
                read(u,*) nphotons
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


            select case(choice)
                case("omg")
                    sdfarray = setup_omg_sdf()
                case("scat_test")
                    if(.not. present(tau))error stop "No tau provided"
                    sdfarray = setup_scat_test(tau)
                case("fresnel_test")
                    sdfarray = setup_fresnel_test()
                case("skin")
                    sdfarray = setup_skin_model()
                case("interior")
                    sdfarray = interior_test()
                case("sphere")
                    sdfarray = setup_sphere()
                case("exp")
                    sdfarray = setup_exp()
                case("jacques")
                    sdfarray = setup_jacques()
                case default
                    error stop "no such routine"
            end select

        end subroutine setup_simulation


        function setup_jacques() result(array)

            use sdfs, only : container, box
            use vector_class

            implicit none

            type(container)         :: array(2)
            type(box), target, save :: medium, bbox
            real                 :: hgg

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

        function setup_sphere() result(array)

            use sdfs, only : sphere, box, container
            use vector_class

            implicit none

            type(container), allocatable :: array(:)
            type(sphere),   target, save :: sph
            type(box),      target, save :: bbox
 
            real :: mus, mua, n, hgg
            integer :: i

            mus = 0.; mua = 0.; hgg = 0.; n = 1.;
            bbox = box(2., mus, mua, hgg, n, 2)
            
            mus = 0.; mua = .5; hgg = .0; n = 1.5;
            sph = sphere(.5, mus, mua, hgg, n, 1)

            allocate(array(2))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox)

            array(1)%p => sph
            array(2)%p => bbox

        end function setup_sphere


        function interior_test() result(array)

            use sdfs, only : container, model, union, sphere, box, translate, model_init
            use vector_class

            implicit none


            type(container) :: array(1), cnta(3)
            type(sphere), target, save :: sph
            type(box),    target, save :: boxes(2), bbox
            type(model),  target, save :: m

            integer :: i
            real :: t(4, 4)


            bbox = box(2., 0., 0., 0., 1., 2)

            t = invert(translate(vector(.25, .0, 0.)))
            boxes(1) = box(vector(.5, 1., 1.), 0., 0., 0., 1., 1, transform=t)
            t = invert(translate(vector(0., .75, 0.)))
            boxes(2) = box(vector(1., .5, 1.), 0., 0., 0., 1., 1, transform=t)

            t = invert(translate(vector(.5, 0., 0.)))
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


        function setup_exp() result(array)
            
            use sdfs,  only : container, box, cylinder, rotate_y
            use utils, only : deg2rad
            use vector_class

            implicit none

            type(container), allocatable :: array(:)
            type(cylinder), target, save :: cyl(2)
            type(box),      target, save :: bbox

            real    :: mus, mua, hgg, n, t(4, 4)
            integer :: i

            n = 1.d0
            hgg = 0.d0
            mua = 1.d-17
            mus = 0.d0

            ! t = rotate_y(deg2rad(90.))
            ! cyl(1) = cylinder(10., 1.55, mus, mua, hgg, 1.3, 1)
            ! cyl(2) = cylinder(10., 1.75, mus, mua, hgg, 1.5, 2)
            error stop "need to fix cylinders"
            bbox  = box(2., mus, mua, hgg, n, 3)

            allocate(array(3))
            allocate(array(1)%p, source=cyl(1))
            allocate(array(2)%p, source=cyl(2))
            allocate(array(3)%p, source=bbox)

            array(1)%p => cyl(1)
            array(2)%p => cyl(2)
            array(3)%p => bbox

        end function setup_exp


        function setup_fresnel_test() result(array)

            use sdfs, only : container, box
            use vector_class

            implicit none

            type(container), allocatable :: array(:)

            type(box), target, save :: bbox, fibre

            integer :: i, layer
            real    :: mus, mua, hgg, n

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

        function setup_scat_test(tau) result(array)

            use sdfs, only : container, sphere, box
            use vector_class

            implicit none

            type(container), allocatable :: array(:)
            real, intent(IN) :: tau

            type(sphere), target, save :: sph
            type(box), target, save :: bbox

            integer :: i, layer
            real :: mus, mua, hgg, n

            layer = 1
            n = 1.d0
            hgg = 0.d0
            mua = 1d-17
            mus = tau

            sph = sphere(1., mus, mua, hgg, n, layer)
            bbox = box(2., 0.d0, mua, hgg, n, 2)

            allocate(array(2))
            allocate(array(1)%p, source=sph)
            allocate(array(2)%p, source=bbox)

            array(1)%p => sph
            array(2)%p => bbox

        end function setup_scat_test

        function setup_omg_sdf() result(array)
            
            use sdfs,      only : container, cylinder, torus, model, box, smoothunion, rotate_y, model_init, translate
            use constants, only : pi
            
            use vector_class

            implicit none

            type(container), allocatable :: array(:)

            type(container) :: cnta(10)
            type(cylinder), target, save :: m(4), g(5)!save needed as these
            type(torus),    target, save :: sph       !are deallocated on sub exit
            type(box),      target, save :: boxy      !should be safe as this sub only called once
            type(model),    target, save :: omg_sdf
            
            type(vector) :: a, b
            real         :: t(4, 4), mus, mua, hgg, n
            integer      :: j, layer

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
            t = invert(translate(vector(0., 0., -0.7)))
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
            boxy = box(2., 0.d0, 10., 0.d0, 1.0, 2)

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

            omg_sdf = model_init(cnta, smoothunion)

            array(1)%p => omg_sdf ! model
            array(2)%p => boxy    ! bbox

        end function setup_omg_sdf


        function setup_skin_model() result(array)

            use sdfs, only : container, model, cylinder, box, translate, model_init, smoothunion
            use utils, only : deg2rad
            use vector_class

            implicit none

            type(container), allocatable :: array(:), cnta(:)

            type(cylinder), target, save :: cyls(12)
            type(box),      target, save :: skin(3)
            type(model),    target, save :: vessels

            integer :: i
            real    :: mus, mua, hgg, n, t(4, 4)

            real :: mus_epi, mus_derm, mua_epi, mua_derm, n_epi, n_derm, hgg_epi, hgg_derm
            real :: mus_b, mua_b, hgg_b, n_b
            type(vector) :: a, b, c

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
            
            a = vector(0., -0.05, 0.)
            b = vector(0., -0.03, 0.)
            cyls(1) = cylinder(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0., -0.03-.005, 0.)
            b = vector(0.03, -0.02, 0.)
            cyls(2) = cylinder(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.03, -0.02-.005, 0.)
            b = vector(0., 0.01, 0.)
            cyls(3) = cylinder(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0., 0.01-.005, 0.)
            b = vector(0., 0.03, 0.)
            cyls(4) = cylinder(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0., 0.03-.005, 0.)
            b = vector(0.03, 0.05, 0.)
            cyls(5) = cylinder(a, b, .005, mus_b, mua_b, hgg_b, n_b, 1)

            !branch 1
            a = vector(-0.0025, -0.04, 0.)
            b = vector(-0.03, -0.03, 0.)
            cyls(6) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03, -0.03, 0.)
            b = vector(-0.03, -0.02, 0.)
            cyls(7) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(-0.03, -0.02, 0.)
            b = vector(0.0025, 0.02, 0.)
            cyls(8) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)


            !branch 2
            a = vector(0.015, -0.025, 0.)
            b = vector(0.048, -0.02, 0.)
            cyls(9) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.048, -0.02, 0.)
            b = vector(0.04, 0.01, 0.)
            cyls(10) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.04, 0.01, 0.)
            b = vector(0.03, 0.025, 0.)
            cyls(11) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)

            a = vector(0.03, 0.025, 0.)
            b = vector(-0.0025, 0.02, 0.)
            cyls(12) = cylinder(a, b, .001, mus_b, mua_b, hgg_b, n_b, 1)


            allocate(array(4))
            ! do i = 1, size(array)
            allocate(array(1)%p, source=vessels)
            allocate(array(2)%p, source=skin(1))
            allocate(array(3)%p, source=skin(2))
            allocate(array(4)%p, source=skin(3))

            ! end do

            allocate(cnta(size(cyls)))
            do i = 1, size(cnta)
                allocate(cnta(i)%p, source=cyls(i))
                cnta(i)%p => cyls(i)
            end do

            vessels = model_init(cnta, smoothunion)
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
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir = trim(homedir)//'res/'

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
                if(time >= 60. .and. id == 0)then
                   print*, floor((time)/60.) + mod(time, 60.)/100.
                else
                   print*, 'time taken ~',time,'s'
                end if
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