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

            integer,                  intent(OUT) :: nphotons
            type(cart_grid),          intent(OUT) :: grid
            type(container), allocatable,         intent(OUT) :: sdfarray(:)
            character(*), intent(IN) :: choice
            real, optional, intent(IN) :: tau

            real    :: xmax, ymax, zmax, albedo(3), hgg(3), g2(3),kappa(3), n1, n2
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
                read(u,*) n1
                read(u,*) n2
            close(u)

            !allocate and set arrays to 0
            call alloc_array(nxg, nyg, nzg)
            call zarray

            ! Set up grid
            grid = cart_grid(nxg, nyg, nzg, xmax, ymax, zmax, n1, n2)
            call init_opt1(kappa, albedo, hgg, g2)

            if(choice == "omg")then
                sdfarray = setup_omg_sdf()
            elseif(choice == "scat_test")then
                if(.not. present(tau))error stop "No tau provided"
                sdfarray = setup_scat_test(tau)
            elseif(choice == "fresnel_test")then
                sdfarray = setup_fresnel_test()
            else
                error stop "no such routine"
            end if

        end subroutine setup_simulation

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
            do i = 1, size(array)
                allocate(array(i)%p)
            end do

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
            do i = 1, size(array)
                allocate(array(i)%p)
            end do

            array(1)%p => sph
            array(2)%p => bbox

        end function setup_scat_test

        function setup_omg_sdf() result(array)
            
            use sdfs,      only : container, cylinder, torus, model, box, smoothunion, rotate_y, model_init
            use constants, only : pi
            
            use vector_class

            implicit none

            type(container), allocatable :: array(:)

            type(container) :: cnta(10)
            type(cylinder), target, save :: m(4), g(5)!save needed as these
            type(torus),    target, save :: sph       !are deallocated on sub exit
            type(box),      target, save :: boxy      !should be safe as this sub only called once
            type(model),    target, save :: omg_sdf
            real    :: t(3, 3), mus, mua, hgg, n
            integer :: j, layer

            mus = 10.
            mua = 0.
            hgg = 0.0d0
            n = 2.65
            layer = 1

            !O letter
            sph = torus(.2, 0.05, mus, mua, hgg, n, layer, c=vector(0., 0., -0.75))
            
            !M letter
            m(1) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,-.25))
            
            t = rotate_y(-25.*pi/180.)
            m(2) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,-.125), transform=t)
            
            t = rotate_y(25.*pi/180.)
            m(3) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,.125), transform=t)

            m(4) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,.25))

            !G letter
            g(1) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,.6))

            t = rotate_y(90.*pi/180.)
            g(2) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(0.25,0.,.775),transform=t)

            t = rotate_y(90.*pi/180.)
            g(3) = cylinder(.25, .05, mus, mua, hgg, n, layer, c=vector(-0.25, 0.,.775),transform=t)

            g(4) = cylinder(.125, .05, mus, mua, hgg, n, layer, c=vector(-0.1,0.,1.))

            t = rotate_y(90.*pi/180.)
            g(5) = cylinder(.125, .05, mus, mua, hgg, n, layer, c=vector(0.,0.,.9),transform=t)

            !bbox
            boxy = box(2., 0.d0, 100., 0.d0, 1.0, 2)

            allocate(array(2))
            do j = 1, size(array)
                allocate(array(j)%p)
            end do

            do j = 1, size(cnta)
                allocate(cnta(j)%p)
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