module testsPhotonMod

    use photonMod
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, new_testsuite, testsuite_type, context_t
    use constants, only : wp

    implicit none

    private
    public :: photon_suite

    contains

    subroutine photon_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test Sources", collect_suite1, context),&
                      new_testsuite("Photon Misc. routines", collect_suite2, context)&
                     ]

    end subroutine photon_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Uniform_source", Uniform_src), &
                new_unittest("Pencil_source", pencil_src), &
                new_unittest("Point_source", point_src), &
                new_unittest("Circular_source", circular_src), &
                new_unittest("SLM_source", SLM_src) &
                ]
    end subroutine collect_suite1

    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [new_unittest("Test init photon", photon_init_test) &
                    ]
    end subroutine collect_suite2
    
    subroutine photon_init_test(error)

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet

        packet = photon(10._wp)

        call check(error, packet%nxp, 10._wp)
        if(allocated(error))return
        call check(error, packet%sint, 10._wp)
        if(allocated(error))return
        call check(error, packet%wavelength, 10._wp)
        if(allocated(error))return
        call check(error, packet%tflag, .true.)
        if(allocated(error))return
    end subroutine photon_init_test

    subroutine Uniform_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use piecewiseMod
        use gridMod

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        type(toml_table) :: dict
        type(spectrum_t) :: spectrum
        integer :: i
        type(constant) :: const

        const = constant(500._wp)
        state%grid = init_grid(200, 200, 200, 7.5_wp, 7.5_wp, 7.5_wp)

        pos = vector(0.0_wp,0.0_wp,0.0_wp)
        dir = vector(1.0_wp,0.0_wp,0.0_wp)

        call set_photon(pos, dir)
        packet = photon("uniform")
        dict = toml_table()

        call set_value(dict, "pos1%x", -7.5_wp)
        call set_value(dict, "pos1%y", -1.0_wp)
        call set_value(dict, "pos1%z", -1.0_wp)

        call set_value(dict, "pos2%x", 0.0_wp)
        call set_value(dict, "pos2%y", 2.0_wp)
        call set_value(dict, "pos2%z", 0.0_wp)

        call set_value(dict, "pos3%x", 0.0_wp)
        call set_value(dict, "pos3%y", 0.0_wp)
        call set_value(dict, "pos3%z", 2.0_wp)

        allocate(spectrum%p, source=const)

        do i = 1, 10000
            call packet%emit(spectrum, dict)
            if(packet%pos%x /= -7.5_wp)then
                call test_failed(error, "Uniform X coordinate failed!", "Pos%x="//str(packet%pos%x)//" should be -7.5")
                return
            end if
            if(packet%pos%y > 1.0_wp .or. packet%pos%y < -1._wp)then
                call test_failed(error, "Uniform Y coordinate failed!", "Pos%x="//str(packet%pos%y)//" should be in range [-1,1]")
                return
            end if
            if(packet%pos%z > 1.0_wp.or. packet%pos%z < -1.0_wp)then
                call test_failed(error, "Uniform Z coordinate failed!", "Pos%x="//str(packet%pos%z)//" should be in range [-1,1]")
                return
            end if
        end do
    end subroutine Uniform_src

    subroutine Point_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod
        use piecewiseMod

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        integer :: i
        type(spectrum_t) :: spectrum
        type(constant) :: const

        const = constant(500._wp)
        allocate(spectrum%p, source=const)

        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(0.0_wp,0.5_wp,-0.25_wp)
        dir = vector(1.0_wp,0.0_wp,0.0_wp)

        call set_photon(pos, dir)
        packet = photon("point")

        do i = 1, 1000
            call packet%emit(spectrum)
            if(packet%pos%x /= 0.0_wp)then
                call test_failed(error, "Uniform X coordinate failed!", "Pos%x="//str(packet%pos%x)//" should be 0.0")
                return
            end if
            if(packet%pos%y /= 0.5_wp)then
                call test_failed(error, "Uniform Y coordinate failed!", "Pos%x="//str(packet%pos%y)//" should be in 0.5")
                return
            end if
            if(packet%pos%z /= -0.25_wp)then
                call test_failed(error, "Uniform Z coordinate failed!", "Pos%x="//str(packet%pos%z)//" should be in -0.25")
                return
            end if
        end do
    end subroutine Point_src

    subroutine circular_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod
        use piecewiseMod
        use random, only : init_rng

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        integer :: i
        real(kind=wp) :: r_pos
        type(toml_table) :: dict
        type(spectrum_t) :: spectrum
        type(constant) :: const

        call init_rng(spread(12345678, 1, 8), .false.)

        const = constant(500._wp)
        allocate(spectrum%p, source=const)

        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(1.0_wp, 1.0_wp, 1.0_wp)
        dir = vector(0.0_wp, 0.0_wp, -1.0_wp)

        call set_photon(pos, dir)
        dict = toml_table()
        call set_value(dict, "radius", 2.5_wp)
        packet = photon("circular")

        do i = 1, 10000
            call packet%emit(spectrum, dict)
            r_pos = sqrt((packet%pos%x+1)**2 + (packet%pos%y+1)**2+ (packet%pos%z+1)**2)
            if(r_pos > 2.5_wp .and. packet%pos%z /= -1._wp)then
                call test_failed(error, "Circle Source failed!", "radial position="//str(r_pos, 5)//" should be less than 2.5")
                return
            end if
        end do
    end subroutine circular_src

    subroutine Pencil_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod
        use piecewiseMod

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        integer :: i
        type(spectrum_t) :: spectrum
        type(constant), target :: const

        const = constant(500._wp)
        allocate(spectrum%p, source=const)
        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(0.0_wp,0.5_wp,-0.25_wp)
        dir = vector(1.0_wp,0.0_wp,0.0_wp)

        call set_photon(pos, dir)
        packet = photon("pencil")

        do i = 1, 1000
            call packet%emit(spectrum)
            if(packet%pos%x /= 0.0_wp)then
                call test_failed(error, "Uniform X coordinate failed!", "Pos%x="//str(packet%pos%x)//" should be 0.0")
                return
            end if
            if(packet%pos%y /= 0.5_wp)then
                call test_failed(error, "Uniform Y coordinate failed!", "Pos%x="//str(packet%pos%y)//" should be in 0.5")
                return
            end if
            if(packet%pos%z /= -0.25_wp)then
                call test_failed(error, "Uniform Z coordinate failed!", "Pos%x="//str(packet%pos%z)//" should be in -0.25")
                return
            end if
        end do
    end subroutine Pencil_src

    subroutine SLM_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod
        use piecewiseMod
        use random, only : init_rng
        use stb_image_mod
        use, intrinsic :: iso_c_binding

        type(error_type), allocatable, intent(out) :: error

        type(piecewise2D) :: img
        type(spectrum_t)  :: spectrum
        real(kind=wp)     :: cell_height, cell_width, sum_dif
        type(photon)      :: packet
        type(vector)      :: pos, dir
        integer           :: i, err, height, width, n_channels, idx, idy
        real(kind=wp),    allocatable :: array(:,:), out(:,:)
        character(len=:), allocatable :: sfile
        integer,          allocatable :: image(:,:,:)

        cell_height = 2./200.
        cell_width = 2./200.

        call init_rng(spread(12345678, 1, 8), .false.)
        
        ! Read in image
        sfile = "test/parse/test.png"
        err = stbi_info(trim(sfile)//c_null_char, width, height, n_channels)
        if(err == 0)then
            print'(2a,1x,a)', "Error reading file: ", trim(sfile),stbi_failure_reason()
            stop 1
        end if
        image = stbi_load(trim(sfile)//c_null_char, width, height, n_channels, 0)
        allocate(array(size(image, 1), size(image, 2)))
        array = image(:,:,1)
        ! set all values to 1 in image
        where(array > 0.)
            array = 1.
        end where

        img = piecewise2D(cell_width, cell_height, array)
        allocate(spectrum%p, source=img)
        allocate(out(size(image, 1), size(image, 2)))
        out = 0.
        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(0.0_wp, 0.0_wp, 1.0_wp)
        dir = vector(0.0_wp, 0.0_wp, -1.0_wp)
        call set_photon(pos, dir)

        packet = photon("slm")

        do i = 1, 1000000
            call packet%emit(spectrum)
            !why is it always +2?
            idx = nint((packet%pos%x + 1.) / (2./200.))+2
            idy = nint((packet%pos%y + 1.) / (2./200.))+2
            if(idx < 1 .or. idy < 1)cycle
            out(idx,idy) = out(idx,idy) + 1.
        end do

        out = out /maxval(out)
        sum_dif = sum(abs(array - out)) / (200**2)
        call check(error, sum_dif, 0.0_wp, thr=6e-2_wp)
        if(allocated(error))return
    end subroutine SLM_src
end module testsPhotonMod