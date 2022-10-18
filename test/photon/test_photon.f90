module testsPhotonMod

    use photonMod
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Uniform_source", Uniform_src), &
                new_unittest("Pencil_source", pencil_src), &
                new_unittest("Point_source", point_src) &
                ! new_unittest("Circular_source", circular_src) &
                ]
    end subroutine collect_suite1

    subroutine Uniform_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        type(toml_table) :: dict
        type(cart_grid) :: grid
        integer :: i

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

        do i = 1, 10000
            call packet%emit(dict)
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

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        type(toml_table) :: dict
        type(cart_grid) :: grid
        integer :: i

        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(0.0_wp,0.5_wp,-0.25_wp)
        dir = vector(1.0_wp,0.0_wp,0.0_wp)

        call set_photon(pos, dir)
        packet = photon("point")

        do i = 1, 1000
            call packet%emit()
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

    subroutine Pencil_src(error)

        use vector_class, only : vector
        use tomlf, only : set_value, toml_table
        use sim_state_mod, only : state
        use utils, only : str
        use gridMod

        type(error_type), allocatable, intent(out) :: error

        type(photon) :: packet
        type(vector) :: pos, dir
        type(toml_table) :: dict
        type(cart_grid) :: grid
        integer :: i

        state%grid = init_grid(200, 200, 200, 1._wp, 1._wp, 1._wp)

        pos = vector(0.0_wp,0.5_wp,-0.25_wp)
        dir = vector(1.0_wp,0.0_wp,0.0_wp)

        call set_photon(pos, dir)
        packet = photon("pencil")

        do i = 1, 1000
            call packet%emit()
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
end module testsPhotonMod
program test_photon

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testsPhotonMod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("Suite: Test bounds", collect_suite1) &
    !               new_testsuite("Suite: Vector .op. scalar", collect_suite2), &
    !               new_testsuite("Suite: Vector functions", collect_suite3), &
    !               new_testsuite("Suite: Vector .op. matrix", collect_suite4) &
                  ]

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat)
    end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_photon