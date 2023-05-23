module testsSDFMod

    use constants, only : wp
    use opticalProperties, only : mono, opticalProp_t
    use sdfs
    use testdrive, only : new_unittest, unittest_type, error_type, check, testsuite_type, new_testsuite, context_t
    use vector_class

    implicit none

    contains

    subroutine SDF_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("SDF", collect_suite1, context)&
                     ]

    end subroutine SDF_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("sphere_test", test_sphere), &
                new_unittest("box_test", test_box), &
                new_unittest("cylinder_test", test_cylinder), &
                new_unittest("torus_test", test_torus), &
                ! new_unittest("triprism_test", test_triprism), &
                ! new_unittest("cone_test", test_cone), &
                ! new_unittest("capsule_test", test_capsule), &
                ! new_unittest("plane_test", test_plane), &
                new_unittest("segment_test", test_segment), &
                new_unittest("egg_test", test_egg)&
                ]
    end subroutine collect_suite1

    subroutine test_sphere(error)

        type(error_type), allocatable, intent(out) :: error

        type(sphere)  :: sph
        type(opticalProp_t) :: opt
        type(vector) :: pos

        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        sph = sphere(1.0_wp, opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, sph%evaluate(pos), -1._wp)
        if(allocated(error))return

        pos = vector(0._wp, 1._wp, 0._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, 1._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(1._wp, 0._wp, 0._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, -1._wp, 0._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, -1._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(-1._wp, 0._wp, 0._wp)
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(sqrt(1./3._wp), sqrt(1./3._wp), sqrt(1./3._wp))
        call check(error, sph%evaluate(pos), 0._wp)
        if(allocated(error))return

    end subroutine test_sphere

    subroutine test_box(error)

        type(error_type), allocatable, intent(out) :: error

        type(box)  :: bbox
        type(opticalProp_t) :: opt
        type(vector) :: pos

        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        bbox = box(vector(2.0_wp, 2.0_wp, 2.0_wp), opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, bbox%evaluate(pos), -1._wp)
        if(allocated(error))return

        pos = vector(0._wp, 1._wp, 0._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, 1._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(1._wp, 0._wp, 0._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, -1._wp, 0._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, -1._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(-1._wp, 0._wp, 0._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(1._wp, 1._wp, 1._wp)
        call check(error, bbox%evaluate(pos), 0._wp)
        if(allocated(error))return

    end subroutine test_box

    subroutine test_cylinder(error)

        type(error_type), allocatable, intent(out) :: error

        type(cylinder)  :: cyl
        type(opticalProp_t) :: opt
        type(vector) :: pos, a, b

        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        a = vector(0._wp, 0._wp, -1._wp)
        b = vector(0._wp, 0._wp, 1._wp)
        cyl = cylinder(a, b, 1._wp, opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, cyl%evaluate(pos), -1._wp)
        if(allocated(error))return

        pos = vector(0._wp, 1._wp, 0._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, 1._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(1._wp, 0._wp, 0._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, -1._wp, 0._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, -1._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(-1._wp, 0._wp, 0._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

        pos = vector(sqrt(1./2._wp), sqrt(1./2._wp), 0._wp)
        call check(error, cyl%evaluate(pos), 0._wp)
        if(allocated(error))return

    end subroutine test_cylinder

    subroutine test_torus(error)

        type(error_type), allocatable, intent(out) :: error

        type(torus)  :: tor
        type(opticalProp_t) :: opt
        type(vector) :: pos

        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        tor = torus(0.5_wp, 1.0_wp, opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, tor%evaluate(pos), -.5_wp)
        if(allocated(error))return

        pos = vector(1.5_wp, 0._wp, 0._wp)
        call check(error, tor%evaluate(pos), 0._wp)
        if(allocated(error))return

    end subroutine test_torus

    subroutine test_segment(error)

        type(error_type), allocatable, intent(out) :: error

        type(segment)  :: seg
        type(opticalProp_t) :: opt
        type(vector) :: pos, a, b

        a = vector(-1._wp, 0., 0._wp)
        b = vector(1._wp, 0., 0._wp)
        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        seg = segment(a, b, opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, seg%evaluate(pos), -0.1_wp)
        if(allocated(error))return

        pos = vector(-1._wp, 0._wp, 0._wp)
        call check(error, seg%evaluate(pos), -0.1_wp)
        if(allocated(error))return

        pos = vector(1._wp, 0._wp, 0._wp)
        call check(error, seg%evaluate(pos), -0.1_wp)
        if(allocated(error))return

        pos = vector(1._wp, 1.1_wp, 0._wp)
        call check(error, seg%evaluate(pos), 1._wp)
        if(allocated(error))return

        pos = vector(0._wp, 1.1_wp, 0._wp)
        call check(error, seg%evaluate(pos), 1._wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, 1.1_wp)
        call check(error, seg%evaluate(pos), 1._wp)
        if(allocated(error))return

    end subroutine test_segment

    subroutine test_egg(error)

        type(error_type), allocatable, intent(out) :: error

        type(egg)  :: eggy
        type(opticalProp_t) :: opt
        type(vector) :: pos
        real(kind=wp) :: r1, r2, h

        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        ! makes a Moss egg. https://www.shadertoy.com/view/WsjfRt
        ! R1 controls "fatness" of the egg. Actually controls the base circle radius.
        ! R2 contorls the pointiness of the egg. Actually controls radius of top circle.
        ! h controls the height of the egg. Actually controls y position of top circle.
        r1 = 2.5_wp
        r2 = 0.75_wp
        h = 1.5_wp
        eggy = egg(r1, r2, h, opt, 1)

        pos = vector(0._wp, 0._wp, 0._wp)
        call check(error, eggy%evaluate(pos), -r1)
        if(allocated(error))return

        pos = vector(r1, 0._wp, 0._wp)
        call check(error, eggy%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(r1, 0._wp, 0._wp)
        call check(error, eggy%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, r1+2._wp*r2, 0._wp)
        call check(error, eggy%evaluate(pos), 0.0_wp, thr=1e-5_wp)
        if(allocated(error))return

        pos = vector(r1, r1, 0._wp)
        call check(error, eggy%evaluate(pos), 0.630294_wp, thr=1e-5_wp)
        if(allocated(error))return

    end subroutine test_egg
end module testsSDFMod