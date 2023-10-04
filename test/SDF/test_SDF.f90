module testsSDFMod

    use constants, only : wp
    use opticalProperties, only : mono, opticalProp_t
    use sdfs
    use sdfHelpers
    use testdrive, only : new_unittest, unittest_type, error_type, check, testsuite_type, new_testsuite, context_t
    use vector_class

    implicit none

    contains

    subroutine SDF_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("SDF Shapes", collect_suite1, context),&
                    !   new_testsuite("SDF Modifiers", collect_suite2, context),&
                      new_testsuite("SDF Helpers", collect_suite3, context)&
                    ]

    end subroutine SDF_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("sphere_test", test_sphere), &
                new_unittest("box_test", test_box), &
                new_unittest("cylinder_test", test_cylinder), &
                new_unittest("torus_test", test_torus), &
                new_unittest("triprism_test", test_triprisim), &
                new_unittest("cone_test", test_cone), &
                new_unittest("capsule_test", test_capsule), &
                new_unittest("plane_test", test_plane), &
                new_unittest("segment_test", test_segment), &
                new_unittest("egg_test", test_egg)&
                ]
    end subroutine collect_suite1

    ! subroutine collect_suite2(testsuite)

    !     type(unittest_type), allocatable, intent(out) :: testsuite(:)

    !     testsuite = [ &
    !             new_unittest("union_test", test_union) &
                ! new_unittest("intersection_test", test_intersection), &
                ! new_unittest("subtraction_test", test_subtraction), &
                ! new_unittest("displacement_test", test_displacement), &
                ! new_unittest("bend_test", test_bend), &
                ! new_unittest("twist_test", test_twist), &
                ! new_unittest("elongate_test", test_elongate), &
                ! new_unittest("repeat_test", test_repeat), &
                ! new_unittest("extrude_test", test_extrude), &
                ! new_unittest("revolution_test", test_revolution), &
                ! new_unittest("onion_test", test_onion), &
                ! new_unittest("translate_test", test_translate) &
                ! ]

    ! end subroutine collect_suite2


    subroutine collect_suite3(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("rotate x test", test_rotate_x), &
                new_unittest("rotate y test", test_rotate_y), &
                new_unittest("rotate z test", test_rotate_z), &
                new_unittest("translate test", test_translate), &
                new_unittest("rotate matrix test", test_rotmat), &
                new_unittest("identity test", test_identity), &
                new_unittest("skewsymmetric test", test_skewsymm), &
                new_unittest("rotation align test", test_rotationalign) &
                ]

    end subroutine collect_suite3

    subroutine test_rotationalign(error)

        type(error_type), allocatable, intent(out) :: error
        real(kind=wp) :: id(4, 4)
        type(vector) :: a, b, c

        a = vector(0., 0., 1.)
        b = vector(1., 0., 0.)

        id = rotationAlign(a, b)
        c = a .dot. id
        call check(error, c%x, b%x)
        if(allocated(error))return

        call check(error, c%y, b%y)
        if(allocated(error))return

        call check(error, c%z, b%z)
        if(allocated(error))return


        a = vector(1., 2., 1.)
        b = vector(1., 4., 5.)

        a = a%magnitude()
        b = b%magnitude()

        id = rotationAlign(a, b)
        c = a .dot. id
        call check(error, c%x, b%x)
        if(allocated(error))return

        call check(error, c%y, b%y)
        if(allocated(error))return

        call check(error, c%z, b%z)
        if(allocated(error))return

    end subroutine test_rotationalign

    subroutine test_rotmat(error)

        use utils, only : deg2rad

        type(error_type), allocatable, intent(out) :: error
        real(kind=wp) :: id(4, 4), angle
        type(vector) :: axis

        axis = vector(0._wp, 0._wp, 1._wp)
        angle = 45._wp
        id = rotmat(axis, angle)
        call check(error, all(id==rotate_z(angle)), .true.)
        if(allocated(error))return

        axis = vector(0._wp, 1._wp, 0._wp)
        angle = 45._wp
        id = rotmat(axis, angle)
        call check(error, all(id==rotate_y(angle)), .true.)
        if(allocated(error))return

        axis = vector(1._wp, 0._wp, 0._wp)
        angle = 45._wp
        id = rotmat(axis, angle)
        call check(error, all(id==rotate_x(angle)), .true.)
        if(allocated(error))return

    end subroutine test_rotmat

    subroutine test_rotate_x(error)
        
        use utils, only : deg2rad

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: id(4, 4), angle

        angle = 45._wp
        id = rotate_x(angle)
        angle = deg2rad(angle)

        call check(error, id(1, 1), 1.0_wp)
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), 0._wp)
        if(allocated(error))return

        call check(error, id(1, 2), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 2), cos(angle))
        if(allocated(error))return
        call check(error, id(3, 2), -sin(angle))
        if(allocated(error))return
        
        call check(error, id(1, 3), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 3), sin(angle))
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return
        angle = 90._wp
        id = rotate_x(angle)
        angle = deg2rad(angle)
        call check(error, id(1, 1), 1.0_wp)
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), 0._wp)
        if(allocated(error))return

        call check(error, id(1, 2), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 2), cos(angle))
        if(allocated(error))return
        call check(error, id(3, 2), -sin(angle))
        if(allocated(error))return
        
        call check(error, id(1, 3), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 3), sin(angle))
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return


        angle = 0._wp
        id = rotate_x(angle)
        angle = deg2rad(angle)
        call check(error, id(1, 1), 1.0_wp)
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), 0._wp)
        if(allocated(error))return

        call check(error, id(1, 2), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 2), cos(angle))
        if(allocated(error))return
        call check(error, id(3, 2), -sin(angle))
        if(allocated(error))return
        
        call check(error, id(1, 3), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 3), sin(angle))
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

        angle = 64.45_wp
        id = rotate_x(angle)
        angle = deg2rad(angle)
        call check(error, id(1, 1), 1.0_wp)
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), 0._wp)
        if(allocated(error))return

        call check(error, id(1, 2), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 2), cos(angle))
        if(allocated(error))return
        call check(error, id(3, 2), -sin(angle))
        if(allocated(error))return
        
        call check(error, id(1, 3), 0._wp)
        if(allocated(error))return
        call check(error, id(2, 3), sin(angle))
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

    end subroutine test_rotate_x

    subroutine test_rotate_y(error)
        
        use utils, only: deg2rad

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: id(4, 4), angle

        angle = 45._wp
        id = rotate_y(angle)
        angle = deg2rad(angle)

        call check(error, id(1, 1), cos(angle))
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), sin(angle))
        if(allocated(error))return

        call check(error, id(1, 2), 0.0_wp)
        if(allocated(error))return
        call check(error, id(2, 2), 1.0_wp)
        if(allocated(error))return
        call check(error, id(3, 2), 0.0_wp)
        if(allocated(error))return
        
        call check(error, id(1, 3), -sin(angle))
        if(allocated(error))return
        call check(error, id(2, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

        angle = 90._wp
        id = rotate_y(angle)
        angle = deg2rad(angle)

        call check(error, id(1, 1), cos(angle))
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), sin(angle))
        if(allocated(error))return

        call check(error, id(1, 2), 0.0_wp)
        if(allocated(error))return
        call check(error, id(2, 2), 1.0_wp)
        if(allocated(error))return
        call check(error, id(3, 2), 0.0_wp)
        if(allocated(error))return
        
        call check(error, id(1, 3), -sin(angle))
        if(allocated(error))return
        call check(error, id(2, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

        angle = 0._wp
        id = rotate_y(angle)
        angle = deg2rad(angle)
        call check(error, id(1, 1), cos(angle))
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), sin(angle))
        if(allocated(error))return

        call check(error, id(1, 2), 0.0_wp)
        if(allocated(error))return
        call check(error, id(2, 2), 1.0_wp)
        if(allocated(error))return
        call check(error, id(3, 2), 0.0_wp)
        if(allocated(error))return
        
        call check(error, id(1, 3), -sin(angle))
        if(allocated(error))return
        call check(error, id(2, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

        angle = 64.45_wp
        id = rotate_y(angle)
        angle = deg2rad(angle)
        call check(error, id(1, 1), cos(angle))
        if(allocated(error))return
        call check(error, id(2, 1), 0._wp)
        if(allocated(error))return
        call check(error, id(3, 1), sin(angle))
        if(allocated(error))return

        call check(error, id(1, 2), 0.0_wp)
        if(allocated(error))return
        call check(error, id(2, 2), 1.0_wp)
        if(allocated(error))return
        call check(error, id(3, 2), 0.0_wp)
        if(allocated(error))return
        
        call check(error, id(1, 3), -sin(angle))
        if(allocated(error))return
        call check(error, id(2, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(3, 3), cos(angle))
        if(allocated(error))return

    end subroutine test_rotate_y

    subroutine test_rotate_z(error)
        
        use utils, only: deg2rad

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: id(4, 4), angle

        angle = 45._wp
        id = rotate_z(angle)
        angle = deg2rad(angle)

        call check(error, id(1, 1), cos(angle))
        if(allocated(error))return
        call check(error, id(2, 1), -sin(angle))
        if(allocated(error))return
        call check(error, id(3, 1), 0.0_wp)
        if(allocated(error))return

        call check(error, id(1, 2), sin(angle))
        if(allocated(error))return
        call check(error, id(2, 2), cos(angle))
        if(allocated(error))return
        call check(error, id(3, 2), 0.0_wp)
        if(allocated(error))return
        
        call check(error, id(1, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(2, 3), 0.0_wp)
        if(allocated(error))return
        call check(error, id(3, 3), 1.0_wp)
        if(allocated(error))return

    end subroutine test_rotate_z


    subroutine test_identity(error)
        
        use stdlib_linalg, only : eye

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: id(4, 4)
        integer :: i(4, 4)

        i = eye(4, 4)
        id = identity()
        call check(error, all(id == real(i, kind=wp)), .true.)
        if(allocated(error))return

    end subroutine test_identity

    subroutine test_translate(error)
        
        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: t(4, 4)
        type(vector)  :: pos
        
        pos = vector(1., 2., 3.)
        t = translate(pos)

        call check(error, t(4, 1), pos%x)
        if(allocated(error))return
        call check(error, t(4, 2), pos%y)
        if(allocated(error))return
        call check(error, t(4, 3), pos%z)
        if(allocated(error))return

    end subroutine test_translate

    subroutine test_skewsymm(error)
        
        use stdlib_linalg, only : is_skew_symmetric

        type(error_type), allocatable, intent(out) :: error

        real(kind=wp) :: id(4, 4)
        type(vector) :: pos

        pos = vector(1._wp, 2._wp, 3._wp)
        id = skewSymm(pos)
        call check(error, is_skew_symmetric(id), .true.)
        if(allocated(error))return

        pos = vector(-1._wp, -2._wp, -3._wp)
        id = skewSymm(pos)
        call check(error, is_skew_symmetric(id), .true.)
        if(allocated(error))return

    end subroutine test_skewsymm

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

    subroutine test_triprisim(error)

        type(error_type), allocatable, intent(out) :: error

        type(triprism)  :: tri
        type(opticalProp_t) :: opt
        type(vector)  :: pos
        real(kind=wp) :: h1, h2

        ! height
        h1 = 1.0_wp
        ! length
        h2 = 5.0_wp
        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        tri = triprism(h1, h2, opt, 1)

        pos = vector(0.0_wp, 0._wp, 5._wp)
        call check(error, tri%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, 1._wp, 0._wp)
        call check(error, tri%evaluate(pos), 0.0_wp)
        if(allocated(error))return

    end subroutine test_triprisim

    subroutine test_capsule(error)

        type(error_type), allocatable, intent(out) :: error

        type(capsule)  :: cap
        type(opticalProp_t) :: opt
        type(vector)  :: pos, a, b
        real(kind=wp) :: r

        ! start
        a = vector(-1.0_wp, 0.0_wp, 0.0_wp)
        ! end
        b = vector(1.0_wp, 0.0_wp, 0.0_wp)
        ! radius
        r = 1.0_wp
        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        cap = capsule(a, b, r, opt, 1)

        pos = vector(0.0_wp, 0._wp, 0._wp)
        call check(error, cap%evaluate(pos), -1.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, 1._wp, 0._wp)
        call check(error, cap%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(2.0_wp, 0._wp, 0._wp)
        call check(error, cap%evaluate(pos), 0.0_wp)
        if(allocated(error))return

    end subroutine test_capsule

    subroutine test_plane(error)

        type(error_type), allocatable, intent(out) :: error

        type(plane)  :: plan
        type(opticalProp_t) :: opt
        type(vector)  :: pos, a

        ! normal to plane
        a = vector(0.0_wp, 0.0_wp, 1.0_wp)
        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        plan = plane(a, opt, 1)

        pos = vector(0.0_wp, 0._wp, 0._wp)
        call check(error, plan%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, 1._wp, 0._wp)
        call check(error, plan%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(2.0_wp, 0._wp, 0._wp)
        call check(error, plan%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, 0._wp, -1._wp)
        call check(error, plan%evaluate(pos), -1.0_wp)
        if(allocated(error))return

        pos = vector(0.0_wp, 0._wp, 1._wp)
        call check(error, plan%evaluate(pos), 1.0_wp)
        if(allocated(error))return


    end subroutine test_plane

    subroutine test_cone(error)

        type(error_type), allocatable, intent(out) :: error

        type(cone)  :: con
        type(opticalProp_t) :: opt
        type(vector)  :: pos, a, b
        real(kind=wp) :: ra, rb

        ! radius of base
        ra = 5.0_wp
        ! radius of tip
        rb = 0.0_wp
        opt = mono(0._wp, 0._wp, 0._wp, 0._wp)
        a = vector(0.0, 0.0, 0.0)
        b = vector(0.0, 0.0, 1.0)
        con = cone(a, b, ra, rb, opt, 1)

        pos = vector(0.0_wp, 0._wp, 1._wp)
        call check(error, con%evaluate(pos), 0.0_wp)
        if(allocated(error))return

        pos = vector(1.0_wp, 1._wp, 0._wp)
        call check(error, con%evaluate(pos), 0.0_wp)
        if(allocated(error))return

    end subroutine test_cone

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