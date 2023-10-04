module testsgeometryMod

    use geometry
    use vector_class
    use testdrive, only : new_unittest, unittest_type, error_type, check, testsuite_type, new_testsuite, context_t
    use constants,    only : wp

    implicit none

    private
    public :: geometry_suite

    contains

    subroutine geometry_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Geometry intersections and routines", collect_suite1, context)&
                    ]

    end subroutine geometry_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("intersect sphere", test_sphere), &
                new_unittest("intersect plane", test_plane), &
                new_unittest("intersect cone", test_cone), &
                new_unittest("intersect cylinder", test_cylinder), &
                new_unittest("intersect ellipse", test_ellipse), &
                new_unittest("intersect circle", test_circle)&
                ]
    end subroutine collect_suite1

    subroutine test_sphere(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: pos, dir, centre, orig, new_pos
        real(kind=wp) :: t, radius

        centre = vector(0.0, 0.0, 0.0)
        radius = 10._wp

        orig = vector(-100., 0., 0.)
        dir = vector(1., 0., 0.)
        flag = intersectSphere(orig, dir, t, centre, radius)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, -10.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, 0.0_wp)
        if(allocated(error))return

        orig = vector(-100., 0., 0.)
        dir = vector(0., 0., 1.)
        flag = intersectSphere(orig, dir, t, centre, radius)

        call check(error, flag, .false.)
        if(allocated(error))return

    end subroutine test_sphere

    subroutine test_cylinder(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: pos, dir, centre, orig, new_pos
        real(kind=wp) :: t, radius

        centre = vector(0.0, 0.0, 0.0)
        radius = 10._wp

        orig = vector(-20.0, 0., 0.)
        dir = vector(1., 0., 0.)
        flag = intersectCylinder(orig, dir, t, centre, radius)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, -10.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, 0.0_wp)
        if(allocated(error))return

        orig = vector(-20.0, 0., 0.)
        dir = vector(-1., 0., 0.)
        flag = intersectCylinder(orig, dir, t, centre, radius)

        call check(error, flag, .false.)
        if(allocated(error))return
    end subroutine test_cylinder

    subroutine test_ellipse(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: pos, dir, centre, orig, new_pos
        real(kind=wp) :: t, semia, semib

        centre = vector(0.0, 0.0, 0.0)
        semia = 5._wp
        semib = 10._wp

        orig = vector(0.0, 0., -20.)
        dir = vector(0.0, 0., 1.)
        flag = intersectEllipse(orig, dir, t, centre, semia, semib)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, -5.0_wp, thr=1e-4_wp)
        if(allocated(error))return

        orig = vector(0.0, 0., -20.)
        dir = vector(0., 0., -1.)
        flag = intersectEllipse(orig, dir, t, centre, semia, semib)

        call check(error, flag, .false.)
        if(allocated(error))return
    end subroutine test_ellipse

    subroutine test_cone(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: pos, dir, centre, orig, new_pos
        real(kind=wp) :: t, radius, height

        centre = vector(0.0, 0.0, 0.0)
        radius = 10._wp
        height = 5._wp

        orig = vector(0.0, 0., -20.)
        dir = vector(0.0, 0., 1.)
        flag = intersectCone(orig, dir, t, centre, radius, height)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, 5.0_wp)
        if(allocated(error))return

        orig = vector(0.0, 0., -20.)
        dir = vector(0., 0., -1.)
        flag = intersectCone(orig, dir, t, centre, radius, height)

        call check(error, flag, .false.)
        if(allocated(error))return
    end subroutine test_cone

    subroutine test_plane(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: n, centre, orig, dir, new_pos
        real(kind=wp) :: t

        centre = vector(0.0, 0.0, 0.0)
        n = vector(0.0, 0.0, 1.0)
        orig = vector(0.0, 0.0, -10.0)
        dir = vector(0.0, 0.0, 1.0)

        flag = intersectPlane(n, centre, orig, dir, t)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, 0.0_wp)
        if(allocated(error))return

        orig = vector(0.0, 0., -20.)
        dir = vector(0., 0., -1.)
        flag = intersectPlane(n, centre, orig, dir, t)

        call check(error, flag, .false.)
        if(allocated(error))return
    end subroutine test_plane

    subroutine test_circle(error)

        type(error_type), allocatable, intent(out) :: error
        logical :: flag
        type(vector) :: n, centre, orig, dir, new_pos
        real(kind=wp) :: t, radius

        radius = 10.0_wp
        centre = vector(0.0, 0.0, 0.0)
        n = vector(0.0, 0.0, 1.0)

        orig = vector(0.0, 0.0, -10.0)
        dir = vector(0.0, 0.0, 1.0)

        flag = intersectCircle(n, centre, radius, orig, dir, t)
        new_pos = orig + t * dir

        call check(error, flag, .true.)
        if(allocated(error))return

        call check(error, new_pos%x, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%y, 0.0_wp)
        if(allocated(error))return

        call check(error, new_pos%z, 0.0_wp)
        if(allocated(error))return

        orig = vector(0.0, 0., -20.)
        dir = vector(0., 0., -1.)
        flag = intersectCircle(n, centre, radius, orig, dir, t)

        call check(error, flag, .false.)
        if(allocated(error))return
    end subroutine test_circle
end module testsgeometryMod