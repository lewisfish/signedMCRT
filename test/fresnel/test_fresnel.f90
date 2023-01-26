module testsFresnelMod

    use surfaces
    use testdrive, only : new_unittest, unittest_type, error_type, check, testsuite_type, new_testsuite
    use constants, only : wp

    implicit none

    contains

    subroutine Fresnel_suite(testsuites)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)

        testsuites = [new_testsuite("Fresnel", collect_suite1)&
                     ]

    end subroutine Fresnel_suite


    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                     new_unittest("simple_refract", simple_refract), &
                     new_unittest("simple_reflect", simple_reflect) &
                ]
    end subroutine collect_suite1

    subroutine simple_refract(error)

        use vector_class
        use surfaces
        use utils, only: rad2deg, deg2rad
        use constants, only : pi

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: I, N
        logical :: rflag
        real(kind=wp) :: ri, n1, n2, theta, real_value, phi


        theta = deg2rad(180._wp + 45._wp)
        phi = 0._wp

        I%x = abs(sin(theta) * cos(phi))
        I%y = sin(theta) * sin(phi)
        I%z = cos(theta)
        I = I%magnitude()
        
        N = vector(0._wp, 0._wp, 1._wp)
        n1 = 1._wp
        n2 = 1.33_wp

        call reflect_refract(I, N, n1, n2, rflag, ri)

        theta = pi - acos(I%z)
        real_value = asin(n1/ n2 * sin(deg2rad(45._wp)))
        call check(error, theta, real_value, thr=1.0e-10_wp)
        if(allocated(error))return

    end subroutine simple_refract

    subroutine simple_reflect(error)

        use vector_class
        use surfaces
        use utils, only: deg2rad

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: I, N, I_before
        logical :: rflag
        real(kind=wp) :: ri, n1, n2, theta, phi

        theta = deg2rad(50._wp)
        phi = 0._wp

        I%x = sin(theta) * cos(phi)
        I%y = sin(theta) * sin(phi)
        I%z = cos(theta)
        I = I%magnitude()

        I_before = I

        N = vector(0._wp, 0._wp, 1._wp)
        n1 = 1.33_wp
        n2 = 1.00_wp

        call reflect_refract(I, N, n1, n2, rflag, ri)

        call check(error, I%x, I_before%x)
        if(allocated(error))return
        call check(error, I%y, I_before%y)
        if(allocated(error))return
        call check(error, I%z, -1._wp*I_before%z)
        if(allocated(error))return

    end subroutine simple_reflect

end module testsFresnelMod