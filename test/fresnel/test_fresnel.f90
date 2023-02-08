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
                     new_unittest("simple_reflect", simple_reflect), &
                     new_unittest("complex_refract", complex_refract), &
                     new_unittest("complex_reflect", complex_reflect) &
                ]
    end subroutine collect_suite1

    subroutine complex_refract(error)

        use vector_class
        use surfaces
        use utils, only: rad2deg, deg2rad
        use constants, only : pi
        use random, only : init_rng

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: I, N
        logical :: rflag
        real(kind=wp) :: ri, n1, n2, theta, real_value, phi, computed_theta
        integer :: j, counter

        call init_rng(input_seed=spread(123456789, 1, 8))

        counter = 0
    
        theta = deg2rad(180._wp + 45.0_wp)
        phi = 0._wp
        
        N = vector(0._wp, 0._wp, 1._wp)
        n1 = 1._wp
        n2 = 1.33_wp
        real_value = asin(n1/ n2 * sin(deg2rad(45._wp)))
        
        ! repeat calculation as relies on random chance
        ! check it falls with allowed range (less than 6e-2%)
        do j = 1, 1000000
            I%x = abs(sin(theta) * cos(phi))
            I%y = sin(theta) * sin(phi)
            I%z = cos(theta)
            I = I%magnitude()
            call reflect_refract(I, N, n1, n2, rflag, ri)
            computed_theta = pi - acos(I%z)
            if(abs(computed_theta - real_value) < 1.0e-10_wp)counter = counter+1
        end do
        call check(error, counter/1000000.0_wp, 1._wp-ri, thr=5e-4_wp)
        if(allocated(error))return

    end subroutine complex_refract

    subroutine complex_reflect(error)

        use vector_class
        use surfaces
        use utils, only: deg2rad

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: I, N, I_before
        logical :: rflag
        real(kind=wp) :: ri, n1, n2, theta, phi
        integer :: j, counter

        counter = 0

        theta = deg2rad(45._wp)
        phi = 0._wp

        I%x = sin(theta) * cos(phi)
        I%y = sin(theta) * sin(phi)
        I%z = cos(theta)
        I = I%magnitude()

        I_before = I

        N = vector(0._wp, 0._wp, 1._wp)
        n1 = 1.33_wp
        n2 = 1.00_wp

        do j = 1, 1000000
            I%x = abs(sin(theta) * cos(phi))
            I%y = sin(theta) * sin(phi)
            I%z = cos(theta)
            I = I%magnitude()
            call reflect_refract(I, N, n1, n2, rflag, ri)
            if(I%z == -1._wp*I_before%z .and. rflag)counter = counter + 1
        end do
        call check(error, counter/1000000.0_wp, ri, thr=5e-4_wp)
        if(allocated(error))return
    end subroutine complex_reflect

    subroutine simple_refract(error)

        use vector_class
        use surfaces
        use utils, only: rad2deg, deg2rad
        use constants, only : pi

        type(error_type), allocatable, intent(out) :: error

        type(vector) :: I, N
        logical :: rflag
        real(kind=wp) :: ri, n1, n2, theta, real_value, phi

        ! normal incidence so no chance of reflection
        theta = deg2rad(180._wp)
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
        real_value = 0.0_wp

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