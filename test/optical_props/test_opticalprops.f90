module testsOpticalPropMod

    use opticalProperties
    use testdrive, only : new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type, context_t, test_failed
    use constants, only : wp

    implicit none

    private
    public :: OpticalProp_suite

    contains

    subroutine OpticalProp_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test Spectral", collect_suite1, context),&
                      new_testsuite("Test Mono", collect_suite2, context)&
                     ]

    end subroutine OpticalProp_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Spectral", test_spectral)&
                ]

    end subroutine collect_suite1

    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Mono", test_mono)&
                ]

    end subroutine collect_suite2

    subroutine test_mono(error)

        type(mono) :: optProp
        type(error_type), allocatable, intent(out) :: error
        real(kind=wp) :: mus, mua, hgg, n, wave

        mus = 10.0_wp
        mua = 0.1_wp
        hgg = 0.9_wp
        n = 1.35_wp

        optProp = mono(mus, mua, hgg, n)

        call optProp%update(wave)

        mus = optProp%mus
        call check(error, mus, 10._wp, thr=0.05_wp)
        if(allocated(error))return

        mua = optProp%mua
        call check(error, mua, 0.1_wp, thr=0.05_wp)
        if(allocated(error))return

        hgg = optProp%hgg
        call check(error, hgg, 0.9_wp, thr=0.05_wp)
        if(allocated(error))return

        n = optProp%n
        call check(error, n, 1.35_wp, thr=0.05_wp)
        if(allocated(error))return

        call check(error, wave, 0.0_wp, thr=0.05_wp)
        if(allocated(error))return

    end subroutine test_mono

    subroutine test_spectral(error)

        use random, only : init_rng
        use utils,  only : str

        type(error_type), allocatable, intent(out) :: error
        type(spectral) :: optProp
        integer ::  u, i
        real(kind=wp), allocatable :: flux(:,:), mua_a(:, :), hgg_a(:, :), n_a(:, :), mus_a(:, :)
        real(kind=wp) :: wave, hgg, mua, n, mus
        
        call init_rng(spread(1234569, 1, 8), .false.)

        allocate(mua_a(10, 2))
        mua_a(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        mua_a(:, 2) = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

        allocate(flux(10, 2))
        flux(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        flux(:, 2) = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        allocate(n_a(10, 2))
        n_a(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        n_a(:, 2) = [1.0, 1.5, 1.5, 1.0, 1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

        allocate(mus_a(10, 2))
        mus_a(:, 1) = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        mus_a(:, 2) = [0.0, 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.0]

        allocate(hgg_a(3, 2))
        hgg_a(:, 1) = [100, 450, 900]
        hgg_a(:, 2) = [0.9, 0.9, 0.9]

        optProp = spectral(mus_a, mua_a, hgg_a, n_a, flux)

        do i = 1, 10000
            call optProp%update(wave)
            if(wave < 100 .or. wave > 1000)then
                call test_failed(error, "Spectral check failed!", "Expected a wavelength between [100, 1000]! Got "//str(wave))
            end if

            hgg = optProp%hgg
            call check(error, hgg, 0.9_wp, thr=0.05_wp)
            if(allocated(error))return
    
            mua = optProp%mua
            call check(error, mua, 0.1_wp, thr=0.05_wp)
            if(allocated(error))return

            n = optProp%n
            if(n < 1.0 .or. n > 2.2)then
                call test_failed(error, "Spectral check failed!", "Expected a refractive between [1, 2.2]! Got "//str(n))
            end if
            if(allocated(error))return

            mus = optProp%mus
            if(mus < 0.0 .or. mus > 4.0)then
                call test_failed(error, "Spectral check failed!", "Expected mus between [0.0, 4.0]! Got "//str(mus))
            end if
            if(allocated(error))return
        end do
    end subroutine test_spectral
end module testsOpticalPropMod