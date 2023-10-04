module testsParseMod

    use parse_mod
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, new_testsuite, testsuite_type, context_t
    use constants, only : wp

    implicit none

    private
    public :: parse_suite

    contains

    subroutine parse_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test input file", collect_suite1, context)&
                     ]

    end subroutine parse_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Test parsing of detectors", detectors), &
                new_unittest("Test parsing of spectra constant", spectra_const),&
                new_unittest("Test parsing of spectra spectral 1D", spectra_spectral_1D),&
                new_unittest("Test parsing of spectra spectral 2D", spectra_spectral_2D)&
                ]
    end subroutine collect_suite1

    subroutine spectra_const(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum

        filename = "res/test_spectra_const.toml"
        call parse_params(filename, packet, dects, spectrum, dict)

        select type(ptr => spectrum%p)
            class is(constant)    
                call check(error, ptr%value, 500.0_wp, "Wrong wavelength value in source spectrum_constant type!")
                if(allocated(error))return
        end select

    end subroutine spectra_const

    subroutine spectra_spectral_1D(error)

        use photonMod, only : photon
        use piecewiseMod, only : spectrum_t, piecewise1D
        use detectors, only : dect_array
        use tomlf, only: toml_table

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum

        filename = "res/test_spectra_const.toml"
        call parse_params(filename, packet, dects, spectrum, dict)
        select type(ptr => spectrum%p)
            class is(piecewise1D)    
                call check(error, size(ptr%cdf), 375, "Wrong CDF length in source spectrum_1D type!")
                if(allocated(error))return
                call check(error, size(ptr%array, 1), 375, "Wrong PDF length in source spectrum_1D type!")
                if(allocated(error))return
                call check(error, size(ptr%array, 2), 2, "Wrong PDF length in source spectrum_1D type!")
                if(allocated(error))return
        end select

    end subroutine spectra_spectral_1D

    subroutine spectra_spectral_2D(error)

        use detectors,    only : dect_array
        use photonMod,    only : photon
        use piecewiseMod, only : spectrum_t, piecewise2D
        use tomlf,        only : toml_table

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum

        filename = "res/test_spectra_const.toml"
        call parse_params(filename, packet, dects, spectrum, dict)
        select type(ptr => spectrum%p)
            class is(piecewise2D)    
                call check(error, size(ptr%cdf), 200*200, "Wrong CDF length in source spectrum_2D type!")
                if(allocated(error))return
                call check(error, ptr%cell_height, 0.5_wp, "Wrong cell height in source spectrum_2D type!")
                if(allocated(error))return
                call check(error, ptr%cell_width, 0.5_wp, "Wrong cell width in source spectrum_2D type!")
                if(allocated(error))return
        end select

    end subroutine spectra_spectral_2D

    subroutine detectors(error)

        use photonMod, only : photon
        use detectors, only : dect_array, circle_dect, camera, annulus_dect
        use piecewiseMod, only : spectrum_t
        use tomlf, only: toml_table

        type(error_type), allocatable, intent(out) :: error
        character(len=:),allocatable                  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum

        integer :: i

        filename = "res/test_dects.toml"
        call parse_params(filename, packet, dects, spectrum, dict)

        do i = 1, size(dects)
            select type(ptr => dects(i)%p)
            class is(circle_dect)
                call check(error, ptr%layer, 4, "Wrong layer value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%radius, 0.5_wp, "Wrong radius value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%nbins, 11, "Wrong nbins value in Circle dect!")
                if(allocated(error))return

                call check(error, ptr%pos%x, -1.0_wp, "Wrong pos%x value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%pos%y, 0.0_wp, "Wrong pos%y value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%pos%z, 0.0_wp, "Wrong pos%z value in Circle dect!")
                if(allocated(error))return

                call check(error, ptr%dir%x, -1.0_wp, "Wrong dir%x value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%dir%y, 0.0_wp, "Wrong dir%y value in Circle dect!")
                if(allocated(error))return
                call check(error, ptr%dir%z, 0.0_wp, "Wrong dir%z value in Circle dect!")
                if(allocated(error))return

            class is(annulus_dect)
                call check(error, ptr%layer, 3, "Wrong layer value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%r1, 0.5_wp, "Wrong r1 value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%r2, 1.0_wp, "Wrong r2 value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%nbins, 11, "Wrong nbins value in annulus dect!")
                if(allocated(error))return
                
                call check(error, ptr%pos%x, -1.0_wp, "Wrong pos%x value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%pos%y, 0.0_wp, "Wrong pos%y value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%pos%z, 0.0_wp, "Wrong pos%z value in annulus dect!")
                if(allocated(error))return

                call check(error, ptr%dir%x, -1.0_wp, "Wrong dir%x value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%dir%y, 0.0_wp, "Wrong dir%y value in annulus dect!")
                if(allocated(error))return
                call check(error, ptr%dir%z, 0.0_wp, "Wrong dir%z value in annulus dect!")
                if(allocated(error))return
            class is(camera)
                call check(error, ptr%layer, 2, "Wrong layer value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%pos%x, -1.0_wp, "Wrong pos%x value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%pos%y, -1.0_wp, "Wrong pos%y value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%pos%z, -1.0_wp, "Wrong pos%z value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%p2%x, 0.0_wp, "Wrong p2%x value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%p2%y, 2.0_wp, "Wrong p2%y value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%p2%z, 0.0_wp, "Wrong p2%z value in camera dect!")
                if(allocated(error))return
                call check(error, ptr%nbinsX, 11, "Wrong nbins value in camera dect!")
                if(allocated(error))return
            end select
        end do

    end subroutine detectors
end module testsParseMod