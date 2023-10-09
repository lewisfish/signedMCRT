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

        testsuites = &![new_testsuite("Test input file: Success", collect_suite1, context),&
                      [new_testsuite("Test input file: Fail", collect_suite2, context)&
                     ]

    end subroutine parse_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Test parsing of detectors", detector_test), &
                new_unittest("Test parsing of spectra constant", spectra_const),&
                new_unittest("Test parsing of spectra spectral 1D", spectra_spectral_1D),&
                new_unittest("Test parsing of spectra spectral 2D", spectra_spectral_2D)&
                ]
    end subroutine collect_suite1


    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Test failing: No file", test_non_existant), &
                new_unittest("Test failing: No Source Table", test_no_source_table), &
                new_unittest("Test failing: Not valid Source Table 1", test_non_valid_src_table_1), &
                new_unittest("Test failing: Not valid Source Table 2", test_non_valid_src_table_2), &
                new_unittest("Test failing: Not valid Source Table 3", test_non_valid_src_table_3) &
                ! new_unittest("Test failing: Not valid Source Table 4", test_non_valid_src_table_4), &
                ! new_unittest("Test failing: Not valid Source Table 5", test_non_valid_src_table_5), &
                ! new_unittest("Test failing: Not valid Source Table 6", test_non_valid_src_table_6), &
                ! new_unittest("Test failing: Not valid Source Table 7", test_non_valid_src_table_7), &
                ! new_unittest("Test failing: Not valid Source Table 8", test_non_valid_src_table_8), &
                ! new_unittest("Test failing: Not valid Source Table 9", test_non_valid_src_table_9), &
                ! new_unittest("Test failing: Not valid Source Table 10", test_non_valid_src_table_10), &
                ! new_unittest("Test failing: Not valid Source Table 12", test_non_valid_src_table_11), &
                ! new_unittest("Test failing: Not valid Source Table 11", test_non_valid_src_table_12), &
                ! new_unittest("Test failing: Not valid detector type", test_non_valid_dect), &
                ! new_unittest("Test failing: Not valid Annulus dect", test_non_valid_annulus), &
                ! new_unittest("Test failing: Not valid spectrum type", test_non_valid_spectrum)&
                ]
    end subroutine collect_suite2


    subroutine test_non_existant(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        filename = "test/parse/does_not_exist.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)
        call check(error, allocated(err), .true.)
        if (allocated(error))return

    end subroutine test_non_existant

    subroutine test_non_valid_src_table_1(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err
        integer :: u

        filename = "test/parse/test_fail1.toml"
        open(newunit=u, file=filename)
        write(u,'(a)') "[source]"
        write(u,'(a)') "position=[0.0, 0.0]"
        close(u)

        call parse_params(filename, packet, dects, spectrum, dict, err)
        call check(error, allocated(err), .true.)
        if (allocated(error))return

    end subroutine test_non_valid_src_table_1

    subroutine test_non_valid_src_table_2(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err
        integer :: u

        filename = "test/parse/test_fail1.toml"
        open(newunit=u, file=filename)
        write(u,'(a)') "[source]"
        write(u,'(a)') "name='point'"
        close(u)

        call parse_params(filename, packet, dects, spectrum, dict, err)
        call check(error, allocated(err), .true.)
        if (allocated(error))return

    end subroutine test_non_valid_src_table_2

    subroutine test_non_valid_src_table_3(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err
        integer :: u

        filename = "test/parse/test_fail1.toml"
        open(newunit=u, file=filename)
        write(u,'(a)') "[source]"
        write(u,'(a)') "name='uniform'"
        write(u,'(a)') "position=[0.0, 0.0, 0.0]"
        write(u,'(a)') "direction=[0.0, 0.0]"
        close(u)

        call parse_params(filename, packet, dects, spectrum, dict, err)
        call check(error, allocated(err), .true.)
        if (allocated(error))return

    end subroutine test_non_valid_src_table_3

    subroutine test_no_source_table(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        filename = "test/parse/test_fail1.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)
        call check(error, allocated(err), .true.)
        if (allocated(error))return

    end subroutine test_no_source_table

    subroutine spectra_const(error)

        use photonMod, only : photon
        use detectors, only : dect_array
        use piecewiseMod, only : spectrum_t, constant
        use tomlf, only: toml_table, toml_error
        
        type(error_type), allocatable, intent(out) :: error
        
        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        filename = "res/test_spectra_const.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)

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
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        filename = "res/test_spectra_1D.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)
        select type(ptr => spectrum%p)
            class is(piecewise1D)    
                call check(error, size(ptr%cdf), 376, "Wrong CDF length in source spectrum_1D type!")
                if(allocated(error))return
                call check(error, size(ptr%array, 1), 376, "Wrong PDF length in source spectrum_1D type!")
                if(allocated(error))return
                call check(error, size(ptr%array, 2), 2, "Wrong PDF length in source spectrum_1D type!")
                if(allocated(error))return
        end select

    end subroutine spectra_spectral_1D

    subroutine spectra_spectral_2D(error)

        use constants,    only : resdir
        use detectors,    only : dect_array
        use photonMod,    only : photon
        use piecewiseMod, only : spectrum_t, piecewise2D
        use tomlf,        only : toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error

        character(len=:),allocatable  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        resdir = ""
        filename = "res/test_spectra_2D.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)
        select type(ptr => spectrum%p)
            class is(piecewise2D) 
                call check(error, size(ptr%cdf), 256*256, "Wrong CDF length in source spectrum_2D type!")
                if(allocated(error))return
                call check(error, ptr%cell_height, 0.5_wp, "Wrong cell height in source spectrum_2D type!")
                if(allocated(error))return
                call check(error, ptr%cell_width, 0.5_wp, "Wrong cell width in source spectrum_2D type!")
                if(allocated(error))return
        end select

    end subroutine spectra_spectral_2D

    subroutine detector_test(error)

        use photonMod, only : photon
        use detectors, only : dect_array, circle_dect, camera, annulus_dect
        use piecewiseMod, only : spectrum_t
        use tomlf, only: toml_table, toml_error

        type(error_type), allocatable, intent(out) :: error
        character(len=:),allocatable                  :: filename
        type(toml_table)              :: dict
        type(photon)                  :: packet
        type(dect_array), allocatable :: dects(:)
        type(spectrum_t)              :: spectrum
        type(toml_error), allocatable :: err

        integer :: i

        filename = "res/test_dects.toml"
        call parse_params(filename, packet, dects, spectrum, dict, err)

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

    end subroutine detector_test
end module testsParseMod