module parseSpectrumMod

    use constants, only : wp
    use parseHelpers
    use tomlf
    use tomlf_error, only: make_error
    use vector_class

    implicit none
    
    private
    public :: parse_spectrum

contains
    
    subroutine parse_spectrum(table, spectrum, dict, context, error)
    !! Parse spectrums to be used
    ! TODO seperate out each case to seperate functions.
    ! TODO add spectra type to input optical properties
    ! handle all possible errors
    ! document code and update config.md
        use piecewiseMod
        use stdlib_io, only: loadtxt
        use constants, only : resdir, sp

        use stb_image_mod
        use, intrinsic :: iso_c_binding

        type(toml_table),              pointer       :: table
        !> Polymorphic spectrum type
        type(spectrum_t),              intent(out)   :: spectrum
        !> Dictionary of metadata
        type(toml_table),              intent(inout) :: dict
        !> Context handle for error reporting.
        type(toml_context),            intent(inout) :: context
        !> Error messages that are raised in the routine
        type(toml_error), allocatable, intent(out)   :: error


        type(toml_array), pointer :: children
        integer :: origin, nlen, i, err, width, height, n_channels,u
        integer, allocatable :: image(:,:,:)
        type(constant), save, target :: const
        type(piecewise1D), save, target :: OneD
        type(piecewise2D), save, target :: TwoD
        character(len=:), allocatable :: stype, sfile, filetype
        real(kind=wp) :: wavelength, cellsize(2)
        real(kind=wp), allocatable :: array(:,:)
        real(kind=sp), allocatable :: array_sp(:,:)

        call get_value(table, "spectrum_type", stype, "constant", origin=origin)
        select case(stype)
            case("constant")
                call get_value(table, "wavelength", wavelength, 500.0_wp)
                const = constant(wavelength)
                allocate(spectrum%p, source=const)
                spectrum%p => const
            case("1D")
                allocate(spectrum%p, source=OneD)
                call get_value(table, "spectrum_file", sfile)
                call loadtxt("res/"//sfile, array_sp)
                array = array_sp
                deallocate(array_sp)
                OneD = piecewise1D(array)
                allocate(spectrum%p, source=OneD)
                spectrum%p => OneD
            case("2D")
                allocate(spectrum%p, source=TwoD)
                call get_value(table, "spectrum_file", sfile)

                call get_value(table, "cell_size", children, requested=.true., origin=origin)
                if(associated(children))then
                    nlen = len(children)
                    if(nlen /= 2)then
                        call make_error(error,&
                        context%report("Need a vector of size 2 for cell_size", origin, "expected vector of size 2"), -1)
                        return
                    end if
                    do i = 1, len(children)
                        call get_value(children, i, cellsize(i))
                    end do
                else
                    call make_error(error,&
                    context%report("Need a vector of size 2 for cell_size", origin, "expected vector of size 2"), -1)
                    return
                end if

                filetype = sfile(len(sfile)-2:)
                select case(filetype)
                case("png")
                    err = stbi_info(trim(resdir)//trim(sfile)//c_null_char, width, height, n_channels)
                    if(err == 0)then
                        call make_error(error, "Error reading file: "//trim(sfile)//" "//stbi_failure_reason(), -1)
                        return
                    end if
                    image = stbi_load(trim(resdir)//trim(sfile)//c_null_char, width, height, n_channels, 0)
                    allocate(array(size(image, 1), size(image, 2)))
                    array = image(:,:,1)

                    deallocate(image)

                case("dat")
                    call loadtxt(resdir//trim(sfile), array)
                case("txt")
                    call loadtxt(resdir//trim(sfile), array)
                case default
                    print'(2a)', "Unknown spectrum file type:", filetype
                end select
                TwoD = piecewise2D(cellsize(1), cellsize(2), array)
                allocate(spectrum%p, source=TwoD)
                spectrum%p => TwoD
            case default
                call make_error(error,&
                context%report("Not a valid spectrum type!", origin, "expected one of either ['constant', '1D', '2D']"),-1)
                return
        end select
    end subroutine parse_spectrum
end module parseSpectrumMod