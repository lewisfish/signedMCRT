module parse_sourcesMod

    use constants, only : wp
    use parse_HelpersMod
    use parse_SpectrumMod
    use vector_class

    use tomlf

    implicit none
    
    private
    public :: parse_source

contains
    
    subroutine parse_source(table, packet, dict, spectrum, context, error)
    !! Parse sources
    !! any updates here MUST be reflected in docs/config.md
        use sim_state_mod, only : state
        use photonmod
        use piecewiseMod
        use tomlf_error
        
        !> Input Toml table
        type(toml_table),  intent(inout) :: table
        !> Dictonary used to store metadata
        type(toml_table),  intent(inout) :: dict
        !> Photon packet. Used to store information to save computation
        type(photon),      intent(out)   :: packet
        !> Spectrum type.
        type(spectrum_t),  intent(out)   :: spectrum
        !> Context handle for error reporting
        type(toml_context) :: context
        !> Error message
        type(toml_error), allocatable, intent(out) :: error

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        type(vector) :: poss, dirr
        real(kind=wp) :: dir(3), pos(3), corners(3, 3), radius, beta, rlo, rhi
        integer :: i, nlen, origin
        character(len=1) :: axis(3)
        character(len=:), allocatable :: direction, annulus_type

        axis = ["x", "y", "z"]
        pos = 0._wp
        dir = 0._wp
        corners = reshape((/ -1._wp, -1._wp, 1._wp, &
                              2._wp,  0._wp, 0._wp, &
                              0._wp,  2._wp, 0._wp /), &
                           shape(corners), order=[2, 1])

        call get_value(table, "source", child, requested=.false.)
        if(associated(child))then
            call get_value(child, "name", state%source, "point")
            call get_value(child, "nphotons", state%nphotons, 1000000)

            poss = get_vector(child, "position", error, context)
            if(allocated(error))then
                if(state%source == "point")then
                    deallocate(error)
                    call make_error(error, &
                    context%report("Point source needs a position!", origin, "Need vector of size 3 for position"), -1)
                end if
                return
            end if

            children => null()
            
            call get_value(child, "direction", children, requested=.false., origin=origin)
            if(associated(children))then
                if(state%source == "point")then
                    print'(a)',context%report(&
                    "Point source needs no direction!!", origin, level=toml_level%warning)
                end if
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a vector of size 3 for direction", origin, "expected vector of size 3"), -1)
                    return
                end if
                if(state%source == "circular")then
                    print'(a)',context%report(&
                    "Direction not yet fully tested for source type Circular. Results may not be accurate!", origin,&
                     level=toml_level%warning)
                end if
                do i = 1, len(children)
                    call get_value(children, i, dir(i))
                end do
                dirr%x = dir(1)
                dirr%y = dir(2)
                dirr%z = dir(3)
            else
                call get_value(child, "direction", direction, origin=origin)
                if(allocated(direction))then
                    if(state%source == "point")then
                        print'(a)',context%report(&
                        "Point source needs no direction!!", origin, level=toml_level%warning)
                    end if
    
                    select case(direction)
                    case("x")
                        dirr = vector(1._wp, 0._wp, 0._wp)
                    case("-x")
                        dirr = vector(-1._wp, 0._wp, 0._wp)                
                    case("y")
                        dirr = vector(0._wp, 1._wp, 0._wp)
                    case("-y")
                        dirr = vector(0._wp, -1._wp, 0._wp)
                    case("z")
                        dirr = vector(0._wp, 0._wp, 1._wp)
                    case("-z")
                        dirr = vector(0._wp, 0._wp, -1._wp)
                    case default
                        call make_error(error, context%report("Direction needs a cardinal direction i.e x, y, or z", origin, &
                                                "Expected cardinal direction"), -1)
                        return 
                    end select
                elseif(state%source /= "point")then
                    call make_error(error, context%report("Need to specify direction for source type!", origin, &
                                              "No direction specified"), -1)
                    return
                end if
            end if

            children => null()
            
            call get_value(child, "point1", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 1))
                    call set_value(dict, "pos1%"//axis(i), corners(i,1))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point1 variable", origin, "expected point1 variable"), -1)
                    return
                end if
            end if

            call get_value(child, "point2", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 2))
                    call set_value(dict, "pos2%"//axis(i), corners(i,2))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point2 variable", origin, "expected point2 variable"), -1)
                    return
                end if
            end if

            call get_value(child, "point3", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 3))
                    call set_value(dict, "pos3%"//axis(i), corners(i,3))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point3 variable", origin, "expected point3 variable"), -1)
                    return
                end if
            end if
            call get_value(child, "radius", radius, 0.5_wp)
            call set_value(dict, "radius", radius)

            ! parameters for annulus beam type
            call get_value(child, "beta", beta, 5._wp)
            call set_value(dict, "beta", beta)

            call get_value(child, "radius_hi", rhi, 0.6_wp)
            call set_value(dict, "rhi", rhi)

            call get_value(child, "annulus_type", annulus_type, "gaussian")
            call set_value(dict, "annulus_type", annulus_type)

            ! parse spectrum
            call parse_spectrum(child, spectrum, dict, context, error)
            if(allocated(error))return
        else
            call make_error(error, context%report("Simulation needs Source table", origin, "Missing source table"), -1)
            return
        end if

        call set_photon(poss, dirr)
        packet = photon(state%source)
        packet%pos = poss
        packet%nxp = dirr%x
        packet%nyp = dirr%y
        packet%nzp = dirr%z

    end subroutine parse_source
end module parse_sourcesMod