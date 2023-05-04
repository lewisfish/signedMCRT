module parse_mod

    use tomlf
    use constants, only : wp
    use vector_class

    implicit none

    private
    public :: parse_params

    contains

    subroutine parse_params(filename, packet, dects, spectrum, dict)
    
        use photonmod
        use detector_mod, only : dect_array
        use piecewiseMod
    
        character(*),      intent(IN)    :: filename
        type(toml_table),  intent(INOUT) :: dict
        type(photon),      intent(OUT)   :: packet
        type(dect_array), allocatable, intent(out) :: dects(:)
        type(spectrum_t), intent(out) :: spectrum

        type(toml_table), allocatable :: table
        type(toml_context) :: context
        type(toml_error), allocatable :: error

        call toml_load(table, trim(filename), context=context, error=error)
        if(allocated(error))then
            print'(a)',error%message
            stop 1
        end if

        call parse_source(table, packet, dict, spectrum, context)
        call parse_grid(table, dict)
        call parse_geometry(table, dict)
        call parse_detectors(table, dects, context)
        call parse_output(table, dict)
        call parse_simulation(table)

    end subroutine parse_params
    
    subroutine parse_detectors(table, dects, context)

        use detector_mod
        use sim_state_mod, only : state

        type(toml_table) :: table
        type(dect_array), allocatable :: dects(:)
        type(toml_context), intent(in) :: context

        type(toml_array), pointer :: array
        type(toml_table), pointer :: child
        character(len=:), allocatable :: dect_type
        type(circle_dect), target, save, allocatable :: dect_c(:)
        type(annulus_dect), target, save, allocatable :: dect_a(:)
        type(camera), target, save, allocatable :: dect_cam(:)
        integer :: i, c_counter, a_counter, cam_counter, j, k,origin

        c_counter = 0
        a_counter = 0
        cam_counter = 0
        call get_value(table, "detectors", array)
        allocate(dects(len(array)))

        do i = 1, len(array)
            call get_value(array, i, child)
            call get_value(child, "type", dect_type, origin=origin)
            select case(dect_type)
            case default
                print'(a)',context%report("Invalid detector type. Valid types are [circle, annulus, camera]", &
                           origin, "expected valid detector type")
                stop 1
            case("circle")
                c_counter = c_counter + 1
            case("annulus")
                a_counter = a_counter + 1
            case("camera")
                cam_counter = cam_counter + 1
            end select
        end do

        if(c_counter > 0)allocate(dect_c(c_counter))
        if(a_counter > 0)allocate(dect_a(a_counter))
        if(cam_counter > 0)allocate(dect_cam(cam_counter))
        c_counter = 1
        a_counter = 1
        cam_counter = 1
        state%trackHistory=.false.
        do i = 1, len(array)
            call get_value(array, i, child)
            call get_value(child, "type", dect_type)
            call get_value(child, "historyFileName", state%historyFilename, "photPos.obj")
            select case(dect_type)
            case("circle")
                call handle_circle_dect(child, dect_c, c_counter, context)
            case("annulus")
                call handle_annulus_dect(child, dect_a, a_counter, context)
            case("camera")
                call handle_camera(child, dect_cam, cam_counter, context)
            end select
        end do

        do i = 1, c_counter-1
            allocate(dects(i)%p, source=dect_c(i))
            dects(i)%p => dect_c(i)
        end do

        do j = 1, a_counter-1
            allocate(dects(j+i-1)%p, source=dect_a(j))
            dects(j+i-1)%p => dect_a(j)
        end do

        do k = 1, cam_counter-1
            allocate(dects(j+i+k-2)%p, source=dect_cam(k))
            dects(j+i+k-2)%p => dect_cam(k)
        end do

        if(.not. allocated(state%historyFilename))state%historyFilename="photPos.obj"

    end subroutine parse_detectors

    subroutine handle_camera(child, dects, counts, context)

        use detector_mod
        use sim_state_mod, only : state

        type(toml_table), pointer, intent(in)    :: child
        type(camera),              intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins
        real(kind=wp) :: maxval
        type(vector)  :: p1, p2, p3
        logical       :: trackHistory

        p1 = get_vector(child, "p1", default=vector(-1.0, -1.0, -1.0), context=context)
        p2 = get_vector(child, "p2", default=vector(2.0, 0.0, 0.0), context=context)
        p3 = get_vector(child, "p3", default=vector(0.0, 2.0, 0.0), context=context)
 
        call get_value(child, "layer", layer, 1)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)error stop "Track history currently incompatable with OpenMP!"
#endif
        dects(counts) = camera(p1, p2, p3, layer, nbins, maxval, trackHistory)
        counts = counts + 1

    end subroutine handle_camera

    subroutine handle_circle_dect(child, dects, counts, context)

        use detector_mod
        use sim_state_mod, only : state

        type(toml_table), pointer, intent(in)    :: child
        type(circle_dect),         intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins
        real(kind=wp) :: maxval, radius
        type(vector)  :: pos, dir
        logical       :: trackHistory

        pos = get_vector(child, "position", context=context)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context)
        dir = dir%magnitude()
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)error stop "Track history currently incompatable with OpenMP!"
#endif
        dects(counts) = circle_dect(pos, dir, layer, radius, nbins, maxval, trackHistory)
        counts = counts + 1

    end subroutine handle_circle_dect

    subroutine handle_annulus_dect(child, dects, counts, context)

        use detector_mod
        use sim_state_mod, only : state

        type(toml_table), pointer, intent(in)    :: child
        type(annulus_dect),        intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins, origin
        real(kind=wp) :: maxval, radius1, radius2
        type(vector)  :: pos
        logical       :: trackHistory

        pos = get_vector(child, "position", context=context)
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius1)
        call get_value(child, "radius2", radius2, origin=origin)
        if(radius2 <= radius1)then
            print'(a)',context%report("Radii are invalid", origin, "Expected radius2 > radius 1")
            stop 1
        end if
            call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)error stop "Track history currently incompatable with OpenMP!"
#endif
        dects(counts) = annulus_dect(pos, layer, radius1, radius2, nbins, maxval, trackHistory)
        counts = counts + 1
    end subroutine handle_annulus_dect

    subroutine parse_spectrum(table, spectrum, dict, context)
    ! TODO seperate out each case to seperate functions.
    ! handle all possible errors
    ! document code and update config.md
        use piecewiseMod
        use stdlib_io, only: loadtxt
        use constants, only : fileplace

        type(toml_table),  intent(INOUT) :: dict
        type(toml_table), pointer :: table

        type(toml_context) :: context
        type(spectrum_t), intent(out) :: spectrum

        type(toml_array), pointer :: children
        integer :: origin, nlen, i
        type(constant), save, target :: const
        type(piecewise1D), save, target :: OneD
        type(piecewise2D), save, target :: TwoD
        character(len=:), allocatable :: stype, sfile
        real(kind=wp) :: wavelength, cellsize(2)
        real(kind=wp), allocatable :: array(:,:)

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
                call loadtxt("res/"//sfile, array)
                OneD = piecewise1D(array)
                allocate(spectrum%p, source=OneD)
                spectrum%p => OneD
            case("2D")
                allocate(spectrum%p, source=TwoD)
                call get_value(table, "spectrum_file", sfile)

                call get_value(table, "cell_size", children, requested=.false., origin=origin)
                if(associated(children))then
                    nlen = len(children)
                    if(nlen /= 2)then
                        print'(a)',context%report("Need a vector of size 2 for cell_size", origin, "expected vector of size 2")
                        stop 1
                    end if
                    do i = 1, len(children)
                        call get_value(children, i, cellsize(i))
                    end do
                end if

                call loadtxt(fileplace//sfile, array)
                TwoD = piecewise2D(cellsize(1), cellsize(2), array)
                allocate(spectrum%p, source=TwoD)
                spectrum%p => TwoD
            case default
                print'(a)',context%report("Not a valid spectrum type!", origin, "expected one of either ['constant', '1D', '2D']")
                stop 1
        end select


    end subroutine parse_spectrum

    subroutine parse_source(table, packet, dict, spectrum, context)
    ! Parse source table
    ! any updates here MUST be reflected in docs/config.md
        use sim_state_mod, only : state
        use photonmod
        use piecewiseMod
        
        type(toml_table),  intent(INOUT) :: table, dict
        type(photon),      intent(OUT)   :: packet
        type(spectrum_t), intent(out) :: spectrum
        type(toml_context) :: context

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

            call get_value(child, "position", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a vector of size 3 for position", origin, "expected vector of size 3 1")
                    stop 1
                end if
                do i = 1, len(children)
                    call get_value(children, i, pos(i))
                end do
            else
                if(state%source == "point")then
                    print'(a)',context%report(&
                    "Point source needs a position!", origin, "Need vector of size 3 for position")
                    stop 1
                end if
            end if
            poss = vector(pos(1), pos(2), pos(3))

            children => null()
            
            call get_value(child, "direction", children, requested=.false., origin=origin)
            if(associated(children))then
                if(state%source == "point")then
                    print'(a)',context%report(&
                    "Point source needs no direction!!", origin, level=toml_level%warning)
                end if
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a vector of size 3 for direction", origin, "expected vector of size 3 2")
                    stop 1
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
                        print'(a)',context%report("Direction needs a cardinal direction i.e x, y, or z", origin, &
                                                "Expected cardinal direction")
                        stop 1 
                    end select
                elseif(state%source /= "point")then
                    print'(a)',context%report("Need to specify direction for source type!", origin, &
                                              "No direction specified")
                    stop 1
                end if
            end if

            children => null()
            
            call get_value(child, "point1", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a matrix row for points", origin, "expected matrix row of size 3")
                    stop 1
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 1))
                    call set_value(dict, "pos1%"//axis(i), corners(i,1))
                end do
            else
                if(state%source == "uniform")then
                    print'(a)',context%report("Uniform source requires point1 variable", origin, "expected point1 variable")
                    stop 1
                end if
            end if

            call get_value(child, "point2", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a matrix row for points", origin, "expected matrix row of size 3")
                    stop 1
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 2))
                    call set_value(dict, "pos2%"//axis(i), corners(i,2))
                end do
            else
                if(state%source == "uniform")then
                    print'(a)',context%report("Uniform source requires point2 variable", origin, "expected point2 variable")
                    stop 1
                end if
            end if

            call get_value(child, "point3", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a matrix row for points", origin, "expected matrix row of size 3")
                    stop 1
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 3))
                    call set_value(dict, "pos3%"//axis(i), corners(i,3))
                end do
            else
                if(state%source == "uniform")then
                    print'(a)',context%report("Uniform source requires point3 variable", origin, "expected point3 variable")
                    stop 1
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
            call parse_spectrum(child, spectrum, dict, context)
        else
            print'(a)',context%report("Simulation needs Source table", origin, "Missing source table")
            stop 1
        end if

        call set_photon(poss, dirr)
        packet = photon(state%source)
        packet%pos = poss
        packet%nxp = dirr%x
        packet%nyp = dirr%y
        packet%nzp = dirr%z

    end subroutine parse_source

    subroutine parse_grid(table, dict)

        use sim_state_mod, only : state
        use gridMod,       only : init_grid 
        
        type(toml_table),  intent(INOUT) :: table, dict

        type(toml_table), pointer     :: child
        integer                       :: nxg, nyg, nzg
        real(kind=wp)                 :: xmax, ymax, zmax
        character(len=:), allocatable :: units

        call get_value(table, "grid", child)

        if(associated(child))then
            call get_value(child, "nxg", nxg, 200)
            call get_value(child, "nyg", nyg, 200)
            call get_value(child, "nzg", nzg, 200)
            call get_value(child, "xmax", xmax, 1.0_wp)
            call get_value(child, "ymax", ymax, 1.0_wp)
            call get_value(child, "zmax", zmax, 1.0_wp)
            call get_value(child, "units", units, "cm")
            call set_value(dict, "units", units)
        else
            error stop "Need grid table in input param file"
        end if

        state%grid = init_grid(nxg, nyg, nzg, xmax, ymax, zmax)

    end subroutine parse_grid

    subroutine parse_geometry(table, dict)

        use sim_state_mod, only : state
                
        type(toml_table),  intent(INOUT) :: table, dict
        
        type(toml_table), pointer :: child
        real(kind=wp)             :: tau, musb, musc, muab, muac, hgg
        integer                   :: num_spheres

        call get_value(table, "geometry", child)

        if(associated(child))then
            call get_value(child, "geom_name", state%experiment, "sphere")
            call get_value(child, "tau", tau, 10._wp)
            call set_value(dict, "tau", tau)

            call get_value(child, "num_spheres", num_spheres, 10)
            call set_value(dict, "num_spheres", num_spheres)

            call get_value(child, "musb", musb, 0.0_wp)
            call set_value(dict, "musb", musb)
            call get_value(child, "muab", muab, 0.01_wp)
            call set_value(dict, "muab", muab)
            call get_value(child, "musc", musc, 0.0_wp)
            call set_value(dict, "musc", musc)
            call get_value(child, "muac", muac, 0.01_wp)
            call set_value(dict, "muac", muac)
            call get_value(child, "hgg", hgg, 0.7_wp)
            call set_value(dict, "hgg", hgg)
        else
            error stop "Need geometry table in input param file"
        end if

    end subroutine parse_geometry

    subroutine parse_output(table, dict)

        use sim_state_mod, only : state
        
        type(toml_table),  intent(INOUT) :: table, dict

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children
        integer :: i, nlen

        call get_value(table, "output", child)

        if(associated(child))then
            call get_value(child, "fluence", state%outfile, "fluence.nrrd")
            call get_value(child, "absorb", state%outfile_absorb, "absorb.nrrd")
            call get_value(child, "render", state%renderfile, "geom_render.nrrd")
            call get_value(child, "render_geom", state%render_geom, .false.)
            

            call get_value(child, "render_size", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    error stop "Need a vector of size 3 for render_size."
                end if
                do i = 1, len(children)
                    call get_value(children, i, state%render_size(i))
                end do
            else
                state%render_size = [200, 200, 200]
            end if

            call get_value(child, "overwrite", state%overwrite, .false.)
        else
            error stop "Need output table in input param file"
        end if

    end subroutine parse_output

    subroutine parse_simulation(table)

        use sim_state_mod, only : state
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child

        call get_value(table, "simulation", child)

        if(associated(child))then
            call get_value(child, "iseed", state%iseed, 123456789)
            call get_value(child, "tev", state%tev, .false.)
            call get_value(child, "absorb", state%absorb, .false.)
        else
            error stop "Need simulation table in input param file"
        end if

    end subroutine parse_simulation

    type(vector) function get_vector(child, key, default, context)

    type(toml_table),   pointer,  intent(in) :: child
    character(*),                 intent(in) :: key
    type(vector),       optional, intent(in) :: default
    type(toml_context), optional, intent(in) :: context

        type(toml_array), pointer  :: arr => null()
        real(kind=wp) :: tmp(3)
        integer :: j, origin

        call get_value(child, key, arr, origin=origin)
        if (associated(arr))then
            if(len(arr) /= 3)then
                print'(a)',context%report("Expected vector of size 3", origin, "Wrong vector size")
                stop 1    
            end if
            do j = 1, len(arr)
                call get_value(arr, j, tmp(j))
            end do
            get_vector = vector(tmp(1), tmp(2), tmp(3))
        else
            get_vector = default
        end if

    end function get_vector
end module parse_mod