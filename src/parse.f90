module parse_mod

    use tomlf
    use constants, only : wp
    use vector_class

    implicit none

    private
    public :: parse_params

    contains

    subroutine parse_params(filename, packet, dects, dict)
    
        use photonmod
        use detector_mod, only : dect_array

        implicit none
    
        character(*),      intent(IN)    :: filename
        type(toml_table),  intent(INOUT) :: dict
        type(photon),      intent(OUT)   :: packet
        type(dect_array), allocatable, intent(out) :: dects(:)

        type(toml_table), allocatable :: table
        type(toml_context) :: context
        type(toml_error), allocatable :: error

        call toml_load(table, trim(filename), context=context, error=error)
        if(allocated(error))then
            print'(a)',error%message
            stop 1
        end if

        call parse_source(table, packet, dict, context)
        call parse_grid(table, dict)
        call parse_geometry(table, dict)
        call parse_detectors(table, dects, context)
        call parse_output(table, dict)
        call parse_simulation(table)

    end subroutine parse_params
    
    subroutine parse_detectors(table, dects, context)

        use detector_mod

        type(toml_table) :: table
        type(dect_array), allocatable :: dects(:)
        type(toml_context), intent(in) :: context

        type(toml_array), pointer :: array
        type(toml_table), pointer :: child
        character(len=:), allocatable :: dect_type
        type(circle_dect), target, save, allocatable :: dect_c(:)
        type(annulus_dect), target, save, allocatable :: dect_a(:)
        type(camera), target, save, allocatable :: dect_cam(:)
        integer :: i, c_counter, a_counter, cam_counter, j, origin, k

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
                print'(a)',context%report("Invalid detector type. Valid types are [circle, annulus]", origin, "expected valid detector type")
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
        do i = 1, len(array)
            call get_value(array, i, child)
            call get_value(child, "type", dect_type)
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

    end subroutine parse_detectors

    subroutine handle_camera(child, dects, counts, context)

        use detector_mod

        type(toml_table), pointer, intent(in)    :: child
        type(camera),              intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins
        real(kind=wp) :: maxval
        type(vector)  :: p1, p2, p3

        p1 = get_vector(child, "p1", context)
        p2 = get_vector(child, "p2", context)
        p3 = get_vector(child, "p3", context)
 
        call get_value(child, "layer", layer, 1)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        dects(counts) = camera(p1, p2, p3, layer, nbins, maxval)
        counts = counts + 1

    end subroutine handle_camera

    subroutine handle_circle_dect(child, dects, counts, context)

        use detector_mod

        type(toml_table), pointer, intent(in)    :: child
        type(circle_dect),         intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins
        real(kind=wp) :: maxval, radius
        type(vector)  :: pos

        pos = get_vector(child, "position", context)
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        dects(counts) = circle_dect(pos, layer, radius, nbins, maxval)
        counts = counts + 1

    end subroutine handle_circle_dect

    subroutine handle_annulus_dect(child, dects, counts, context)

        use detector_mod

        type(toml_table), pointer, intent(in)    :: child
        type(annulus_dect),        intent(inout) :: dects(:)
        integer,                   intent(inout) :: counts
        type(toml_context),        intent(in) :: context

        integer       :: layer, nbins, origin
        real(kind=wp) :: maxval, radius1, radius2
        type(vector)  :: pos

        pos = get_vector(child, "position", context)
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius1)
        call get_value(child, "radius2", radius2, origin=origin)
        if(radius2 <= radius1)then
            print'(a)',context%report("Radii are invalid", origin, "Expected radius2 > radius 1")
            stop 1
        end if
            call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        dects(counts) = annulus_dect(pos, layer, radius1, radius2, nbins, maxval)
        counts = counts + 1
    end subroutine handle_annulus_dect

    subroutine parse_source(table, packet, dict, context)

        use sim_state_mod, only : state
        use photonmod

        implicit none
        
        type(toml_table),  intent(INOUT) :: table, dict
        type(photon),      intent(OUT)   :: packet
        type(toml_context) :: context

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        type(vector) :: poss, dirr
        real(kind=wp) :: dir(3), pos(3), corners(3, 3), radius
        integer :: i, nlen, origin
        character(len=1) :: axis(3)
        character(len=:), allocatable :: direction

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
                    print'(a)',context%report("Need a vector of size 3 for position", origin, "expected vector of size 3")
                    stop 1
                end if
                do i = 1, len(children)
                    call get_value(children, i, pos(i))
                    call dict%set(key("pos%"//axis(i)), value=pos(i))
                end do
            else
                if(state%source == "point")then
                    pos = [0._wp, 0._wp, 0._wp]
                    call dict%set(key("pos%x"), value=pos(1))
                    call dict%set(key("pos%y"), value=pos(2))
                    call dict%set(key("pos%z"), value=pos(3))
                end if
            end if

            children => null()
            
            call get_value(child, "direction", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    print'(a)',context%report("Need a vector of size 3 for direction", origin, "expected vector of size 3")
                    stop 1
                end if
                if(state%source == "circular")then
                    print'(a)',context%report("Direction not yet fully implmented for source type Circular. Results may not be accurate!", origin, level=toml_level%warning)
                end if
                do i = 1, len(children)
                    call get_value(children, i, dir(i))
                end do
                dirr%x = dir(1)
                dirr%y = dir(2)
                dirr%z = dir(3)
            else
                if(state%source == "uniform" .or. state%source == "circular")then
                    print'(a)',context%report("Uniform source needs vector direction", origin, "expected vector of size 3")
                    stop 1
                end if
                call get_value(child, "direction", direction, "-z")
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
                    error stop "fucked it!"
                end select
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
        else
            print'(a)',context%report("Simulation needs Source table", origin, "Missing source table")
            stop 1
        end if

        call set_photon(poss, dirr)
        packet = photon(state%source)

    end subroutine parse_source

    subroutine parse_grid(table, dict)

        use sim_state_mod, only : state
        use gridMod,       only : init_grid 

        implicit none
        
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
        
        implicit none
        
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
        logical :: overwrite
        integer :: i, nlen

        call get_value(table, "output", child)

        if(associated(child))then
            call get_value(child, "fluence", state%outfile, "fluence.nrrd")
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

            call get_value(child, "overwrite", overwrite, .false.)
            call set_value(dict, "overwrite", overwrite)
        else
            error stop "Need output table in input param file"
        end if

    end subroutine parse_output

    subroutine parse_simulation(table)

        use sim_state_mod, only : state

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child

        call get_value(table, "simulation", child)

        if(associated(child))then
            call get_value(child, "seed", state%iseed, 123456789)
            call get_value(child, "tev", state%tev, .false.)

        else
            error stop "Need simulation table in input param file"
        end if

    end subroutine parse_simulation

    type(vector) function get_vector(child, key, context)

        type(toml_context),        intent(in) :: context
        type(toml_table), pointer, intent(in) :: child
        character(*),              intent(in) :: key

        type(toml_array), pointer  :: arr
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
        end if

    end function get_vector

end module parse_mod
