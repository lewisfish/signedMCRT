module parse_detectorsMod

    use tomlf
    use tomlf_error, only : make_error
    use constants, only : wp
    use vector_class
    use parseHelpers

    implicit none

    private
    public :: parse_detectors

contains
    
    subroutine parse_detectors(table, dects, context, error)
        !! parse the detectors

        use detectors,     only : dect_array, circle_dect, annulus_dect, camera
        use sim_state_mod, only : state

        !> Input Toml table
        type(toml_table),               intent(inout) :: table
        !> Detector array to be filled.
        type(dect_array), allocatable,  intent(out)   :: dects(:)
        !> Context handle for error reporting.
        type(toml_context),             intent(in)    :: context
        !> Error message
        type(toml_error),  allocatable, intent(out)   :: error

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
                call make_error(error, &
                    context%report("Invalid detector type. Valid types are [circle, annulus, camera]", &
                    origin, "expected valid detector type"), -1)
                return
            case("circle")
                c_counter = c_counter + 1
            case("annulus")
                a_counter = a_counter + 1
            case("camera")
                cam_counter = cam_counter + 1
            end select
        end do

        if(c_counter > 0)then
            if(allocated(dect_c))deallocate(dect_c)
            allocate(dect_c(c_counter))
        end if
        if(a_counter > 0)then
            if(allocated(dect_a))deallocate(dect_a)
            allocate(dect_a(a_counter))
        end if
        if(cam_counter > 0)then
            if(allocated(dect_cam))deallocate(dect_cam)
            allocate(dect_cam(cam_counter))
        end if
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
                call handle_circle_dect(child, dect_c, c_counter, context, error)
                if(allocated(error))return
            case("annulus")
                call handle_annulus_dect(child, dect_a, a_counter, context, error)
                if(allocated(error))return
            case("camera")
                call handle_camera(child, dect_cam, cam_counter, context, error)
                if(allocated(error))return
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

    subroutine handle_camera(child, dects, counts, context, error)
        !! Read in Camera settings and initalise variable
        use detectors,     only : camera
        use sim_state_mod, only : state

        !> Detector table
        type(toml_table), pointer,     intent(in)    :: child
        !> Array of cameras
        type(camera),                  intent(inout) :: dects(:)
        !> Number of cameras to create
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error

        integer       :: layer, nbins
        real(kind=wp) :: maxval
        type(vector)  :: p1, p2, p3
        logical       :: trackHistory

        p1 = get_vector(child, "p1", default=vector(-1.0, -1.0, -1.0), context=context, error=error)
        p2 = get_vector(child, "p2", default=vector(2.0, 0.0, 0.0), context=context, error=error)
        p3 = get_vector(child, "p3", default=vector(0.0, 2.0, 0.0), context=context, error=error)

        call get_value(child, "layer", layer, 1)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
        dects(counts) = camera(p1, p2, p3, layer, nbins, maxval, trackHistory)
        counts = counts + 1

    end subroutine handle_camera

    subroutine handle_circle_dect(child, dects, counts, context, error)
        !! Read in Circle_detector settings and initalise variable
        use detectors,     only : circle_dect
        use sim_state_mod, only : state

        !> Detector table
        type(toml_table), pointer,     intent(in)    :: child
        !> Array ofcircle_dects
        type(circle_dect),             intent(inout) :: dects(:)
        !> Number of circle_dects to create
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error

        integer       :: layer, nbins
        real(kind=wp) :: maxval, radius
        type(vector)  :: pos, dir
        logical       :: trackHistory

        pos = get_vector(child, "position", context=context, error=error)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context, error=error)
        dir = dir%magnitude()
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
        dects(counts) = circle_dect(pos, dir, layer, radius, nbins, maxval, trackHistory)
        counts = counts + 1

    end subroutine handle_circle_dect

    subroutine handle_annulus_dect(child, dects, counts, context, error)
        !! Read in Annulus_detector settings and initalise variable
        
        use detectors,     only : annulus_dect
        use sim_state_mod, only : state
        use utils,         only : str

        !> Detector Table
        type(toml_table), pointer,     intent(in)    :: child
        !> Array of annulus_dects
        type(annulus_dect),            intent(inout) :: dects(:)
        !> Number of anulluar dects to create
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error

        integer       :: layer, nbins, origin
        real(kind=wp) :: maxval, radius1, radius2
        type(vector)  :: pos, dir
        logical       :: trackHistory

        pos = get_vector(child, "position", context=context, error=error)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context, error=error)
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius1)
        call get_value(child, "radius2", radius2, origin=origin)
        
        if(radius2 <= radius1)then
            call make_error(error,&
            context%report("Radii are invalid", origin, &
            "Expected radius2 ("//str(radius2,6)//") > radius 1 ("//str(radius1,6)//")"), -1)
            return
        end if
        
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
        dects(counts) = annulus_dect(pos, dir, layer, radius1, radius2, nbins, maxval, trackHistory)
        counts = counts + 1
    end subroutine handle_annulus_dect
end module parse_detectorsMod