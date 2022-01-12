module parse_mod

    use tomlf
    use constants, only : wp

    implicit none

    private
    public :: parse_params

    contains

    subroutine parse_params(filename, packet, dict)
    
        use fhash, only : fhash_tbl_t
        use photonmod

        implicit none
    
        character(*),      intent(IN)    :: filename
        type(fhash_tbl_t), intent(INOUT) :: dict
        type(photon),      intent(OUT)   :: packet

        type(toml_table), allocatable :: table

        integer :: u

        open(newunit=u, file=trim(filename))
        call toml_parse(table, u)
        
        call parse_source(table, packet, dict)
        call parse_grid(table)
        call parse_geometry(table, dict)
        call parse_output(table)
        call parse_simulation(table)

    end subroutine parse_params
    
    subroutine parse_source(table, packet, dict)

        use fhash, only : fhash_tbl_t, key=>fhash_key
        use sim_state_mod, only : state
        use photonmod

        implicit none
        
        type(toml_table),  intent(INOUT) :: table
        type(fhash_tbl_t), intent(INOUT) :: dict
        type(photon),      intent(OUT)   :: packet

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        real(kind=wp) :: dir(3), pos(3)
        integer :: i, nlen
        character(len=1) :: axis(3)
        character(len=:), allocatable :: direction

        axis = ["x", "y", "z"]

        call get_value(table, "source", child)

        if(associated(child))then
            call get_value(child, "name", state%source, "point")
            call get_value(child, "nphotons", state%nphotons, 1000000)

            call get_value(child, "position", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    error stop "Need a vector of size 3 for position."
                end if
                do i = 1, len(children)
                    call get_value(children, i, pos(i))
                    call dict%set(key("pos"//axis(i)), value=pos(i))
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
            
            call get_value(child, "direction", children, requested=.false.)
            if(associated(children))then
                if(state%source == "uniform")error stop "Source uniform cant have vector direction!"
                nlen = len(children)
                if(nlen < 3)then
                    error stop "Need a vector of size 3 for direction."
                end if
                do i = 1, len(children)
                    call get_value(children, i, dir(i))
                    call dict%set(key("dir%"//axis(i)), value=dir(i))
                end do
            else
                call get_value(child, "direction", direction, "-z")
                call dict%set(key("dir"), value=direction)
            end if
        else
            error stop "Need source table in input param file"
        end if

        packet = photon(state%source)

    end subroutine parse_source

    subroutine parse_grid(table)

        use sim_state_mod, only : state
        use gridMod,       only : init_grid 

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child
        integer                   :: nxg, nyg, nzg
        real(kind=wp)             :: xmax, ymax, zmax

        call get_value(table, "grid", child)

        if(associated(child))then
            call get_value(child, "nxg", nxg, 200)
            call get_value(child, "nyg", nyg, 200)
            call get_value(child, "nzg", nzg, 200)
            call get_value(child, "xmax", xmax, 1.0_wp)
            call get_value(child, "ymax", ymax, 1.0_wp)
            call get_value(child, "zmax", zmax, 1.0_wp)
        else
            error stop "Need grid table in input param file"
        end if

        state%grid = init_grid(nxg, nyg, nzg, xmax, ymax, zmax)

    end subroutine parse_grid

    subroutine parse_geometry(table, dict)

        use fhash,         only : fhash_tbl_t, key=>fhash_key
        use sim_state_mod, only : state
        
        implicit none
        
        type(toml_table),  intent(INOUT) :: table
        type(fhash_tbl_t), intent(INOUT)    :: dict
        
        type(toml_table), pointer :: child
        real(kind=wp)             :: tau

        call get_value(table, "geometry", child)

        if(associated(child))then
            call get_value(child, "geom_name", state%experiment, "sphere")
            call get_value(child, "tau", tau, 10._wp)
            call dict%set(key("tau"), value=tau)
        else
            error stop "Need geometry table in input param file"
        end if

    end subroutine parse_geometry

    subroutine parse_output(table)

        use sim_state_mod, only : state

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child

        call get_value(table, "output", child)

        if(associated(child))then
            call get_value(child, "fluence", state%outfile, "fluence.nrrd")
            call get_value(child, "render", state%renderfile, "geom_render.nrrd")
            call get_value(child, "render_geom", state%render_geom, .false.)
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
        else
            error stop "Need simulation table in input param file"
        end if

    end subroutine parse_simulation
end module parse_mod
