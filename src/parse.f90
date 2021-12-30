module parse_mod

    use tomlf
    use constants, only : wp

    implicit none

    private
    public :: parse_params

    contains

    subroutine parse_params(filename)
    
        implicit none
    
        character(*), intent(IN) :: filename
    
        type(toml_table), allocatable :: table

        integer :: u

        open(newunit=u, file=trim(filename))
        call toml_parse(table, u)
        
        call parse_source(table)
        call parse_grid(table)
        call parse_geomtry(table)
        call parse_output(table)
        call parse_simulation(table)

    end subroutine parse_params
    
    subroutine parse_source(table)

        use sim_state_mod, only : state

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child
        real(kind=wp) :: dir(3)

        call get_value(table, "source", child)

        if(associated(child))then
            call get_value(child, "name", state%source, "point")
            call get_value(child, "nphotons", state%nphotons, 1000000)
            call get_value(child, "direction", dir, [0._wp, 0._wp, -1._wp])
        else
            error stop "Need source table in input param file"
        end if

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

    subroutine parse_geomtry(table)

        use sim_state_mod, only : state

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child

        call get_value(table, "geomtry", child)

        if(associated(child))then
            call get_value(child, "geom_name", state%experiment, "sphere")
        else
            error stop "Need geomtry table in input param file"
        end if

    end subroutine parse_geomtry

    subroutine parse_output(table)

        use sim_state_mod, only : state

        implicit none
        
        type(toml_table), intent(INOUT) :: table

        type(toml_table), pointer :: child

        call get_value(table, "output", child)

        if(associated(child))then
            call get_value(child, "fluence", state%outfile, "fluence.nrrd")
            call get_value(child, "render", state%renderfile, "geom_render.nrrd")
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
