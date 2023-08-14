module sim_state_mod

    use gridMod,   only : cart_grid

    implicit none
    
    type :: settings_t
        integer :: nphotons, iseed
        integer :: render_size(3)
        character(len=:), allocatable :: experiment, outfile, renderfile, source, historyFilename, outfile_absorb
        type(cart_grid) :: grid
        logical :: render_geom, tev, overwrite, trackHistory, absorb
    end type settings_t

    type(settings_t) :: state

    private
    public :: settings_t, state
end module sim_state_mod