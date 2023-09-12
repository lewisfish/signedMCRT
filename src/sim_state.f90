module sim_state_mod

    !! This module defines the setting_t type which holds simulation metadata:

    !! - nphotons. Number of photons to run
    !! - iseed. initial seed
    !! - render_size. Size of voxel grid to render SDFs to 
    !! - experiment. Name of experiment/simulation
    !! - outfile. Name of fluence output file
    !! - renderfile. Name of voxel render file
    !! - source. Light source used
    !! - historyFilename. Name of photon history file
    !! - grid. Cart_grid type
    !! - render_geom. Boolean to indicate wether to render SDF to voxels or not.
    !! - tev. Boolean to indicate wether to use TEV as debug viewer.


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