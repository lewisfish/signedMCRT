module sim_state_mod
!! This module defines the setting_t type which holds simulation metadata:

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
        !> Number of photons to run
        integer :: nphotons
        !> initial seed for random number generator
        integer :: iseed
        !> Size of the voxel grid to render SDFs to
        integer :: render_size(3)
        !> Name of experiment/simulation
        character(len=:), allocatable :: experiment
        !> Name of fluence output file
        character(len=:), allocatable :: outfile
        !> Name of voxel render file
        character(len=:), allocatable :: renderfile
        !> Light source used
        character(len=:), allocatable :: source
        !> Name of photon history file
        character(len=:), allocatable :: historyFilename
        !> Name of absoprtion output file
        character(len=:), allocatable :: outfile_absorb
        !> Cart_grid type
        type(cart_grid) :: grid
        !> Boolean to indicate whether to render SDF to voxels or not.
        logical :: render_geom
        !> Boolean to indicate whether to use TEV as debug viewer.
        logical :: tev
        !> Boolean to indicate whether to use overwrite datafiles if they have the same name.
        logical :: overwrite
        !> Boolean to indicate whether to store history of photons positions
        logical :: trackHistory
        !> Boolean to indicate whether to store absoption data.
        logical :: absorb
    end type settings_t

    !> global var that stores simulation state
    type(settings_t) :: state

    private
    public :: settings_t, state
end module sim_state_mod