module sim_state_mod

    use gridMod,   only : cart_grid

    implicit none
    
    type :: settings_t
        integer :: nphotons, iseed
        character(len=:), allocatable :: experiment, outfile, renderfile, source
        type(cart_grid) :: grid
    end type settings_t

    type(settings_t) :: state

    private
    public :: settings_t, state    

end module sim_state_mod