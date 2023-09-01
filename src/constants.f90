module constants
    !! Module contains the various constants used in the simulation.
    !! Also defines the working precision and variables for folder operations.
    use iso_fortran_env, only : real64, real32

    implicit none

    !> current working precision
    integer,           parameter :: wp = real64 !can change this to other precision, not tested for lower or higher precisions.
    !> single precision variable.
    integer,           parameter :: sp = real32
    !> \[\pi\]
    real(kind=wp),     parameter :: PI=4._wp*atan(1._wp)
    !> \[2 \pi\]
    real(kind=wp),     parameter :: TWOPI=2._wp*PI
    !> Weight threshold for roulette
    real(kind=wp),     parameter :: THRESHOLD = 0.01_wp
    !> Proportion of packet that survive roulette
    real(kind=wp),     parameter :: CHANCE = 0.1_wp
    !> root directory
    character(len=255)          :: homedir
    !> place where output files are saved
    character(len=255)          :: fileplace
    !> directory to input files
    character(len=255)          :: resdir

end module constants