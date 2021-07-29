MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.
!

    implicit none

    real,    parameter :: PI=4.*atan(1.), TWOPI=2.*PI
    character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants