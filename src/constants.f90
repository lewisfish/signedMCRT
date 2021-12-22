module constants
!
! Module containing constants:
!         PI, TWOPI
!         and the vars for the filepaths.
!         resdir holds the path to the directory that holds the parameter and other associated input files
!         cwd is the current working directory
!         homedir is the root directory of this code
!         fileplace is the data folder directory
    implicit none

    real,    parameter :: PI=4.*atan(1.), TWOPI=2.*PI
    character(len=255) :: cwd, homedir, fileplace, resdir

end module constants