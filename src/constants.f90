module constants

    use iso_fortran_env, only : real32, real64

! Module containing constants:
!         PI, TWOPI
!         and the vars for the filepaths.
!         resdir holds the path to the directory that holds the parameter and other associated input files
!         cwd is the current working directory
!         homedir is the root directory of this code
!         fileplace is the data folder directory
    implicit none

    integer,          parameter :: wp = real64 !can change this to real32 if need be
    real(kind=wp),    parameter :: PI=4._wp*atan(1._wp), TWOPI=2._wp*PI
    character(len=255)          :: cwd, homedir, fileplace, resdir


end module constants