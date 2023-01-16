module constants

    use iso_fortran_env, only : real64

! Module containing constants:
!         PI, TWOPI
!         and the vars for the filepaths.
!         resdir holds the path to the directory that holds the parameter and other associated input files
!         homedir is the root directory of this code
!         fileplace is the data folder directory
    implicit none

    integer,          parameter :: wp = real64 !can change this to other precision, not tested for lower or higher precisions.
    real(kind=wp),    parameter :: PI=4._wp*atan(1._wp), TWOPI=2._wp*PI
    character(len=255)          :: homedir, fileplace, resdir

end module constants