module iarray
!
!  Contains all array var names.
!
    use constants, only : wp

    implicit none

    real(kind=wp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
end module iarray