module iarray
!
!  Contains all array var names.
!
    use constants, only : wp

    implicit none

    real(kind=wp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:) ! Fluence
    real(kind=wp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:) ! Absorbed energy
end module iarray