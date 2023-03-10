module iarray
!
!  Contains all array var names.
!
    use constants, only : wp

    implicit none

    !real(kind=wp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    complex(kind=wp), allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
    real(kind=wp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    real(kind=wp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:)
end module iarray