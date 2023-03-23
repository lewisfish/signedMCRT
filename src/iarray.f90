module iarray
!
!  Contains all array var names.
!
    use constants, only : sp

    implicit none

    complex(kind=sp), allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
    real(kind=sp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    real(kind=sp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:)
end module iarray