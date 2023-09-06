module iarray
!!  Contains all array var names.

    use constants, only : sp

    implicit none
    !> phase data array
    complex(kind=sp), allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
    !> fluence data array
    real(kind=sp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    !> absorption data array
    real(kind=sp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:)
end module iarray