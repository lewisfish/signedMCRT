module iarray

!!    The iarray module contains the variables that record the fluence. These are 3D arrays, with roughly the same dimensions as the cart_grid type.
!!    Jmean is the *local* fluence. JmeanGLOBAL is the *global* fluence grid. The global version is the one that is written to disk at the simulations end.

    use constants, only : sp

    implicit none

    complex(kind=sp), allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
    real(kind=sp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    real(kind=sp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:)
end module iarray