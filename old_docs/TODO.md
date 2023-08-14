---
title: todos
---
# TODOs

List of ToDo's for SignedMCRT.

## Additional Features
### Finished Features
- [x] Make CI run tests
- [x] Add Code Coverage reports
- [x] Remove spurious implicit nones
- [x] Make sure all optical properties are the same for a model instance (SDF)
- [x] Add "Scattering" kernels
    - [x] path length counter method
    - [x] Weight method
- [x] Add documentation on piecewise
    - [x] Constant
    - [x] 1D
    - [x] 2D
- [x] Finish new SDF API
    - [x] add all SDFs
    - [x] add adjustment functions (twist, union etc)
    - [x] propagate to subs.f90

### Minor Features
- [ ] Add more saved state to photon_origin to save compute time
- [ ] Finish Circular, focus, and annulus source types
    - [ ] Circular
    - [ ] Focus
    - [x] Annulus. Partially done via control of Beta parameter. 
- [ ] Add Direction component to rest of Detectors
    - [x] Circle
    - [x] Camera
    - [ ] Annulus
- [ ] Add photon trajectory history tracking
    - [ ] Add to each detector separately 
    - [ ] Fix openMP troubles
    - [ ] Fix speed issues
### Major Features
- [ ] Make code work on Windows
- [ ] Automate benchmarking so we can catch performance regressions
- [ ] Add voxel geometry
- [ ] Add mesh geometry
- [ ] Improve performance of SDF intersection
    - [ ] Implement KD trees in 2D
    - [ ] Implement KD trees in 3D
- [ ] Make code serializable so that we can checkpoint simulations
     - [ ] Save input toml file
     - [ ] photons run
     - [ ] Save output data
        - [ ] Detectors
        - [ ] Fluence
        - [ ] Absorb
        - [ ] NScatt
- [ ] Add optics to Camera type
- [x] Add phase tracking (https://github.com/lewisfish/signedMCRT/pull/2).
    - [ ] Add phase screen detector to camera
    - [ ] Add refractive index accounting
- [ ] Compress output data (https://github.com/aras-p/float_compr_tester/blob/main/src/compression_helpers.cpp)
- [ ] Add more error handling for spectrums in parse.f90
- [ ] Add optical property type, to allow for multi-spectral input.
    - [x] base optical property type
        - [ ] function defined
        - [x] Tabulated
    - [x] propagate to SDFs
    - [x] propagate to subs.f90
    - [ ] Document optical properties
    - [x] Change API to match that of SDFs, i.e easier to use
- [ ] Add MPI + openMP mode (e.g. run openMP on N nodes with minimal communication)

## Testing

- [x] Vec3 class
- [x] Matrix class
- [x] Vec4 Class
- [ ] SDF Class
- [ ] Detector Class
    - [ ] Circle
    - [ ] Camera
    - [ ] Annulus
- [ ] Photon class
    - [x] Uniform
    - [x] Point
    - [x] Pencil
    - [ ] Circular
    - [ ] Annulus
    - [ ] Scattering
     - [ ] Isotropic
     - [ ] Henyey-Greenstein
     - [ ] Importance sampling biased scattering
- [ ] Photon movement code
- [ ] History Stack Class
- [ ] I/O
- [ ] Random Numbers
- [x] Fresnel reflections
    - [x] Simple reflect
    - [x] Simple refract
    - [x] Complex reflect
    - [x] Complex refract
- [ ] End to End tests
    - [x] Scattering Test
    - [ ] Others
- [ ] test phase
    - [ ] double slit
    - [ ] square aperture
    - [ ] gaussian beam
    - [ ] bessel beam

## Bugs
- [ ] Fix CI so that build on Macos run, and builds using Intel run.
    - [ ] Macos
    - [ ] Intel
- [ ] Can't operate trackHistory in parallel
    - [ ] Make each thread write to tmp file and finish method collate results
- [x] Added default array option to get_vector in parse.f90