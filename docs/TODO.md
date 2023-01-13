# TODOs

List of ToDo's for SignedMCRT.

## Additional Features
### Finished Features
- [x] Make CI run tests
- [x] Add Code Coverage reports
- [x] Remove spurious implicit nones
### Minor Features
- [ ] Add more saved state to photon_origin to save compute time
- [ ] Finish Circular, focus, and annulus source types
    - [ ] Circular
    - [ ] Focus
    - [x] Annulus. Partially done via control of Beta parameter. 
- [ ] Add Direction component to rest of Detectors
    - [x] **Circle**
    - [x] **Camera**
    - [ ] **Annulus**
- [ ] Add photon tracjectory history tracking
    - [ ] Fix openmp troubles
    - [ ] Fix speed issues
- [ ] Add "Scattering" kernels
    - [x] path length counter method
    - [ ] Weight method
- [ ] Add phase tracking. To be done by Masters Student?
- [ ] Make sure all optical properties are the same for a model instance (SDF)
### Major Features
- [ ] Make code work on Windows
- [ ] Automate benchmarking so we can catch performance regressions
- [ ] Improve performance of SDF intersection
    - [ ] Implement KD trees in 2D
    - [ ] Implement KD trees in 3D
- [ ] Make code serializable so that we can checkpoint simulations
     - [ ] Save input toml file
     - [ ] photons run
     - [ ] Save output data
        - [ ] Detectors
        - [ ] Fluence
        - [ ] Nscatt
- [ ] Add optics to **Camera** type

## Testing

- [x] Vec3 class
- [ ] Matrix class
 - [ ] Need Mat mult Mat test
- [x] Vec4 Class
- [ ] SDF Class
- [ ] Detector Class
    - [ ] Circle
    - [ ] Camera
    - [ ] Annulus
- [ ] Photon class
    - [ ] Uniform
    - [ ] Point
    - [ ] Pencil
    - [ ] Circular
    - [ ] Annulus
- [ ] Stack Class
- [ ] String Utils
- [ ] Utils
- [ ] I/O
- [ ] Random Numbers
- [ ] Scattering
- [ ] Fresnel reflections
- [ ] End to End tests
    - [x] Scattering Test
    - [ ] Others
- [ ] Implment a fuzzer to test input space
- [ ] Add a input schema to validate input toml
    - [ ] Add Even better TOML?
    - [ ] Learn from https://json-schema.org/learn/getting-started-step-by-step.html

## Bugs
- [ ] Fix CI so that build on Macos run, and builds using Intel run.
    - [ ] Macos
    - [ ] Intel
- [ ] Can't operate trackHistory in parallel
    - [ ] Make each thread write to tmp file and finish method collate results
- [x] Always runs history%finish leading to error messages
- [x] photon history filename not set when no detectors are used. Causes crash.
    - [x] Also causes crash on write out if no detectors used.
