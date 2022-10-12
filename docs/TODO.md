# TODOs

List of ToDo's for SignedMCRT.

## Additional Features

- [x] Make CI run tests
- [x] Add Code Coverage reports
- [ ] Add Direction component to rest of Detectors
    - [x] **Circle**
    - [x] **Camera**
    - [ ] **Annulus**
- [ ] Add optics to **Camera** type
- [ ] Add "Scattering" kernels
    - [ ] path length counter method
    - [ ] Weight method
- [ ] Make sure all optical properties are the same for a model instance (SDF)
- [ ] Remove spurious implicit nones
- [ ] Improve performance of SDF intersection
    - [ ] Implement KD trees in 2D
    - [ ] Implement KD trees in 3D
- [ ] Make code serializable so that we can checkpoint simulations


## Testing

- [x] Vec3 class
- [x] Matrix class
 - [ ] Need Mat mult Mat test
- [ ] Vec4 Class
- [ ] SDF Class
- [ ] Detector Class
- [ ] Photon class
- [ ] Stack Class
- [ ] String Utils
- [ ] Utils
- [ ] I/O
- [ ] Random Numbers
- [ ] Scattering
- [ ] Fresnel reflections
- [ ] End to End tests
    - [ ] Scattering Test
    - [ ] Others
- [ ] Implment a fuzzer to test input space
- [ ] Add a input schema to validate input toml
    - [ ] Add Even better TOML?
    - [ ] Learn from https://json-schema.org/learn/getting-started-step-by-step.html

## Bugs

- [ ] Can't operate trackHistory in parallel
    - [ ] Make each thread write to tmp file and finish method collate results
- [x] Always runs history%finish leading to error messages
- [x] photon history filename not set when no detectors are used. Causes crash.
    - [x] Also causes crash on write out if no detectors used.