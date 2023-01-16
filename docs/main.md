# Documentation

This document is the incomplete documentation of **signedMCRT**.

## Build system

To build signedMCRT, the only current method is using [FPM](https://fpm.fortran-lang.org/en/index.html).
FPM can be easily installed on any platform, and is simple to use to pull all dependencies, and build and compile signedMCRT.
We also provide several commands via FPM response file ((found here)[fpm.rsp]), to enable the use of OpenMP, other compliers, and various debug modes.

## Dependencies

Below is the current list of dependencies:

* [test drive](https://github.com/fortran-lang/test-drive)
* [Fortran TEV Bindings](https://github.com/lewisfish/fortran_tev_bindings)
* [Fortran Utilities](https://github.com/lewisfish/fortran_utilities)

Test drive is used to run all tests.
Fortran [TEV](https://github.com/Tom94/tev/) Bindings is used to interface with TEV, to show live slices of fluences as the simulation is run, which is handy for debugging purposes.
Finally, Fortran Utilities is my personal collection of useful fortran utilities such as mathematical functions, or progress bars.

## Config file
signedMCRT uses TOML as it's configuration file format.
Documentation of the input file format can be found in [here](config.md)


## Code Structure

See below for a brief description of each source code file.
See (this)[relations.pdf] for rat's nest diagram of the relationships and dependencies between each source code file. 
Blue is a module with functions.
Grey is a external dependency.
Red is a module with only data components.

### constants.f90

This module contains mathematical constants and strings that contain the various directories used by the program.
Math constants:
- PI
- 2 PI
- wp (working precision of the whole program). Default is double precision (64bit floats)
Directories:
- homedir. Root directory of this code
- fileplace. data folder directory
- resdir. holds the path to the directory that holds the parameter and other associated input files

Source code can be found [here](/src/constants.f90)

### detectorMod.f90

Currently in development.

Source code can be found [here](/src/detectorMod.f90)

### geometeryMod.f90

Defines a set of functions for intersecting a line and a surface.

- Circle
- Plane
- Cone
- Cylinder
- Ellipse
- Sphere

Source code can be found [here](/src/geometeryMod.f90)

### grid.f90

This module defines the cartesian grid type (cart_grid) and associated routines.

The cart_grid type contains information related to the grid used to record the fluence. This includes the number of voxels in each cardinal direction (nxg, nyg, nzg), the **half** size of the grid in each direction (xmax, ymax, zmax), and the locations of the voxels walls in each direction (xface, yface, zface).
The type-bound function get_voxel takes a position (vector) and returns the voxel the position falls in. 

Init_grid initialises a cart_grid instance.

Source code can be found [here](/src/grid.f90)

### historyStack.f90

Currently in development.

Source code can be found [here](/src/historyStack.f90)

### iarray.f90

The iarray module contains the variables that record the fluence. These are 3D arrays, with roughly the same dimensions as the cart_grid type.
Jmean is the *local* fluence. JmeanGLOBAL is the *global* fluence grid. The global version is the one that is written to disk at the simulations end.

Source code can be found [here](/src/iarray.f90)

### inttau2.f90

inttau2 is the heart of the MCRT simulation. It moves the photons though the simulated media.
tauint2 is the only public function here and is the main function that moves the photon.
Changes should only be made here if bugs are discovered or new methods of tracking photons (i.e phase tracking) or moving photons (i.e new geometry method) is needed.

Source code can be found [here](/src/inttau2.f90)

### kernelsMod.f90

Contains the main program and scattering loop. Calls all other routine to setup, run and break down the simulation.

Source code can be found [here](/src/kernelsMod.f90)

### mat_class.f90

Matrix class module. Defines a matrix type (4x4 matrix) and associated operations on matrices and other types.

Source code can be found [here](/src/mat_class.f90)

### photon.f90

This source file contains the photon type, all the photon launch routines for different light sources, and the scattering code.

Below are the current types of light sources available. Check [here](config.md) for parameters needed for each light source.
- uniform
- pencil
- annulus
- focus
- point
- circular

Source code can be found [here](/src/photon.f90)

### parse.f90

This file contains all the logic for parsing the toml input file using the toml-f library.
parse_params is the only public function in this module. This function returns a dictionary (toml_table) of parameters needed to setup the simulation and an array of detectors if any are defined.

Any additions to the toml input file should be added in here.

Source code can be found [here](/src/parse.f90)

### random_mod.f90

This module defines a set of functions that return random numbers in different distributions.

- ran2. Returns a single float uniformly in the range [0, 1)
- ranu. Return a single float uniformly in the range [a, b)
- randint. Returns a single integer uniformly in the range [a, b)
- rang. Returns a single float from a Gaussian distribution with mean *avg* and std *sigma*.
- init_rng. Seeds the internal random number generator with a reproducible seed.

Source code can be found [here](/src/random_mod.f90)

### sdfsMod.f90

This module defines the signed distance function (SDF) abstract type and all types that inherit from it.
The SDF abstract type defines the optical properties of an SDF (mus, mua, kappa, albedo, hgg, g2,and n), as well as a transform (4x4 matrix), and the layer ID code of the SDF. The SDF abstract type also provides an abstract interface (evaluate) which each inheriting function must implement. This evaluate function is the heart of the SDF implementation. Each individual evaluate is the direct implementation of that SDF, e.g. that function defines the mathematical SDF.

- cylinder
- sphere
- box
- torus
- cone
- triprism (triangular prism)
- capsule
- plane
- moon (not tested)
- neural (requires the use of separate python scripts to generate SDF. Current model is the Stanford bunny.)
- segment
- egg

This module also defines tow meta-container for dealing with multiple SDFs at the same time:

- model. This joins multiple SDFs into one SDF type.
- container. This is essentially an array type for storing multiple SDFs at once.

This module also defines transforms that can be applied to each SDF:

- Union
- Intersection
- Subtraction
- SmoothUnion
- Rotate_{x,y,z}
- Translate
- RotationAlign (not tested)
- RotMat (not tested)
- Displacement
- Bend
- Twist
- Elongate
- Repeat
- Extrude
- Revolution
- Onion

Source code can be found [here](/src/sdfsMod.f90)

### sim_state.f90

This module defines the setting_t type which holds simulation metadata:

- nphotons. Number of photons to run
- iseed. initial seed
- render_size. Size of voxel grid to render SDFs to 
- experiment. Name of experiment/simulation
- outfile. Name of fluence output file
- renderfile. Name of voxel render file
- source. Light source used
- historyFilename. Name of photon history file
- grid. Cart_grid type
- render_geom. Boolean to indicate wether to render SDF to voxels or not.
- tev. Boolean to indicate wether to use TEV as debug viewer.

Source code can be found [here](/src/sim_state.f90)

### subs.f90

This file sets up some simulations variables and initialises the geometry for the simulation.
Add any new geometry setups here.

Source code can be found [here](/src/subs.f90)

### surfaces.f90

Contains the routines that handle reflection, and refraction via the Fresnel equations.

Source code can be found [here](/src/surfaces.f90)

### vec4_class.f90

Vector4 class module. Defines a vector4 type (x, y, z, p) and associated operations on vectors and other types.

Source code can be found [here](/src/vec4_class.f90)

### vector_class.f90

Vector class module. Defines a vector type (x, y, z) and associated operations on vectors and other types.

Source code can be found [here](/src/vector_class.f90)

### writer.f90

This module defines all functions that write simulation data to the disk or pre-process data before writing.

- normalise_fluence. Normalises fluence by number of photons run and size of each voxel. **!Does not normalise by power!**
- write_fluence. Write out fluence in either raw or nrrd format. Default is nrrd.
- write_detected_photons. Write out photons detected by detectors.

Changes should only be made here if there is a bug or new data types need to be written to disk (phase information) or new file format is needed.

Source code can be found [here](/src/writer.f90)

## Monte Carlo Radiation Transfer (MCRT) method

Please see my [thesis](main.pdf) for an overview of the MCRT method

## Citation

SignedMCRT has so far been used in 2 papers:

+ MESHLESS MONTE CARLO RADIATION TRANSFER METHOD FOR CURVED GEOMETRIES USING SIGNED DISTANCE FUNCTIONS
L. McMillan, G. D. Bruce, K. Dholakia, [J. Biomed. Opt. 27(8), 083003 (2022)](https://doi.org/10.1364/OE.451496)/[arXiv:2112.08035 (2021)](https://arxiv.org/abs/2112.08035)
+ TO FOCUS-MATCH OR NOT TO FOCUS-MATCH INVERSE SPATIALLY OFFSET RAMAN SPECTROSCOPY: A QUESTION OF LIGHT PENETRATION
G.E. Shillito, L. McMillan, G. D. Bruce, K. Dholakia, [Opt. Express 30, 8876 (2022)](https://doi.org/10.1364/OE.451496)/[arXiv:2112.08877](https://arxiv.org/abs/2112.08877)

## TODO's
The current TODO list of planned features and current bugs can be found [here](TODO.md).