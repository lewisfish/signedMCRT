# Documentation

This document is the incomplete documentation of **signedMCRT**.

## Build system

To build signedMCRT, the only current method is using [FPM](https://fpm.fortran-lang.org/en/index.html).
FPM can be easily installed on any platform, and is simple to use to pull all dependencies, and build and compile signedMCRT.
We also provide several commands via FPM response file ([found here](fpm.rsp)), to enable the use of OpenMP, other compliers, and various debug modes.

## Running the code

The code is run using FPM. To run on a single core with no debug flags enabled:

```
fpm run
```

To run on all available threads on current computer with no debug flags:

```
fpm @runmp
```

To run the code on one thread with all debug flags enabled:

```
fpm @debug
```

To run the code on all threads with all debug flags enabled:

```
fpm @debugmp
```

Please see ([here](fpm.rsp)) for other possible options.

## Dependencies

Below is the current list of dependencies:

* [test drive](https://github.com/fortran-lang/test-drive)
* [Fortran TEV Bindings](https://github.com/lewisfish/fortran_tev_bindings)
* [stdlib](https://github.com/fortran-lang/stdlib)
* [stb_image](https://github.com/lewisfish/fortran_stb_image)
* [Fortran Utilities](https://github.com/lewisfish/fortran_utilities)

Test drive is used to run all tests.
Fortran [TEV](https://github.com/Tom94/tev/) Bindings is used to interface with TEV, to show live slices of fluences as the simulation is run, which is handy for debugging purposes.
Stdlib is a collection of routines purposed for inclusion within the Fortran standard. Stdlib is used here for it's loadtxt function to load arbitrary plain text data into arrays. More of stdlib may be used in future.
Fortran_stb_Image is used to load images into arrays. Fortran_stb_image are the Fortran bindings for [stb_image](https://github.com/nothings/stb).
Finally, Fortran Utilities is my personal collection of useful Fortran utilities such as mathematical functions, or progress bars.

## Config file
signedMCRT uses TOML as it's configuration file format.
Documentation of the input file format can be found in [here](config.md)


## Code Structure

See below for a brief description of each source code file.
See [this](relations.pdf) for rat's nest diagram of the relationships and dependencies between each source code file. 
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

Source code can be found [here](/src/geometryMod.f90)

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

### opticalPropertiesMod.f90

Source code can be found [here](/src/opticalProps/opticalProperties.f90)


### photon.f90

This source file contains the photon type, all the photon launch routines for different light sources, and the scattering code.

Below are the current types of light sources available. Check [here](config.md) for parameters needed for each light source.
- uniform
- pencil
- annulus
- focus
- point
- circular
- SLM (2D image source)
- double slit
- square aperture

Source code can be found [here](/src/photon.f90)

### parse.f90

This file contains all the logic for parsing the toml input file using the toml-f library.
parse_params is the only public function in this module. This function returns a dictionary (toml_table) of parameters needed to setup the simulation and an array of detectors if any are defined.

Any additions to the toml input file should be added in here.

Source code can be found [here](/src/parse.f90)

### piecewise.f90

This file contains the piecewise abstract type, for sampling from constants, 1D or 2D arrays. Inspired by [PBRT](https://www.pbr-book.org/) piecewise class.
Currently, the following public types are defined:

* Constant. Used in the case where there is only one value.
* 1D. Used in the case where there is a spectrum
* 2D. Used in the case where SLM or other image based source types are needed.

The piecewise type ensures that there is a method (sample) that can be called on all inherited types, e.g

call 2Dimage%p%sample(x, y)

will return a position (x,y) from where to release a photon.

This class can be used to have multi-spectral or single valued wavelength, or used as a 2D image input source i.e SLMs.
NOTE: optical properties are not currently adjusted on wavelength change.

Source code can be found [here](/src/piecewise.f90)

### random_mod.f90

This module defines a set of functions that return random numbers in different distributions.

- ran2. Returns a single float uniformly in the range [0, 1)
- ranu. Return a single float uniformly in the range [a, b)
- randint. Returns a single integer uniformly in the range [a, b)
- rang. Returns a single float from a Gaussian distribution with mean *avg* and std *sigma*.
- init_rng. Seeds the internal random number generator with a reproducible seed.

Source code can be found [here](/src/random_mod.f90)

### sdf_base.f90

This module defines the signed distance function (SDF) abstract type, sdf_base type, and model type.
The SDF abstract type defines the optical properties of an SDF (mus, mua, kappa, albedo, hgg, g2,and n), as well as a transform (4x4 matrix), and the layer ID code of the SDF. The SDF abstract type also provides an abstract interface (evaluate) which each inheriting function must implement. This evaluate function is the heart of the SDF implementation. Each individual evaluate is the direct implementation of that SDF, e.g. that function defines the mathematical SDF. For more information on SDFs, check out Inigo Quilez's [website](https://iquilezles.org/articles/) from which most of the below SDFs and transforms have been taken.

The sdf_base type acts as a container for all SDFs derived from the abstract SDF type.

Source code can be found [here](/src/sdfs/sdf_base.f90)

### sdfHelpers.f90

Collection of helper functions for SDFs:

This module defines transforms that can be applied to each SDF:

- Rotate_{x,y,z}
- Translate
- RotationAlign (not tested)
- RotMat (not tested)
- Identity
- SkewSymm

Source code can be found [here](/src/sdfs/sdfHelpers.f90)

### sdfModifiers.f90

This module defines transforms that can be applied to each SDF:

- Union
- Intersection
- Subtraction
- Displacement
- Bend
- Twist
- Elongate
- Repeat
- Extrude
- Revolution
- Onion
- Translate

Source code can be found [here](/src/sdfs/sdfModifiers.f90)

### sdfs.f90

This module defines the signed distance function (SDF) abstract type and all types that inherit from it.
The SDF abstract type defines the optical properties of an SDF (mus, mua, kappa, albedo, hgg, g2,and n), as well as a transform (4x4 matrix), and the layer ID code of the SDF. The SDF abstract type also provides an abstract interface (evaluate) which each inheriting function must implement. This evaluate function is the heart of the SDF implementation. Each individual evaluate is the direct implementation of that SDF, e.g. that function defines the mathematical SDF. For more information on SDFs, check out Inigo Quilez's [website](https://iquilezles.org/articles/) from which most of the below SDFs and transforms have been taken.

- cylinder
- sphere
- box
- torus
- cone
- triprism (triangular prism)
- capsule
- plane
- segment
- egg

**This is the module the user should import to other module not sdf_base!**

Source code can be found [here](/src/sdfs.f90)

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

### setup.f90

This file sets up some simulations variables and assigns the geometry for the simulation.

Source code can be found [here](/src/setup.f90)

### setupGeometry.f90

Add any new geometry setups here.
Source code can be found [here](/src/setupGeometry.f90)

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

## Plotting Results

To view the output of simulations you can use [this](https://github.com/lewisfish/data_cube_viewer).
Alternatively to customise the plot you can adjust the following [script](../tools/plot_nrrd.py).

## Monte Carlo Radiation Transfer (MCRT) method

Please see my [thesis](main.pdf) for an overview of the MCRT method

## Citation

SignedMCRT has so far been used in 2 papers:

+ MESHLESS MONTE CARLO RADIATION TRANSFER METHOD FOR CURVED GEOMETRIES USING SIGNED DISTANCE FUNCTIONS
L. McMillan, G. D. Bruce, K. Dholakia, [J. Biomed. Opt. 27(8), 083003 (2022)](https://doi.org/10.1117/1.JBO.27.8.083003)/[arXiv:2112.08035 (2021)](https://arxiv.org/abs/2112.08035)
+ TO FOCUS-MATCH OR NOT TO FOCUS-MATCH INVERSE SPATIALLY OFFSET RAMAN SPECTROSCOPY: A QUESTION OF LIGHT PENETRATION
G.E. Shillito, L. McMillan, G. D. Bruce, K. Dholakia, [Opt. Express 30, 8876 (2022)](https://doi.org/10.1364/OE.451496)/[arXiv:2112.08877](https://arxiv.org/abs/2112.08877)

## TODO's
The current TODO list of planned features and current bugs can be found [here](TODO.md).
