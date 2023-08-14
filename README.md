![SignedMCRT logo made in signedMCRT](https://github.com/lewisfish/signedMCRT/raw/main//images/sMCRT_logo.png)
============
[![DOI](https://zenodo.org/badge/390770167.svg)](https://zenodo.org/badge/latestdoi/390770167) [![codecov](https://codecov.io/github/lewisfish/signedMCRT/branch/main/graph/badge.svg?token=U402PQWWUY)](https://codecov.io/github/lewisfish/signedMCRT)
![workflow](https://github.com/lewisfish/signedMCRT/actions/workflows/build.yml/badge.svg
)


Use of signed distance fields in Monte Carlo Radiative Transfer.
This allows modelling of smooth surfaces with out the need to use triangle or similar meshes.

Paper describing operation of signedMCRT: [https://doi.org/10.1117/1.JBO.27.8.083003](https://doi.org/10.1117/1.JBO.27.8.083003)

Paper describing of an application of signedMCRT: [https://doi.org/10.1364/OE.451496](https://doi.org/10.1364/OE.451496)


## Instructions

Code prerequisites: Fortran 2018 compliant complier e.g gfortran-9+ or intel oneAPI Fortran.
Only tested on Linux and Mac. Unlikely to run on Windows without some changes.

To run the code you can use:
  - The [Fortran package manager](https://fpm.fortran-lang.org/en/index.html) by running; "fpm run" for single core uage
  or "fpm @runmp" for multicore usage.

See [this](docs/main.md) for more detailed instructions and details on the code and build systems used.
  
## Publication
The code in this repo forms the basis for the following publication:

Meshless Monte Carlo Radiation Transfer Method for Curved Geometries using Signed Distance Functions, currently in review.
Data created with this code can be found at: https://doi.org/10.5281/zenodo.5780513

## Example of SDF models produced by signedMCRT

SDF of University of St Andrews crest. SVG of crest was simplified in Inkscape and then exported using svg.f90 to line segments and assigned optical properties before being illuminated uniformly by light.
![Image of SDF model of university crest](https://github.com/lewisfish/signedMCRT/raw/main/images/crest-sdf-svg.png)

Comparison of voxel and SDF models of a glass bottle with scattering contents. The voxel model exhibits un-physical reflections and refractions.
![Comparison of light distribution in voxel and SDF bottle](https://github.com/lewisfish/signedMCRT/raw/main/images/georgie_compare_sdf_vs_voxel.png)

Comparison of voxel interpolated surface normals, mesh based Monte Carlo (MMC), and SDFs model of a sphere being illuminated by a tophat laser beam.
![Comparison of sMCRT to voxel and smoothed surface normal method](https://github.com/lewisfish/signedMCRT/raw/main/images/sdf_vs_mmc_aptran%20(1).png)

## SDF information
List of SDF functions
https://iquilezles.org/www/articles/distfunctions/distfunctions.htm

What can be modelled with SDF's
https://iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm
