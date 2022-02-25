# signedMCRT
Use of signed distance fields in Monte Carlo Radiative Transfer.
This allows modelling of smooth surfaces with out the need to use triangle or similar meshes.

Pre-print describing operation of signedMCRT: https://arxiv.org/abs/2112.08035

Pre-print of an application of signedMCRT: https://arxiv.org/abs/2112.08877


## Instructions

Code prerequisites: fortran 2018 compliant complier e.g gfortran-10 or intel oneAPI fortran.
Only tested on Linux. May run on Mac, but unlikly to run on Windows without some changes.

To run the code you can use:
  - The [fortran package manager](https://fpm.fortran-lang.org/en/index.html) by running; fpm @run
  
## Publication
The code in this repo forms the basis for the following publication:

Meshless Monte Carlo Radiation Transfer Method for Curved Geometries using Signed Distance Functions, currently in review.
Data created with this code can be found at: https://doi.org/10.5281/zenodo.5780513

## Example of SDF models produced by signedMCRT

SDF of University of St Andrews crest. SVG of crest was simplfied in Inkscape and then exported using svg.f90 to line segments and assiagned optical properties before being illuminated uniformly by light.
![Image of SDF model of university crest](https://github.com/lewisfish/signedMCRT/raw/main/crest-sdf-svg.png)

Comparison of voxel and SDF models of a glass bottle with scattering contents. The voxel model exhibts un-physical reflections and refractions.
![Comparison of light distrbution in voxel and SDF bottle](https://github.com/lewisfish/signedMCRT/raw/main/georgie_compare_sdf_vs_voxel.png)

Comparison of voxel, interpolated surface normals, and SDFs model of a sphere being illuminated by a tophat laser beam.
![Comparison of sMCRT to voxel and smoothed surface normal method](https://github.com/lewisfish/signedMCRT/raw/main/sdf_vs_voxel_sphere%20(1).png)

## SDF information
List of SDF functions
https://iquilezles.org/www/articles/distfunctions/distfunctions.htm

What can be modelled with SDF's
https://iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm
