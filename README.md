# signedMCRT
Use of signed distance fields in Monte Carlo Radiative Transfer.
This allows modelling of smooth surfaces with out the need to use triangle or similar meshes.

## Instructions

Code prerequisites: fortran 2018 compliant complier e.g gfortran-10 or intel oneAPI fortran.
Only tested on Linux. May run on Mac, but unlikly to run on Windows without some changes.

To run the code you can use:
  - The [fortran package manager](https://fpm.fortran-lang.org/en/index.html) by running; fpm @run
  - Or the custom install script ./install.sh
  
## Publication
The code in this repo forms the basis for the following publication:

Meshless Monte Carlo Radiation Transfer Method for CurvedGeometries using Signed Distance Functions, currently in review.
Data created with this code can be found at: https://doi.org/10.5281/zenodo.5780513

## Example of SDF model produced by signedMCRT

![Image of SDF model in blender](https://github.com/lewisfish/signedMCRT/raw/main/blender-omg-smooth.png)

## SDF information
List of SDF functions
https://iquilezles.org/www/articles/distfunctions/distfunctions.htm

What can be modelled with SDF's
https://iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm
