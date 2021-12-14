# signedMCRT
Use of signed distance fields in Monte Carlo Radiative Transfer.
This allows modelling of smooth surfaces with out the need to use triangle or similar meshes.

## Publication
The code in this repo forms the basis for the following publication:

Meshless Monte Carlo Radiation Transfer Method for CurvedGeometries using Signed Distance Functions, currently in review.

## SDF information
List of SDF functions
https://iquilezles.org/www/articles/distfunctions/distfunctions.htm

What can be modelled with SDF's
https://iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm

## Example of SDF model produced by signedMCRT

"Ray traced" in signedMCRT where the letters are highly absorbing in a non-scattering/absorbing medium.

![Image of SDF model ray traced in signedMCRT](https://github.com/lewisfish/signedMCRT/raw/main/omg_logo_sdf.png)

Blender render

![Image of SDF model in blender](https://github.com/lewisfish/signedMCRT/raw/main/blender-omg-smooth.png)
