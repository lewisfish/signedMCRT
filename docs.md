---
project: signedMCRT
src_dir: ./src
         ./app
output_dir: ./docs
media_dir: ./images
project_github: https://github.com/lewisfish/signedMCRT
project_download: https://github.com/lewisfish/signedMCRT/releases/latest
summary: ![signedMCRT](|media|/sMCRT_logo.png)<br> Monte Carlo radiation transfer (MCRT) using Signed Distance functions (SDF)
author: Lewis McMillan
github: https://github.com/lewisfish
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
search: true
sort: alpha
fpp_extensions: fpp
preprocess: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
md_extensions: markdown.extensions.toc
               markdown_checklist.extension
               ford.md_striped_table
page_dir: ./old_docs
---

--------------------
[TOC]

Brief description
-----------------

A Monte Carlo radiation transfer code with signed distance functions representing the geometry, written in modern Fortran.
This allows modelling of smooth surfaces with out the need to use triangle or similar meshes.


Installation
------------

To build signedMCRT, the only current method is using [FPM](https://fpm.fortran-lang.org/).
FPM can be easily installed on any platform, and is simple to use to pull all dependencies, and build and compile signedMCRT.
We also provide several commands via FPM response file ([found here](https://github.com/lewisfish/signedMCRT/blob/main/fpm.rsp)), to enable the use of OpenMP, other compliers, and various debug modes.

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


References
----------

SignedMCRT has so far been used in 2 papers:

+ MESHLESS MONTE CARLO RADIATION TRANSFER METHOD FOR CURVED GEOMETRIES USING SIGNED DISTANCE FUNCTIONS
L. McMillan, G. D. Bruce, K. Dholakia, [J. Biomed. Opt. 27(8), 083003 (2022)](https://doi.org/10.1117/1.JBO.27.8.083003)/[arXiv:2112.08035 (2021)](https://arxiv.org/abs/2112.08035)
+ TO FOCUS-MATCH OR NOT TO FOCUS-MATCH INVERSE SPATIALLY OFFSET RAMAN SPECTROSCOPY: A QUESTION OF LIGHT PENETRATION
G.E. Shillito, L. McMillan, G. D. Bruce, K. Dholakia, [Opt. Express 30, 8876 (2022)](https://doi.org/10.1364/OE.451496)/[arXiv:2112.08877](https://arxiv.org/abs/2112.08877)

License
-------

The signedMCRT source code and related files and documentation are
distributed under a permissive free software license (MIT).

