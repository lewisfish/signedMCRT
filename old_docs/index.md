---
title: main
---
# Documentation

This document is the incomplete documentation of **signedMCRT**.

## Build system

To build signedMCRT, the only current method is using [FPM](https://fpm.fortran-lang.org/en/index.html).
FPM can be easily installed on any platform, and is simple to use to pull all dependencies, and build and compile signedMCRT.
We also provide several commands via FPM response file ([found here](https://github.com/lewisfish/signedMCRT/blob/main/fpm.rsp)), to enable the use of OpenMP, other compliers, and various debug modes.

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

Please see ([here](https://github.com/lewisfish/signedMCRT/blob/main/fpm.rsp)) for other possible options.

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
Documentation of the input file format can be found in [here](https://lewismcmillan.com/signedMCRT/page/config.html)

## Plotting Results

To view the output of simulations you can use [this](https://github.com/lewisfish/data_cube_viewer).
Alternatively to customise the plot you can adjust the following [script](https://github.com/lewisfish/signedMCRT/blob/main/tools/plot_nrrd.py).

## Monte Carlo Radiation Transfer (MCRT) method

Please see my [thesis](https://github.com/lewisfish/my_amazing_thesis/blob/master/main.pdf) for an overview of the MCRT method.

## Citation

SignedMCRT has so far been used in 2 papers:

+ MESHLESS MONTE CARLO RADIATION TRANSFER METHOD FOR CURVED GEOMETRIES USING SIGNED DISTANCE FUNCTIONS
L. McMillan, G. D. Bruce, K. Dholakia, [J. Biomed. Opt. 27(8), 083003 (2022)](https://doi.org/10.1117/1.JBO.27.8.083003)/[arXiv:2112.08035 (2021)](https://arxiv.org/abs/2112.08035)
+ TO FOCUS-MATCH OR NOT TO FOCUS-MATCH INVERSE SPATIALLY OFFSET RAMAN SPECTROSCOPY: A QUESTION OF LIGHT PENETRATION
G.E. Shillito, L. McMillan, G. D. Bruce, K. Dholakia, [Opt. Express 30, 8876 (2022)](https://doi.org/10.1364/OE.451496)/[arXiv:2112.08877](https://arxiv.org/abs/2112.08877)

## TODO's
The current TODO list of planned features and current bugs can be found [here](https://lewismcmillan.com/signedMCRT/page/TODO.html).