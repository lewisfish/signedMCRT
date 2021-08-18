#!/bin/bash

cd src
gfortran-10 utils.f90 vector_class.f90  sdfsMod.f90 test.f90 -cpp -fopenmp -freal-4-real-8 && ./a.out 
cd ..
python render.py --file cone.dat
python render.py --file bend.dat
python render.py --file elongate.dat
python render.py --file twist.dat
python render.py --file displacement.dat