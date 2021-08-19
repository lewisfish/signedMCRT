#!/bin/bash
set -e

cd src
gfortran-10 utils.f90 vector_class.f90  sdfsMod.f90 test.f90 -cpp -fopenmp -freal-4-real-8 && ./a.out 
cd ..

echo "Done generating voxels"
python render.py --file cone.dat
python render.py --file triprisim.dat
python render.py --file torus.dat
python render.py --file box.dat
python render.py --file cylinder.dat
python render.py --file sphere.dat
python render.py --file capsule.dat
python render.py --file plane.dat
echo "Done rendering shapes"


python render.py --file bend.dat
python render.py --file elongate.dat
python render.py --file twist.dat
python render.py --file displacement.dat
python render.py --file repeat.dat
echo "Done rendering operations"

python render.py --file union.dat
python render.py --file smoothunion.dat
python render.py --file intersection.dat
python render.py --file subtraction.dat
echo "Done"
