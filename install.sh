#!/bin/bash

function showhelp
{

  echo 'Usage: ./install.sh [option] [option]...'
  echo 'Compile and run mcgrid.'
  echo 
  echo '   -h, --help            Shows this dialog.'
  echo '   -n, --cores           Compiles and run code on n cores. Default is 1 core.'
  echo '   -d, --debug           Compiles and runs code with debug flags on n core.'
  echo '   -m, --make            Compiles the code with warning enabled.'
  echo '   -omp, --openmp         Use openmp as parallel library.'
  echo '   -mpi, --mpi           Use MPI as parallel library.'

}

function makebuild
{
  if [ "$comp" = 'gnu' ] && [ "$openmp" = 1 ];then
      string="FCOMP=gfortran"
  elif [ "$comp" = 'gnu' ] && [ "$mpi" = 1 ];then
    string="FCOMP=/usr/bin/mpifort"
  elif [ "$comp" = 'intel' ] && [ "$openmp" = 1 ];then
      string="FCOMP=ifort"
  elif [ "$comp" = 'intel' ] && [ "$openmp" = 1 ];then
      string="FCOMP=mpiifort"
  else
    string="FCOMP=gfortran"
  fi

  if [ "$debug" = 1 ];then
    make clean && make debug $string
  elif [ "$NUM_CORES" = 0 ];then
    make clean && make build $string
  else
    if [ "$comp" = 'gnu' ];then
      if [ "$openmp" = 1 ];then
        export OMP_NUM_THREADS=$NUM_CORES
        make clean && make mp $string
      else
        make clean && make $string
      fi
    elif [ "$comp" = 'intel' ];then
      make clean&& make $string
    fi
  fi
}

function createdirs
{
  if [ ! -d "data" ]; then
      mkdir "data"
      mkdir "data/jmean"
      mkdir "data/im"
      mkdir "data/deposit"
  fi

  if [ ! -d "build" ]; then
     mkdir "build"
  fi
  cd build
  ndirec="$(pwd)"
  cd ..
  if [ ! -d "bin" ]; then
     mkdir "bin"
  fi
  cd bin
  bdirc="$(pwd)"
  cd ..
  cd src
}

function run
{
  for i in *; do
     if [ "${i}" != "${i%.mod}" ];then
        cp "${i}" "$ndirec"
     fi
     if [ "${i}" != "${i%.o}" ];then
        mv "${i}" "$ndirec"
     fi
  done


  if [ "$NUM_CORES" = "0" ]; then #just make code
      exit 0
  fi

  mv mcgrid "$bdirc" && echo " "&& echo "*****Install complete*****" && echo " "

  # clear
  cd ../bin

  if [ "$NUM_CORES" = "1" ]; then
      ./mcgrid
  else
    if [ $mpi == 1 ];then
      if [ $comp = 'gnu' ];then
        /usr/local/bin/mpirun -n $NUM_CORES ./mcgrid
      elif [ $comp = 'intel' ];then
        mpirun -n $NUM_CORES ./mcgrid
      fi
    else
      ./mcgrid
    fi
  fi
}

#defaults
NUM_CORES=1
debug=0
help=0
openmp=0
mpi=0
comp="gnu"

set -e

createdirs

while [ "$1" != "" ]; do
    case $1 in
        -n | --cores )          NUM_CORES=$2
                                ;;
        -c | --comp )           comp=$2
                                ;;
        -h | --help )           showhelp
                                exit
                                ;;
        -m | --make )           NUM_CORES=0
                                makebuild
                                exit
                                ;;
        -d | --debug )          debug=1
                                ;;
        -omp | --openmp )        openmp=1
                                mpi=0
                                ;;
        -mpi | --mpi )          mpi=1
                                openmp=0
                                ;;
    esac
    shift
done

makebuild
run