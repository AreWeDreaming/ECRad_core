#!/bin/bash

if [ $HOSTNAME == *"mpg.de"* ]
  then
  module purge
  module load texlive
  module load intel
  module load mkl
  module load anaconda/3/2020.02
  module load git
  module load hdf5-serial
  module load netcdf-serial
  setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
else
  echo "Type g for gfortran or anything else for intel"
  read COMPILERINP
  if [ $COMPILERINP == "g" ]
    then
    export COMPILER="g"
  else
    export COMPILER="i"
  fi
fi
rm id
git rev-parse HEAD > id
make clean # Just for good measure
cd fitpack
make clean
make
cd ../odepack
make clean
make
cd ../
make DEBUG=True IDA=True
make IDA=True
make DEBUG=True
make
make OPEN_MP=True
make OPEN_MP=True DEBUG=True
make USE_3D=True DEBUG=True
make USE_3D=True
make USE_3D=True DEBUG=True OPEN_MP=True
make USE_3D=True OPEN_MP=True