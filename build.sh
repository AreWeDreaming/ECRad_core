#!/bin/bash

if [ $SYS == "amd64_sles15" ]
  then
  	module purge
  	module load intel/19.0.3
  	module load git
  	module load anaconda/2/2018.12
  	module load hdf5-serial
  	module load netcdf-serial
fi
echo "Type g for gfortran or anything else for intel"
read COMPILERINP
if [ $COMPILERINP == "g" ]
  then
  export COMPILER="g"
else
  export COMPILER="i"
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