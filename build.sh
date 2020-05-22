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
mkdir $SYS
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
#f2py -c --fcompiler=gnu95 ../src/ECRad_python.f90 --build-dir modiOMP -m ECRad_python -ImodiOMP/ -lECRadOMP --f90flags=-openmp -lgomp -L../../netlib/fitpack/ -lfit
#f2py -c --fcompiler=intelem src/ECRad_python.f90 -m ECRad_python -I$SYS/modiOMPUSE3D/ --f90flags="-qopenmp -fpp -DOMP -DUSE3D -O2" -liomp5 -L../netlib/fitpack/ -L../Mconf/lib/ -L$SYS -lECRadOMPUSE3D -lfit -lmconf64
#f2py -c --fcompiler=intelem src/ECRad_python.f90 -m ECRad_python -I$SYS/modiUSE3D/ --f90flags="-fpp -DUSE3D -g -traceback" -liomp5 -L../netlib/fitpack/ -L../Mconf/lib/ -L$SYS -lECRadUSE3Ddb -lfit -lmconf64