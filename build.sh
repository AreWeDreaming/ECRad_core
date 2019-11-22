#!/bin/tcsh
module purge
if ($SYS == "amd64_sles11") then
  module load intel/17.0
  mpdule load git
else if ($SYS == "amd64_sles12") then
  module load intel/17.0
  mpdule load git
else if ($SYS == "amd64_sles15") then
  module load intel/19.0.3
  module load git
  module load anaconda/2/2018.12
endif
rm id
git rev-parse HEAD > id
make clean # Just for good measure
make DEBUG=True IDA=True
make IDA=True
make DEBUG=True
make
make OPENMP=True
make OPENMP=True DEBUG=True
module load hdf5-serial
module load netcdf-serial
make USE_3D=True DEBUG=True
make USE_3D=True
make USE_3D=True DEBUG=True OPENMP=True
make USE_3D=True OPENMP=True
f2py -c --fcompiler=intelem ../src/ECRad_python.f90 --build-dir modiOMP -m ECRad_python -ImodiOMP/ -lECRadOMP --f90flags=-openmp -lgomp -L../../netlib/fitpack/ -lfit
#f2py -c --fcompiler=intelem ../src/ECRad_python.f90 --build-dir modiOMP -m ECRad_python -ImodiOMP/ -lECRadOMP --f90flags=-qopenmp -liomp5 -L../../netlib/fitpack/ -lfit