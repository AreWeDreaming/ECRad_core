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
#ifort -module /afs/ipp-garching.mpg.de/home/s/sdenk/ECRad_testing/augd_ecrad/amd64_sles15/modiIDA -I/afs/ipp/cips/ipp/data_analysis/lib_std/$SYS/modi -fpp -DOMP -o $SYS/ECRad_python.o -c -O2 -fp-model source -qopenmp -axavx /afs/ipp-garching.mpg.de/home/s/sdenk/ECRad_testing/augd_ecrad/src/ECRad_python.f90
#f2py -c --fcompiler=intelem --ccompiler=intelem --f90exec=/mpcdf/soft/SLE_15/packages/sandybridge/intel/19.0.3/bin/ifort src/ECRad_python.pyf $SYS/ECRad_python.o $SYS/libECRadIDA.a -m ECRad_python -module $SYS/modiIDA/ -L$SYS/ -lECRadIDA
#f2py -c src/ECRad_python.pyf $SYS/ECRad_python.o $SYS/libECRadIDA.a -m ECRad_python -I$SYS/modiIDA/ -L$SYS/ -lECRadIDA -module $SYS/modiIDA/ 
#f2py -c -I$SYS/modiIDA/ -L$SYS/ -lECRadIDA -m ECRad_python  src/ECRad_python.f90
#f2py -c --noarch --f90exec=/mpcdf/soft/SLE_15/packages/sandybridge/intel/19.0.3/bin/ifort --fcompiler=intel src/ECRad_python.f90 $SYS/libECRadIDA.a