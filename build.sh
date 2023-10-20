#!/bin/bash
export F90=gfortran
export FC=gfortran
export F77=gfortran
export COMPILER=GNU
make clean # Just for good measure
cd fitpack
make clean
make
cd ../odepack
make clean
make
cd ../
make lib F2PY_wrapper DEBUG=True COMPILER=GNU
make lib F2PY_wrapper COMPILER=GNU
make lib F2PY_wrapper OPEN_MP=True COMPILER=GNU
make lib F2PY_wrapper OPEN_MP=True DEBUG=True COMPILER=GNU
rm -f id
python3 -m build -n -x
python3 -m pip install --no-deps .
# make USE_3D=True DEBUG=True
# make USE_3D=True
# make USE_3D=True DEBUG=True OPEN_MP=True
# make USE_3D=True OPEN_MP=True