#!/bin/bash
$PYTHON -m build -n -x
$PYTHON -m pip install --no-deps .
make clean # Just for good measure
cd fitpack
make clean
make
cd ../odepack
make clean
make
cd ../
make DEBUG=True
make
make OPEN_MP=True
make OPEN_MP=True DEBUG=True
rm id
git rev-parse HEAD > id
# make USE_3D=True DEBUG=True
# make USE_3D=True
# make USE_3D=True DEBUG=True OPEN_MP=True
# make USE_3D=True OPEN_MP=True