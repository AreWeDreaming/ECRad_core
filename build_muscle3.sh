#!/bin/bash
make clean
make clean -C fitpack
make clean -C odepack

make -C fitpack
make -C odepack
make OPEN_MP=True MUSCLE3=True IMAS=True
make MUSCLE3=True IMAS=True
make OPEN_MP=True
make
