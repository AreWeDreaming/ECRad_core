#!/bin/bash

make clean
make clean -C fitpack
make clean -C odepack

make -C fitpack
make -C odepack
make DEBUG=True IMAS=True

iwrap -f iwrap/ecrad.yaml -i ~/public/PYTHON_ACTORS
