module load intel
module load mkl
module load anaconda/3/2018.12
f2py -c --compiler=intel --fcompiler=intelem src/ECRad_python.f90 -m ECRad_python -I$PWD/$SYS/modiOMP/ --build-dir $SYS --f90flags="-qopenmp -fpp -DOMP -O2 -mkl=parallel" -liomp5 -L$MKL_LIB -L$PWD/fitpack/ -L$PWD/odepack/ -L$PWD/$SYS -lECRadOMP -lFITPack -lODEPack -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
