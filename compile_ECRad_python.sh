module load intel
module load mkl
module load anaconda/3/2018.12
f2py -c --fcompiler=intelem ../src/ECRad_python.f90 -m ECRad_python -I$SYS/modiOMPUSE3D/ --f90flags="-qopenmp -fpp -DOMP -O2" -liomp5 -L../fitpack/ -L../odepack/ -L$SYS -lECRadOMP -lFITPack -lODEPack