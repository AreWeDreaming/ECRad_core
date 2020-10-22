module load intel
module load mkl
module load anaconda/3/2018.12
f2py -c --fcompiler=gnu95  src/ECRad_python.f90 -m ECRad_python --build-dir $SYS -I$SYS/modgOMP/ --f90flags="-O2 -mavx -ffree-form -fpic -ffree-line-length-none" -L/mnt/c/Users/Severin/git/augd_ecrad/fitpack/ -L/mnt/c/Users/Severin/git/augd_ecrad/odepack/ -L/mnt/c/Users/Severin/git/augd_ecrad/$SYS -lECRadOMP -lFITPack -lODEPack