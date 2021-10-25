#!/usr/bin/tcsh
if($HOSTNAME =~ *mpg.de) then
  echo "Current system identified as IPP TOK cluster"
  module purge
  module load texlive
  module load intel
  module load mkl
  module load anaconda/3/2020.02
  module load git
  module load hdf5-serial
  module load netcdf-serial
  setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
  conda env create -f ECRad_env.yml
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
  if ( -f "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/etc/profile.d/conda.csh" ) then
      source "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/etc/profile.d/conda.csh"
  else
      setenv PATH "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/bin:$PATH"
  endif
  conda activate ECRad_conda
  conda env update -f ECRad_env.yml --prune
  set COMPILER = "i"
else
  echo "Type g for gfortran or anything else for intel"
  set COMPILERINP = $<
  if($COMPILERINP == "g") then
  	set COMPILER = "g"
  else
	set COMPILER = "i"
  endif
endif

echo "Compiler is $COMPILER"
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
make DEBUG=True IDA=True
make IDA=True
make DEBUG=True
make
make OPEN_MP=True
make OPEN_MP=True DEBUG=True
# make USE_3D=True DEBUG=True
# make USE_3D=True
# make USE_3D=True DEBUG=True OPEN_MP=True
# make USE_3D=True OPEN_MP=True
#f2py -c --fcompiler=gnu95 ../src/ECRad_python.f90 --build-dir modiOMP -m ECRad_python -ImodiOMP/ -lECRadOMP --f90flags=-openmp -lgomp -L../../netlib/fitpack/ -lfit -L../../magconf/lib/ -lmconf64
#f2py -c --fcompiler=intelem src/ECRad_python.f90 -m ECRad_python -I$SYS/modiOMPUSE3D/ --f90flags="-qopenmp -fpp -DOMP -DUSE3D -O2" -liomp5 -L../netlib/fitpack/ -L../Mconf/lib/ -L$SYS -lECRadOMPUSE3D -lfit -lmconf64
#f2py -c --fcompiler=intelem src/ECRad_python.f90 -m ECRad_python -I$SYS/modiUSE3D/ --f90flags="-fpp -DUSE3D -g -traceback" -liomp5 -L../netlib/fitpack/ -L../Mconf/lib/ -L$SYS -lECRadUSE3Ddb -lfit -lmconf64
