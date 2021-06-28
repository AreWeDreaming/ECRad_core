#!/bin/bash

if [ $HOSTNAME == *"mpg.de"* ]
  then
  module purge
  module load texlive
  module load intel
  module load mkl
  module load anaconda/3/2020.02
  module load git
  module load hdf5-serial
  module load netcdf-serial
  setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
  export COMPILER="i"
else if [$HOSTNAME == *"cm.cluster"*]
  then
  module use /home/software/psfc/modulefiles/
  module load psfc/config
  module load anaconda3/2020.11
  module load psfc/netcdf/intel-17/4.4.1.1
  module load intel
  module load psfc/pgplot/5.2.2
  module load texlive
  module load engaging/git
  # >>> conda initialize >>>
  # !! Contents within this block are managed by 'conda init' !!
  __conda_setup="$('/home/software/anaconda3/2020.11/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
  if [ $? -eq 0 ]; then
      eval "$__conda_setup"
  else
      if [ -f "/home/software/anaconda3/2020.11/etc/profile.d/conda.sh" ]; then
          . "/home/software/anaconda3/2020.11/etc/profile.d/conda.sh"
      else
          export PATH="/home/software/anaconda3/2020.11/bin:$PATH"
      fi
  fi
  unset __conda_setup
  # <<< conda initialize <<<
  conda env create -f ECRad_env.yml
  conda activate ECRad_conda
  setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
  export COMPILER="i"
else
  echo "Type g for gfortran or anything else for intel"
  read COMPILERINP
  if [ $COMPILERINP == "g" ]
    then
    export COMPILER="g"
  else
    export COMPILER="i"
  fi
fi
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
make USE_3D=True DEBUG=True
make USE_3D=True
make USE_3D=True DEBUG=True OPEN_MP=True
make USE_3D=True OPEN_MP=True