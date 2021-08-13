#!/bin/bash

if [[ $HOSTNAME == *"mpg"* ]]
  then
  module purge
  export PYTHONPATH
  module load intel
  module load mkl
  module load texlive
  module load anaconda/3/2020.02
  module load git
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64/
  conda activate ECRad_conda
elif [[ $HOSTNAME == *"cm.cluster"* ]]
  then
  module purge
  export PYTHONPATH
  module use /home/software/psfc/modulefiles/
  module load psfc/config
  module load slurm
  module load gcc
  module load intel/2017-01
  module load psfc/netcdf/intel-17/4.4.1.1
  module load psfc/mkl/17
  module load engaging/git
  module load anaconda3/2020.11
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
  conda activate ECRad_conda
elif [[ $HOSTNAME == *"iter"* ]]; then
  module purge
  export PYTHONPATH
  module load IMAS
  module load texlive
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib/python3.8/site-packages/wx/
fi