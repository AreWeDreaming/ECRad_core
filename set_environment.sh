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
elif [[ $HOSTNAME == *"cm.cluster"* ]] || [[$HOSTNAME == "eofe8"]]
  then
  module purge
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
  export PYTHONPATH
  conda activate ECRad_conda
elif [[ $HOSTNAME == *"iter"* ]]; then
  module purge
  export PYTHONPATH
  module load IMAS
  module load texlive
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib/python3.8/site-packages/wx/
elif [[ $HOSTNAME == *"iris"* ]]; then
  module purge
  module load omfit
  module load gcc-9.2.0
  export SYS=CENTOS
  # >>> conda initialize >>>
  # !! Contents within this block are managed by 'conda init' !!
  __conda_setup="$('/fusion/projects/codes/atom/omfit_v3.x/atom/miniconda3_3.2021.10.9/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
  if [ $? -eq 0 ]; then
      eval "$__conda_setup"
  else
      if [ -f "/fusion/projects/codes/atom/omfit_v3.x/atom/miniconda3_3.2021.10.9/etc/profile.d/conda.sh" ]; then
          . "/fusion/projects/codes/atom/omfit_v3.x/atom/miniconda3_3.2021.10.9/etc/profile.d/conda.sh"
      else
          export PATH="/fusion/projects/codes/atom/omfit_v3.x/atom/miniconda3_3.2021.10.9/bin:$PATH"
      fi
  fi
  unset __conda_setup
  # <<< conda initialize <<<
  conda activate ECRad_conda
else
  conda activate ECRad_conda
fi