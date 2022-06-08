#!/usr/bin/tcsh

if ($DOMAIN =~ *"mpg"* ) then
  module purge
  module load texlive
  module load intel/21.3.0
  module load mkl/2021.3
  module load anaconda/3/2021.11
  module load git
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
  if ( -f "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.11/etc/profile.d/conda.csh" ) then
      source "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.11/etc/profile.d/conda.csh"
  else
      setenv PATH "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.11/bin:$PATH"
  endif
# <<< conda initialize <<<
  conda activate ECRad_conda
  bash -c 'source $INTEL_HOME/setvars.sh ; exec csh'
endif