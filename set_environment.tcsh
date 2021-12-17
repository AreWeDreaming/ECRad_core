#!/usr/bin/tcsh

if ($DOMAIN =~ *"mpg"* ) then
  module purge
  module load intel
  module load mkl
  module load texlive
  module load anaconda/3/2020.02
  module load git
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
  if ( -f "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/etc/profile.d/conda.csh" ) then
      source "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/etc/profile.d/conda.csh"
  else
      setenv PATH "/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2020.02/bin:$PATH"
  endif
# <<< conda initialize <<<
  conda activate ECRad_conda
  if ($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$MKLROOT/lib/intel64/
  else
  	setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
  endif
endif