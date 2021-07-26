#!/bin/tcsh

if ($HOSTNAME =~ *"mpg"* ) then
  module purge
  module load intel
  module load mkl
  module load texlive
  module load anaconda/3/2020.02
  module load git
  if ($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$MKLROOT/lib/intel64/
  else
  	setenv LD_LIBRARY_PATH $MKLROOT/lib/intel64/
  endif
endif