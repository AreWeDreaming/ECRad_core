# Start from clean environment
module purge

# IMAS and iWrap
module load IMAS iWrap

# Compiler
export F90=ifort
export FC=ifort

# For debugging, just in case
module load TotalView

# Actor folder
export ACTOR_FOLDER=~/public/PYTHON_ACTORS
mkdir -p $ACTOR_FOLDER

# Need to remove the stack limit to avoid segmentation fault inside codes
ulimit -Ss unlimited

# Libraries needed for code parameters
module load XMLlib/3.3.1-intel-2020b

module load MUSCLE3

# EXTEND PYTHON PATH AND AVOID DOUBLONS
export PYTHONPATH=$ACTOR_FOLDER:$PYTHONPATH
export PYTHONPATH="$(perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, $ENV{PYTHONPATH}))')"

