COMPILER=i
APPLICATION = ECRad

ROOTDIR=$(CURDIR)
ECRadLIBDir=$(ROOTDIR)/Library/$(SYS)
ECRadLIB=$(APPLICATION)
MODECRad=$(ROOTDIR)/Library/$(SYS)/mod$(COMPILER)
obj=	
#Linux
LINUX26 = $SYS
F90 = ifort
F77 = ifort
FFP_FLAGS = -fpp
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)
# Optimization
#F90_FLAGi$(LINUX26) = -c -O
#F90_FLAGi$(LINUX26) = -c -O2 -fp-model source -axavx
# Profiling
#F90_FLAGi$(LINUX26) = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
# Parallelisation
# F90_FLAGi$(LINUX26) =  -c -O2 -fp-model source -openmp-report -openmp -axavx 
#
# Debugging
F90_FLAGi = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
# Debugging parallelisation
# F90_FLAGi$(LINUX26) = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
#F90_FLAGi$(LINUX26) = -c -g -traceback -check bounds  -check all -u -check uninit -fp-model source
F77_FLAG = $(F90_FLAGi) -C --wide
	# CEC-Library
FITPACK = ../netlib/fitpack/lib.a
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB) $(NAGF90LIB) $(NAGFLIB) $(FITPACK) # -qopenmp
LDFLAGS = -z muldefs
MODULES = $(NAG_MOD) -module $(MODECRad)

ECRad: $(OBJJ) $(APPLICATION).o makefile
	echo $(OBJJ$(SYS))
	$(F90) $(LDFLAGS) $(APPLICATION).o ${OBJJ} $(LIBS) \
	-o $(APPLICATION)

$(APPLICATION).o : $(APPLICATION).f90
	$(F90) -c ${MODULES} $(FFP_FLAGS) $(F90FLAGS) $(APPLICATION).f90

clean:
	 rm -f ECRad && rm ECRad.o
