APPLICATION = ECRad
ROOTDIR=$(CURDIR)
ifndef SYS
		SYS = bin
endif
ECRadLIB=$(APPLICATION)
SRCP=$(ROOTDIR)/src
ifdef PREFIX
	ECRadLIBDir=$(PREFIX)/lib
	CONDALIBS = $(ECRadLIBDir)
	ECRad_pythonDir = $(ROOTDIR)/src/ecrad_core
else
	ECRadLIBDir=$(ROOTDIR)/$(SYS)
	ECRad_pythonDir = $(ECRadLIBDir)
endif


ifeq ($(IDA),True)
	STDP=/afs/ipp/cips/ipp/data_analysis/lib_std/$(SYS)
	MODSTDP=$(STDP)/mod$(COMPILER)
	STDPLIB=$(STDP)/libstd$(COMPILER).a
	IDAFLAG = IDA
else
	STDPLIB = $(SRCP)/std_lib.f90
endif
obj=
F2PYEXT_SUFFIX := $(shell python3-config --extension-suffix)
ifeq ($(COMPILER),GNU)
	F90OPTFLAGS = -O2 -mavx -ffree-form -ffree-line-length-none -fPIC
	F90DBGFLAGS = -g -ffree-form -ffree-line-length-none -fPIC -fbacktrace
	F2PYOPTFLAGS = -O2 -mavx -ffree-form -ffree-line-length-none
	F2PYDBGFLAGS = -g -ffree-form -ffree-line-length-none -fbacktrace
	F90PARFLAGS = -fopenmp
	F90PARLIBFLAGS = -lgomp
	FFPFLAGS = -cpp
	MODULEFLAG = 	
	LIBFLAG = -L$(CONDALIBS) -static-libgcc -lopenblas 
	F2PYLIBFLAGS = -L$(CONDALIBS) -lopenblas
	LDFLAGS = -Wl,-rpath=$(CONDALIBS)/lib
	F2PYCOMPILER = gnu95
else
	COMPILER=INTEL
	FFPFLAGS = -fpp -DINTEL
	F90OPTFLAGS = -O2 -fpic -fp-model source -axavx 
	F90DBGFLAGS = -O0 -g -fpic -traceback -shared-intel  #-DTBB_USE_DEBUG -check all -ftrapuv
	F2PYOPTFLAGS = -O2
	F2PYDBGFLAGS = -g -traceback -DTBB_USE_DEBUG
	F90PARFLAGS = -qopenmp
	F90PARLIBFLAGS = 
	MODULEFLAG = -module
	LIBFLAG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
	INCLUDEFLAGS = -I"${MKLROOT}/include"
	F2PYLIBFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
	LDFLAGS = -Wl,-rpath=${MKLROOT}/lib/intel64
	F2PYCOMPILER = intelem
	CC = icc
endif
MKDIR_P = mkdir -p
ifeq ($(IDA),True)
	FFPFLAGS += -DIDA
endif
ifeq ($(OPEN_MP),True)
	FFPFLAGS += -DOMP
	OMPFLAG = OMP
endif
ifeq ($(USE_3D),True)
	FFPFLAGS += -DUSE_3D
	USE3DFLAG = USE3D
endif
#ifeq ($(IDA)$(USE_3D),TrueFalse)
#	FFPFLAGS += -DIDAUSE_2D
#endif
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)
MODECRad=$(ROOTDIR)/$(SYS)/mod$(COMPILER)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)
# Debugging -> enable if DEBUG==		True
ifeq ($(DEBUG),True)
	F90FLAGS = $(F90DBGFLAGS)
	F2PYFLAGS = $(F2PYDBGFLAGS)
	F2PYDBG = --debug-capi
	DB = db
	# Optimized
else
  F90FLAGS = $(F90OPTFLAGS)
  F2PYFLAGS = $(F2PYOPTFLAGS)
  # Profiling
  #F90FLAGS = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir -fpic
endif
ifeq ($(OPEN_MP),True)
  # Parallelisation
  F90FLAGS += $(F90PARFLAGS)
  F2PYFLAGS += $(F90PARFLAGS)
endif
# Debugging parallelisation
#F90FLAGS = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -warn interfaces -fpe3 -qopenmp -qopenmp-report -DTBB_USE_DEBUG
F77FLAGS = $(F90FLAGS) -C
# Libraries
FITPACK = -L$(ROOTDIR)/fitpack/ -lFITPack
ODEPACK = -L$(ROOTDIR)/odepack/ -lODEPack
F2PYLIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(OMPFLAG)$(USE3DFLAG)$(DB) \
	$(NAGF90LIB) $(NAGFLIB) $(FITPACK) $(ODEPACK)
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(OMPFLAG)$(USE3DFLAG)$(DB) \
	$(NAGF90LIB) $(NAGFLIB) $(FITPACK) $(ODEPACK) \
	$(LIBFLAG)
ifeq ($(USE_3D),True)
 	LIBS += $(ROOTDIR)/MConf/lib/libmconf64.a
#   LIBS += $(ROOTDIR)/../magconf/lib/libmconf64.a
#../Mconf/unix/mconf_matlab64.so
	LIBS += -lpthread -lstdc++
	#LIBS += -L${NETCDF_HOME}/lib/  -lnetcdf_c++4  -lnetcdf 
endif
ifeq ($(OPEN_MP),True)
	LIBS += $(F90PARLIBFLAGS)
	F2PYLIBS += $(F90PARLIBFLAGS)
endif
F2PYLIBS += $(F2PYLIBFLAGS)
#LDFLAGS = -z muldefs
ifeq ($(COMPILER),GNU)
	MODULES = $(MODULEFLAG) -J$(MODECRad)
else
	MODULES = $(MODULEFLAG) $(MODECRad) $(INCLUDEFLAGS)
endif
ifeq ($(IDA),True)
	ifneq ($(USE_3D),True)
		MODULES += -I$(MODSTDP)
	endif
endif
ifeq ($(NAG),True)
	MODULES += $(NAG_MOD)
endif

#Targets for library
ifeq ($(IDA),True)
OBJECTS =
else
OBJECTS = std_lib$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o
endif
ifeq ($(USE_3D),True)
OBJECTS = std_lib$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o
endif
OBJECTS += \
	quadrature$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_contour$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	magconfig3D$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_types$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_interpol$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_abs_Fa$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_fp_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_gene_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_dist$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_abs_Al$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ripple3d$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_raytrace_initialize$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_raytrace$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad_rad_transp$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ECRad$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o
	
ECRad_pythonOBJ = ECRad_python$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o

OBJS := $(addprefix $(MODECRad)/, $(OBJECTS))

# Rules
.PHONY: INFO
ifeq ($(IDA),True)
all: INFO directories lib
else
all: INFO directories lib \
	$(ECRadLIBDir)/ECRad_python$(OMPFLAG)$(USE3DFLAG)$(DB) $(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB)
endif

lib: directories MANIFEST \
	$(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a \
	$(ECRadLIBDir)/ECRad_python$(OMPFLAG)$(USE3DFLAG)$(DB)$(F2PYEXT_SUFFIX)
	
MANIFEST: 
	echo include src/ecrad_core/ECRad_python$(OMPFLAG)$(USE3DFLAG)$(DB)$(F2PYEXT_SUFFIX) >> MANIFEST.in
ifeq ($(COMPILER),GNU)
INFO:
	echo "Assuming GNU toochain"
else
INFO:
	echo "Assuming INTEL toochain"
endif

$(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB): $(OBJJ)  \
	$(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a $(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o makefile
	echo $(OBJJ$(SYS))
	$(F90) $(LDFLAGS) $(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o ${OBJJ} $(LIBS) \
	-o $(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB)
	
$(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o : $(SRCP)/$(APPLICATION).f90
	$(F90) ${MODULES} $(FFPFLAGS) -c $(F90FLAGS) $< -o $@
	
$(ECRadLIBDir)/ECRad_python$(OMPFLAG)$(USE3DFLAG)$(DB)$(F2PYEXT_SUFFIX): $(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a
	cd $(ECRad_pythonDir); \
	python -m numpy.f2py $(F2PYDBG) -c --fcompiler=$(F2PYCOMPILER) $(ROOTDIR)/src/ECRad_python$(OMPFLAG)$(USE3DFLAG).f90 -m ECRad_python$(OMPFLAG)$(USE3DFLAG)$(DB) \
		-I$(MODECRad) --opt='' --f90flags='$(F2PYFLAGS)' $(F2PYLIBS); \
	cd -

#libECRad
$(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a: $(OBJS)
	ar rv $@ $(OBJS)

$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a: $(OBJS)

$(MODECRad)/%$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(SRCP)/%.f90
	 $(F90) $(MODULES) $(FFPFLAGS) -c $(F90FLAGS) $< -o $@
	 
#making the directories
directories:
	${MKDIR_P} ${ECRadLIBDir}
	${MKDIR_P} ${MODECRad} 
	${MKDIR_P} $(ECRad_pythonDir)

#Dependencies

$(MODECRad)/mod_contour$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB)

$(MODECRad)/mod_ECRad_types$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90

$(MODECRad)/mod_ECRad_interpol$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_types.f90

$(MODECRad)/mod_ECRad_abs_Fa$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB)

$(MODECRad)/mod_ECRad_fp_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODECRad)/mod_ECRad_gene_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_interpol.f90

$(MODECRad)/mod_ECRad_dist$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_fp_dist_utils.f90 \
	$(SRCP)/mod_ECRad_gene_dist_utils.f90 \
  $(SRCP)/mod_ECRad_interpol.f90

$(MODECRad)/mod_ECRad_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_interpol.f90 \
	$(SRCP)/mod_ECRad_fp_dist_utils.f90 \
	$(SRCP)/mod_ECRad_gene_dist_utils.f90

$(MODECRad)/mod_ECRad_abs_Al$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_utils.f90 \
	$(SRCP)/mod_ECRad_fp_dist_utils.f90 \
	$(SRCP)/mod_ECRad_gene_dist_utils.f90 \
	$(SRCP)/mod_ECRad_dist.f90 \
	$(SRCP)/mod_ECRad_abs_Fa.f90

$(MODECRad)/mod_ripple3d$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_types.f90 

$(MODECRad)/mod_ECRad_raytrace_initialize$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ECRad_interpol.f90 \
	$(SRCP)/mod_ECRad_utils.f90

$(MODECRad)/mod_ECRad_raytrace$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ECRad_raytrace_initialize.f90 \
	$(SRCP)/mod_ECRad_interpol.f90 \
	$(SRCP)/mod_ECRad_utils.f90

$(MODECRad)/mod_ECRad_rad_transp$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_utils.f90 \
	$(SRCP)/mod_ECRad_abs_Al.f90 \
	$(SRCP)/mod_ECRad_raytrace.f90

$(MODECRad)/mod_ECRad$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: \
	$(SRCP)/mod_ECRad_types.f90 \
	$(SRCP)/mod_ECRad_rad_transp.f90 \
	$(SRCP)/mod_ECRad_abs_Al.f90 \
	$(SRCP)/mod_ECRad_raytrace_initialize.f90 \
	$(SRCP)/mod_ECRad_raytrace.f90 \
	$(SRCP)/mod_ECRad_interpol.f90 \
	$(SRCP)/mod_ECRad_utils.f90

$(MODECRad)/ECRad_python$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: \
	$(SRCP)/mod_ECRad.f90

ifndef PREFIX
clean:
	ifneq ($(ROOTDIR),$(ECRadLIBDir))
		rm -rf $(ECRadLIBDir)
	endif
else
clean:
	rm -rf $(ECRadLIBDir)/*ECRad*
	rm -rf $(ECRadLIBDir)/*ecrad*
	rm -rf 
	rm -rf $(ROOTDIR)/MANIFEST.in
endif