COMPILER=i
APPLICATION = ECRad

ROOTDIR=$(CURDIR)
ECRadLIBDir=$(ROOTDIR)/$(SYS)
ECRadLIB=$(APPLICATION)
SRCP=$(ROOTDIR)/src
ifeq ($(IDA),True)
	ifneq ($(USE_3D),True)
		STDP=/afs/ipp/cips/ipp/data_analysis/lib_std/$(SYS)
		MODSTDP=$(STDP)/mod$(COMPILER)
		STDPLIB=$(STDP)/libstd$(COMPILER).a
	else
		STDPLIB = $(SRCP)/std_lib.f90
	endif
	IDAFLAG = IDA
else
	STDPLIB = $(SRCP)/std_lib.f90
endif
obj=
#Linux
F90 = ifort
F77 = ifort
MKDIR_P = mkdir -p
FFPFLAGS = -fpp
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
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)
MODECRad=$(ROOTDIR)/$(SYS)/mod$(COMPILER)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)
$(shell   mkdir -p $(MODECRad))
# Debugging -> enable if DEBUG==		True
ifeq ($(DEBUG),True)
	F90FLAGS = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check-uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
	DB = db
	# Optimized
else
  F90FLAGS = -c -O2 -fp-model source -axavx
  # Profiling
  #F90FLAGS = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
endif
ifeq ($(OPEN_MP),True)
  # Parallelisation
  F90FLAGS += -qopenmp
endif
# Debugging parallelisation
# F90FLAGS = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
F77FLAGS = $(F90FLAGS) -C
# Libraries
FITPACK = $(ROOTDIR)/../netlib/fitpack/lib.a
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(OMPFLAG)$(USE3DFLAG)$(DB) $(NAGF90LIB) $(NAGFLIB) $(FITPACK) $(MAGCONF)
ifeq ($(USE_3D),True)
 	LIBS += $(ROOTDIR)/../Mconf/lib/libmconf64.a
#   LIBS += $(ROOTDIR)/../magconf/lib/libmconf64.a
#../Mconf/unix/mconf_matlab64.so
	LIBS += -lpthread -lstdc++
	LIBS += -L${NETCDF_HOME}/lib/  -lnetcdf_c++4  -lnetcdf 
endif
ifeq ($(OPEN_MP),True)
	LIBS += -qopenmp
endif
LDFLAGS = -z muldefs
MODULES = -module $(MODECRad)
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
	mod_ecfm_refr_types$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_interpol$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_abs_Fa$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_dist$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_abs_Al$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ripple3d$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	dlsode$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_raytrace$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_rad_transp$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o #\
	Fortran_Stop_Handler$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o \

OBJS := $(addprefix $(MODECRad)/, $(OBJECTS))

# Rules
#ECRad application
# directories
ifeq ($(IDA),True)
all: $(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a
else
all: $(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a $(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB)
endif

$(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB): $(OBJJ) $(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o makefile
	echo $(OBJJ$(SYS))
	$(F90) $(LDFLAGS) $(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o ${OBJJ} $(LIBS) \
	-o $(ECRadLIBDir)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB)
	
$(MODECRad)/$(APPLICATION)$(OMPFLAG)$(USE3DFLAG)$(DB).o : $(SRCP)/$(APPLICATION).f90
	$(F90) ${MODULES} $(FFPFLAGS) $(F90FLAGS) $< -o $@

#libECRad
$(ECRadLIBDir)/lib$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a: makefile $(OBJS) $(SOURCE)
	ar rv $@ $(OBJS)
	ar -t $@ -o $@
$(ECRadLIB)$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).a: $(OBJS)

$(MODECRad)/%$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(SRCP)/%.f90
	 $(F90) ${MODULES} $(FFPFLAGS) $(F90FLAGS) $< -o $@
	 
$(MODECRad)/dlsode$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(SRCP)/dlsode.f
	 $(F77) $(FFPFLAGS) $(F77FLAGS) \
	 $< -o $@
# making the directories
#directories: ${MODECRad}

#${ECRadLIBDir}:
#	${MKDIR_P} ${ECRadLIBDir}

#${MODECRad}:
#	${MKDIR_P} ${MODECRad}

#Dependencies

$(MODECRad)/mod_contour$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_types$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90

$(MODECRad)/mod_ecfm_refr_interpol$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90

$(MODECRad)/mod_ecfm_refr_abs_Fa$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODECRad)/mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_dist$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
  $(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_utils$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODECRad)/mod_ecfm_refr_abs_Al$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Fa.f90

$(MODECRad)/mod_ripple3d$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 

$(MODECRad)/mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_raytrace$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/dlsode.f \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_rad_transp$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90

$(MODECRad)/mod_ecfm_refr$(IDAFLAG)$(OMPFLAG)$(USE3DFLAG)$(DB).o: \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_rad_transp.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

clean:
	rm -r $(ECRadLIBDir)