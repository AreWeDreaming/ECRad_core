COMPILER=i
APPLICATION = ECRad

ROOTDIR=$(CURDIR)
ECRadLIBDir=$(ROOTDIR)/$(SYS)
ECRadLIB=$(APPLICATION)
MODECRad=$(ROOTDIR)/$(SYS)/mod$(COMPILER)
SRCP=$(ROOTDIR)/src
ifeq ($(IDA),True)
	STDP=/afs/ipp/cips/ipp/data_analysis/lib_std/$(SYS)
	MODSTDP=$(STDP)/mod$(COMPILER)
	STDPLIB=$(STDP)/libstd$(COMPILER).a
	IDAFLAG = IDA
	MODECRad=$(ROOTDIR)/$(SYS)/mod$(COMPILER)IDA
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
ifeq ($(OPENMP),True)
	FFPFLAGS += -DOMP
	OMPFLAG = OMP
endif
ifeq ($(USE_3D),True)
	FFPFLAGS += -DUSE_3D
	USE3DFLAG = USE3D
endif
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)


# Debugging -> enable if DEBUG==		True
ifeq ($(DEBUG),True)
	F90FLAGS = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
	DB = db
	# Optimized
  F90FLAGS = -c -O2 -fp-model source -axavx
  # Profiling
  #F90FLAGS = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
endif
ifeq ($(OPENMP),True)
  # Parallelisation
  F90FLAGS += -qopenmp
endif
# Debugging parallelisation
# F90FLAGS = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
F77FLAGS = $(F90FLAGS) -C
# Libraries
FITPACK = $(ROOTDIR)/../netlib/fitpack/lib.a
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(USE3DFLAG)$(DB) $(NAGF90LIB) $(NAGFLIB) $(FITPACK) $(MAGCONF)
ifeq ($(USE_3D),True)
	LIBS += $(ROOTDIR)/../Mconf/lib/libmconf64.a
#../Mconf/unix/mconf_matlab64.so
	LIBS += -lpthread -lstdc++
	LIBS += -L${NETCDF_HOME}/lib/  -lnetcdf_c++4  -lnetcdf 
endif
ifeq ($(OPENMP),True)
	LIBS += -qopenmp
endif
LDFLAGS = -z muldefs
MODULES = -module $(MODECRad)
ifeq ($(NAG),True)
	MODULES += $(NAG_MOD)
endif
ifeq ($(IDA),True)
	MODULES += -I$(MODSTDP)
endif

#Targets for library
ifeq ($(IDA),True)
OBJECTS =
else
OBJECTS =	std_lib$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o
endif
OBJECTS += \
	quadrature$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_contour$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	magconfig3D$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_types$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_interpol$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_abs_Fa$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_dist$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_abs_Al$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ripple3d$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	dlsode$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_raytrace$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr_rad_transp$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	mod_ecfm_refr$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o 

OBJS := $(addprefix $(MODECRad)/, $(OBJECTS))

# Rules
#ECRad application
ifeq ($(IDA),True)
all: directories lib$(ECRadLIB)$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).a
else
all: directories lib$(ECRadLIB)$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).a $(APPLICATION)$(USE3DFLAG)$(DB)
endif

$(APPLICATION)$(USE3DFLAG)$(DB): $(OBJJ) $(APPLICATION)$(USE3DFLAG)$(DB).o makefile
	echo $(OBJJ$(SYS))
	cd $(ECRadLIBDir) && $(F90) $(LDFLAGS) $(APPLICATION)$(USE3DFLAG)$(DB).o ${OBJJ} $(LIBS) \
	-o $(APPLICATION)$(USE3DFLAG)$(DB)
	
$(APPLICATION)$(OPENMP)$(USE3DFLAG)$(DB).o : $(SRCP)/$(APPLICATION).f90
	cd $(ECRadLIBDir) && $(F90) ${MODULES} $(FFPFLAGS) -o $(APPLICATION)$(OPENMP)$(USE3DFLAG)$(DB).o \
	$(F90FLAGS) $(SRCP)/$(APPLICATION).f90

#libECRad
lib$(ECRadLIB)$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).a: makefile $(OBJS) $(SOURCE)
	echo $<
	cd $(MODECRad)&&ar rv $@ $(OBJS)
	mv $(MODECRad)/$@ $(ECRadLIBDir)
	ar -t $(ECRadLIBDir)/$@
$(ECRadLIB)$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).a: $(OBJS)

$(MODECRad)/%$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o: $(SRCP)/%.f90
	 cd $(MODECRad) && echo $(MODULES) &&	$(F90) ${MODULES} $(FFPFLAGS) -o $@ \
	 $(F90FLAGS) $<
# making the directories
directories: ${ECRadLIBDir} ${MODECRad}

${ECRadLIBDir}:
	${MKDIR_P} ${ECRadLIBDir}

${MODECRad}:
	${MKDIR_P} ${MODECRad}

#Dependencies
ifneq ($(IDA),True)
	$(MODECRad)/std_lib$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/std_lib.f90
endif

$(MODECRad)/quadrature$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/quadrature.f90

$(MODECRad)/mod_contour$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_contour.f90 $(STDPLIB)

$(MODECRad)/magconfig3D$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/magconfig3D.f90

$(MODECRad)/mod_ecfm_refr_types$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_types.f90 $(STDPLIB) \
	$(SRCP)/magconfig3D.f90

$(MODECRad)/mod_ecfm_refr_interpol$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_interpol.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90

$(MODECRad)/mod_ecfm_refr_abs_Fa$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_abs_Fa.f90 $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODECRad)/mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_dist$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_dist.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
  $(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_utils$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_utils.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODECRad)/mod_ecfm_refr_abs_Al$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_abs_Al.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Fa.f90

$(MODECRad)/mod_ripple3d$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ripple3d.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 

$(MODECRad)/mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/dlsode$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/dlsode.f
	 cd $(MODECRad)&& $(F77) $(FFPFLAGS) $(F77FLAGS) -o dlsode$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o \
	 $(SRCP)/dlsode.f

$(MODECRad)/mod_ecfm_refr_raytrace$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace.f90 $(STDPLIB) \
	$(SRCP)/magconfig3D.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/dlsode.f \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_rad_transp$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_rad_transp.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90

$(MODECRad)/mod_ecfm_refr$(IDAFLAG)$(OPENMP)$(USE3DFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_rad_transp.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

clean:
	rm -r $(ECRadLIBDir)