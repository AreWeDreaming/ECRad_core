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
else
	ifneq ($(DEBUG),True)
		FFPFLAGS += -DOMP
	endif
endif
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)


# Debugging -> enable if DEBUG==		True
ifeq ($(DEBUG),True)
	F90FLAGS = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
	DB = db
else ifeq ($(IDA),True)
	# Optimized
  F90FLAGS = -c -O2 -fp-model source -axavx
  # Profiling
  #F90FLAGS = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
else
  # Parallelisation
  F90FLAGS = -c -O2 -fp-model source -qopenmp -axavx 
endif	
# Debugging parallelisation
# F90FLAGS = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
F77FLAGS = $(F90FLAGS) -C
# Libraries
FITPACK = $(ROOTDIR)/../netlib/fitpack/lib.a
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(DB) $(NAGF90LIB) $(NAGFLIB) $(FITPACK)
ifneq ($(DEBUG),True)
	ifneq ($(IDA),True)
		LIBS += -qopenmp
	endif
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
OBJECTS =	std_lib$(IDAFLAG)$(DB).o
endif
OBJECTS += \
	quadrature$(IDAFLAG)$(DB).o \
	mod_contour$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_types$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_interpol$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_abs_Fa$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_utils$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_dist$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_abs_Al$(IDAFLAG)$(DB).o \
	mod_ripple3d$(IDAFLAG)$(DB).o \
  mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_em_Hu$(IDAFLAG)$(DB).o \
	dlsode$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_raytrace$(IDAFLAG)$(DB).o \
	mod_ecfm_refr_rad_transp$(IDAFLAG)$(DB).o \
	mod_ecfm_refr$(IDAFLAG)$(DB).o 

OBJS := $(addprefix $(MODECRad)/, $(OBJECTS))

# Rules
#ECRad application
ifeq ($(IDA),True)
all: directories lib$(ECRadLIB)$(IDAFLAG)$(DB).a
else
all: directories lib$(ECRadLIB)$(IDAFLAG)$(DB).a $(APPLICATION)$(DB)
endif

$(APPLICATION)$(DB): $(OBJJ) $(APPLICATION)$(DB).o makefile
	echo $(OBJJ$(SYS))
	cd $(ECRadLIBDir) && $(F90) $(LDFLAGS) $(APPLICATION)$(DB).o ${OBJJ} $(LIBS) \
	-o $(APPLICATION)$(DB)
	
$(APPLICATION)$(DB).o : $(SRCP)/$(APPLICATION).f90
	cd $(ECRadLIBDir) && $(F90) ${MODULES} $(FFPFLAGS) -o $(APPLICATION)$(DB).o \
	$(F90FLAGS) $(SRCP)/$(APPLICATION).f90

#libECRad
lib$(ECRadLIB)$(IDAFLAG)$(DB).a: makefile $(OBJS) $(SOURCE)
	echo $<
	cd $(MODECRad)&&ar rv $@ $(OBJS)
	mv $(MODECRad)/$@ $(ECRadLIBDir)
	ar -t $(ECRadLIBDir)/$@
$(ECRadLIB)$(IDAFLAG)$(DB).a: $(OBJS)

$(MODECRad)/%$(IDAFLAG)$(DB).o: $(SRCP)/%.f90
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
	$(MODECRad)/std_lib$(IDAFLAG)$(DB).o:   $(SRCP)/std_lib.f90
endif
$(MODECRad)/quadrature$(IDAFLAG)$(DB).o:   $(SRCP)/quadrature.f90 

$(MODECRad)/mod_contour$(IDAFLAG)$(DB).o:   $(SRCP)/mod_contour.f90 $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_types$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_types.f90 $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_interpol$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_interpol.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90

$(MODECRad)/mod_ecfm_refr_abs_Fa$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_abs_Fa.f90 $(STDPLIB)

$(MODECRad)/mod_ecfm_refr_fp_dist_utils$(IDAFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODECRad)/mod_ecfm_refr_gene_dist_utils$(IDAFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_dist$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_dist.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
  $(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/mod_ecfm_refr_utils$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_utils.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODECRad)/mod_ecfm_refr_abs_Al$(IDAFLAG)$(DB).o:		$(SRCP)/mod_ecfm_refr_abs_Al.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Fa.f90

$(MODECRad)/mod_ripple3d$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ripple3d.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 

$(MODECRad)/mod_ecfm_refr_raytrace_initialize$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_em_Hu$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_em_Hu.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90

$(MODECRad)/dlsode$(IDAFLAG)$(DB).o:   $(SRCP)/dlsode.f
	 cd $(MODECRad)&& $(F77) $(FFPFLAGS) $(F77FLAGS) -o dlsode$(IDAFLAG)$(DB).o \
	 $(SRCP)/dlsode.f

$(MODECRad)/mod_ecfm_refr_raytrace$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/dlsode.f \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_rad_transp$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr_rad_transp.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_em_Hu.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90

$(MODECRad)/mod_ecfm_refr$(IDAFLAG)$(DB).o:   $(SRCP)/mod_ecfm_refr.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_em_Hu.f90 \
	$(SRCP)/mod_ecfm_refr_rad_transp.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

clean:
	rm -r $(ECRadLIBDir)