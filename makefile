COMPILER=i
APPLICATION = ECRad

ROOTDIR=$(CURDIR)
ECRadLIBDir=$(ROOTDIR)/$(SYS)
ECRadLIB=$(APPLICATION)
MODECRad=$(ROOTDIR)/$(SYS)/mod$(COMPILER)
SRCP=$(ROOTDIR)/src
obj=	
#Linux
F90 = ifort
F77 = ifort
MKDIR_P = mkdir -p
ifeq ($(DEBUG),True)
	FFPFLAGS = -fpp
else
	FFPFLAGS = -fpp -DOMP
endif
OBJJ = $(obj)
NAG_MOD = $(NAGF90MOD)
# Optimized
#F90FLAGS = -c -O2 -fp-model source -axavx
# Profiling
#F90FLAGS = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
# Parallelisation
F90FLAGS = -c -O2 -fp-model source -qopenmp -axavx 
#
# Debugging -> enable if DEBUG==		True
ifeq ($(DEBUG),True)
	F90FLAGS = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
	DB = db
endif	
# Debugging parallelisation
# F90FLAGS = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
F77FLAGS = $(F90FLAGS) -C
# Libraries
FITPACK = $(ROOTDIR)/../netlib/fitpack/lib.a
LIBS = -L$(ECRadLIBDir) -l$(ECRadLIB)$(DB) $(NAGF90LIB) $(NAGFLIB) $(FITPACK)
ifneq ($(DEBUG),True)
	$(LIBS) = $(LIBS) --openmp 
endif
LDFLAGS = -z muldefs
MODULES = $(NAG_MOD) -module $(MODECRad)

#Targets for library
OBJECTS = \
	std_lib$(DB).o \
	quadrature$(DB).o \
	mod_contour$(DB).o \
	ida_type$(DB).o \
	ece_type$(DB).o \
	mod_fit_params$(DB).o \
	mod_ecfm_refr_types$(DB).o \
	mod_ecfm_refr_interpol$(DB).o \
	mod_ecfm_refr_abs_Fa$(DB).o \
	mod_ecfm_refr_fp_dist_utils$(DB).o \
	mod_ecfm_refr_gene_dist_utils$(DB).o \
	libe_utils$(DB).o \
	eqi_type$(DB).o \
	mod_eqi$(DB).o \
	mod_ecfm_refr_utils$(DB).o \
	mod_ecfm_refr_dist$(DB).o \
	mod_ecfm_refr_abs_Al$(DB).o \
	mod_ripple3d$(DB).o \
  mod_ecfm_refr_raytrace_initialize$(DB).o \
	mod_ecfm_refr_em_Hu$(DB).o \
	dlsode$(DB).o \
	mod_ecfm_refr_raytrace$(DB).o \
	mod_ecfm_refr_rad_transp$(DB).o \
	mod_ecfm_refr$(DB).o 

OBJS := $(addprefix $(MODECRad)/, $(OBJECTS))

# Rules
#ECRad application

all: directories lib$(ECRadLIB)$(DB).a $(APPLICATION)$(DB)

$(APPLICATION)$(DB): $(OBJJ) $(APPLICATION)$(DB).o makefile
	echo $(OBJJ$(SYS))
	cd $(ECRadLIBDir) && $(F90) $(LDFLAGS) $(APPLICATION)$(DB).o ${OBJJ} $(LIBS) \
	-o $(APPLICATION)$(DB)
	
$(APPLICATION)$(DB).o : $(SRCP)/$(APPLICATION).f90
	cd $(ECRadLIBDir) && $(F90) ${MODULES} $(FFPFLAGS) -o $(APPLICATION)$(DB).o \
	$(F90FLAGS) $(SRCP)/$(APPLICATION).f90

#libECRad
lib$(ECRadLIB)$(DB).a: makefile $(OBJS) $(SOURCE)
	echo $<
	cd $(MODECRad)&&ar rv $@ $(OBJS)
	mv $(MODECRad)/$@ $(ECRadLIBDir)
	ar -t $(ECRadLIBDir)/$@
$(ECRadLIB)$(DB).a: $(OBJS)

$(MODECRad)/%$(DB).o: $(SRCP)/%.f90
	 cd $(MODECRad) && $(F90) ${MODULES} $(FFPFLAGS) -o $@ \
	 $(F90FLAGS) $<
# making the directories
directories: ${ECRadLIBDir} ${MODECRad}

${ECRadLIBDir}:
	${MKDIR_P} ${ECRadLIBDir}

${MODECRad}:
	${MKDIR_P} ${MODECRad}

#Dependencies

$(MODECRad)/std_lib$(DB).o:   $(SRCP)/std_lib.f90

$(MODECRad)/quadrature$(DB).o:   $(SRCP)/quadrature.f90

$(MODECRad)/mod_contour$(DB).o:   $(SRCP)/mod_contour.f90 $(SRCP)/std_lib.f90

$(MODECRad)/ida_type$(DB).o:   $(SRCP)/ida_type.f90 $(SRCP)/std_lib.f90

$(MODECRad)/ece_type$(DB).o:   $(SRCP)/ece_type.f90 $(SRCP)/std_lib.f90

$(MODECRad)/libe_utils$(DB).o:   $(SRCP)/libe_utils.f90 $(SRCP)/std_lib.f90

$(MODECRad)/mod_fit_params$(DB).o:   $(SRCP)/mod_fit_params.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/ida_type.f90

$(MODECRad)/mod_ecfm_refr_types$(DB).o:   $(SRCP)/mod_ecfm_refr_types.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/ida_type.f90 \
	$(SRCP)/ece_type.f90

$(MODECRad)/mod_ecfm_refr_interpol$(DB).o:   $(SRCP)/mod_ecfm_refr_interpol.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90

$(MODECRad)/mod_ecfm_refr_abs_Fa$(DB).o:   $(SRCP)/mod_ecfm_refr_abs_Fa.f90 $(SRCP)/std_lib.f90

$(MODECRad)/mod_ecfm_refr_fp_dist_utils$(DB).o:		$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODECRad)/mod_ecfm_refr_gene_dist_utils$(DB).o:		$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90

$(MODECRad)/eqi_type$(DB).o:   $(SRCP)/eqi_type.f90 $(SRCP)/std_lib.f90

$(MODECRad)/mod_eqi$(DB).o:   $(SRCP)/mod_eqi.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/eqi_type.f90

$(MODECRad)/mod_ecfm_refr_dist$(DB).o:   $(SRCP)/mod_ecfm_refr_dist.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODECRad)/mod_ecfm_refr_utils$(DB).o:   $(SRCP)/mod_ecfm_refr_utils.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_eqi.f90 \
	$(SRCP)/libe_utils.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODECRad)/mod_ecfm_refr_abs_Al$(DB).o:		$(SRCP)/mod_ecfm_refr_abs_Al.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Fa.f90

$(MODECRad)/mod_ripple3d$(DB).o:   $(SRCP)/mod_ripple3d.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 

$(MODECRad)/mod_ecfm_refr_raytrace_initialize$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_em_Hu$(DB).o:   $(SRCP)/mod_ecfm_refr_em_Hu.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90

$(MODECRad)/dlsode$(DB).o:   $(SRCP)/dlsode.f
	 cd $(MODECRad)&& $(F77) $(FFPFLAGS) $(F77FLAGS) -o dlsode$(DB).o \
	 $(SRCP)/dlsode.f

$(MODECRad)/mod_ecfm_refr_raytrace$(DB).o:   $(SRCP)/mod_ecfm_refr_raytrace.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/dlsode.f \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODECRad)/mod_ecfm_refr_rad_transp$(DB).o:   $(SRCP)/mod_ecfm_refr_rad_transp.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_em_Hu.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90

$(MODECRad)/mod_ecfm_refr$(DB).o:   $(SRCP)/mod_ecfm_refr.f90 $(SRCP)/std_lib.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_em_Hu.f90 \
	$(SRCP)/mod_ecfm_refr_rad_transp.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/ida_type.f90 \
	$(SRCP)/ece_type.f90 \
	$(SRCP)/mod_fit_params.f90

clean:
	rm -r $(ECRadLIBDir)