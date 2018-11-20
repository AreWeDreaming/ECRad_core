# setenv SYS sun4x_59

# lAe 09.11.05, /Fi
# .KEEP_STATE:
# SUFFIXES: .o .f
#set -x
.SUFFIXES: .o  .a    .f90

COMPILER=i


ROOTDIR=$(CURDIR)/..
LIPLIB=$(ROOTDIR)/$(SYS)

SRCP=$(ROOTDIR)/src
MODSTDP=$(SRCP)/std_lib.f90
#Linux
LINUX = $(SYS)
#CC$(LINUX) = /usr/bin/gcc
F90$(LINUX) = ifort
F77$(LINUX) = ifort
FPP_DEFSi$(LINUX) = -fpp #-DAUG # -DNAG #-DRCV
#CFLAGS$(LINUX) = -c -Wall -O -I/afs/ipp/aug/ads/i386_linux24/include -I/afs/ipp/aug/ads/i386_linux24/codes
#F90_FLAGf$(LINUX) = -c -O3 -r8
#F90_FLAGi$(LINUX) = -c -O
#F90_FLAGi$(LINUX) =  -c -O2 -fp-model source -axavx
#F90_FLAGi$(LINUX) =  -c -O2 -fp-model source -openmp-report -openmp -axavx 
# For profiling
#F90_FLAGi$(LINUX) = -c -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir
# -openmp -openmp-report2 -parallel
# For debugging
F90_FLAGi$(LINUX) = -c -g -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -debug all -gen-interfaces -warn interfaces -fpe3
# Debug of parallelisation
#F90_FLAGi$(LINUX) = -c -g -O0 -shared-libgcc -traceback -check bounds  -check all -u -warn all -diag-disable 7712 -check uninit -fp-model source -gen-interfaces -warn interfaces -fpe3 -openmp -openmp-report -DTBB_USE_DEBUG
#F90_FLAGi$(LINUX) = -c -O2 -fp-model source
F77_FLAG$(LINUX) = $(F90_FLAG$(COMPILER)$(SYS)) -w
#F77_FLAG$(LINUX) = -c -g -w
LILIBL$(LINUX) = li_lib$(COMPILER)
LILIB$(LINUX) = $(LIPLIB)/lib$(LILIBL$(LINUX)).a
	# ECE Forwardmodeling Library
#ECFMP		= /afs/ipp/home/s/skr/F90/ECE/lib_ecfm/$(SYS)
#MODECFM	= $(ECFMP)/mod
#LIBP		= /afs/ipp/home/s/skr/F90/LIB/$(SYS)
#LIB			= $(LIBP)/lib.a
#LIBMOD 		= $(LIBP)/mod
#MODULES$(LINUX) = -I$(LIPLIB)/mod$(COMPILER) -I$(MODSTDP) -I$(LIBMOD) -I$(MODECFM) $(NAGF90MOD)
MODULES$(LINUX) = -I$(LIPLIB)/mod$(COMPILER) -I$(MODSTDP) $(NAGF90MOD)
#MODULES$(LINUX) = -I$(MODSTDP) $(NAGF90MOD) -I$(LIPLIB)/mod$(COMPILER)
MODSRCP$(LINUX) = $(LIPLIB)/mod$(COMPILER)
MAKEFILE$(LINUX) = makefile$(COMPILER)
LESTAX$(LINUX) = lestax_rh90
LESTAX_OPO$(LINUX) = lestax_opo_rh90
MATRIX$(LINUX) = matrix_rh90

CC = $(CC$(SYS))
F90 = $(F90$(SYS))
F77C = $(F77$(SYS))
CFLAGS += $(CFLAGS$(SYS))
FFLAGS += -c
FFP_FLAGS=$(FPP_DEFS$(COMPILER)$(SYS))
F77FLAGS=$(F77_FLAG$(SYS))
F90FLAGS=$(F90_FLAG$(COMPILER)$(SYS))
LIBPPLOT =  $(LIBPPLOT$(SYS))
MODULES = $(MODULES$(SYS))
LILIB = $(LILIB$(SYS))
LIBLI = -L$(LIPLIB) -l$(LILIBL$(SYS))
MODSRCP = $(MODSRCP$(SYS))
MAKEFILE = $(MAKEFILE$(SYS))
LESTAX = $(LESTAX$(SYS))
LESTAX_OPO = $(LESTAX_OPO$(SYS))
MATRIX = $(MATRIX$(SYS))

OBJECTS = \
	std_lib.o \
	quadrature.o \
	mod_contour.o \
	ida_type.o \
	ece_type.o \
	mod_fit_params.o \
	mod_ecfm_refr_types.o \
	mod_ecfm_refr_interpol.o \
	mod_ecfm_refr_abs_Fa.o \
	mod_ecfm_refr_fp_dist_utils.o \
	mod_ecfm_refr_gene_dist_utils.o \
	libe_utils.o \
	eqi_type.o \
	mod_eqi.o \
	mod_ecfm_refr_utils.o \
	mod_ecfm_refr_dist.o \
	mod_ecfm_refr_abs_Al.o \
	mod_ripple3d.o \
  mod_ecfm_refr_raytrace_initialize.o \
	mod_ecfm_refr_em_Hu.o \
	dlsode.o \
	mod_ecfm_refr_raytrace.o \
	mod_ecfm_refr_rad_transp.o \
	mod_ecfm_refr.o 

OBJS := $(addprefix $(MODSRCP)/, $(OBJECTS))

#Rules: libli_libi
lib$(LILIBL$(SYS)).a: makefile $(OBJS) $(SOURCE)
	echo $<
	cd $(MODSRCP)&&ar rv $@ $(OBJS)
	mv $(MODSRCP)/$@ .
	ar -t $@
lib$(LILIBL$(SYS)).a: $(OBJS)

$(MODSRCP)/%.o: $(SRCP)/%.f90
	 cd $(MODSRCP)&& $(F90) ${MODULES} $(FFP_FLAGS) $(F90FLAGS) $<

#Dependencies

$(MODSRCP)/quadrature.o:   $(SRCP)/quadrature.f90

$(MODSRCP)/mod_contour.o:   $(SRCP)/mod_contour.f90 $(STDPLIB)

$(MODSRCP)/ida_type.o:   $(SRCP)/ida_type.f90 $(STDPLIB)

$(MODSRCP)/ece_type.o:   $(SRCP)/ece_type.f90 $(STDPLIB)

$(MODSRCP)/libe_utils.o:   $(SRCP)/libe_utils.f90 $(STDPLIB)

$(MODSRCP)/mod_fit_params.o:   $(SRCP)/mod_fit_params.f90 $(STDPLIB) \
	$(SRCP)/ida_type.f90

$(MODSRCP)/mod_ecfm_refr_types.o:   $(SRCP)/mod_ecfm_refr_types.f90 $(STDPLIB) \
	$(SRCP)/ida_type.f90 \
	$(SRCP)/ece_type.f90

$(MODSRCP)/mod_ecfm_refr_interpol.o:   $(SRCP)/mod_ecfm_refr_interpol.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90

$(MODSRCP)/mod_ecfm_refr_abs_Fa.o:   $(SRCP)/mod_ecfm_refr_abs_Fa.f90 $(STDPLIB)

$(MODSRCP)/mod_ecfm_refr_fp_dist_utils.o:		$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_contour.f90

$(MODSRCP)/mod_ecfm_refr_gene_dist_utils.o:		$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90

$(MODSRCP)/eqi_type.o:   $(SRCP)/eqi_type.f90 $(STDPLIB)

$(MODSRCP)/mod_eqi.o:   $(SRCP)/mod_eqi.f90 $(STDPLIB) \
	$(SRCP)/eqi_type.f90

$(MODSRCP)/mod_ecfm_refr_dist.o:   $(SRCP)/mod_ecfm_refr_dist.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODSRCP)/mod_ecfm_refr_utils.o:   $(SRCP)/mod_ecfm_refr_utils.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_eqi.f90 \
	$(SRCP)/libe_utils.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90

$(MODSRCP)/mod_ecfm_refr_abs_Al.o:		$(SRCP)/mod_ecfm_refr_abs_Al.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Fa.f90

$(MODSRCP)/mod_ripple3d.o:   $(SRCP)/mod_ripple3d.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 

$(MODSRCP)/mod_ecfm_refr_raytrace_initialize.o:   $(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODSRCP)/mod_ecfm_refr_em_Hu.o:   $(SRCP)/mod_ecfm_refr_em_Hu.f90 $(STDPLIB) \
	$(SRCP)/quadrature.f90 \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_fp_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_gene_dist_utils.f90 \
	$(SRCP)/mod_ecfm_refr_dist.f90

$(MODSRCP)/dlsode.o:   $(SRCP)/dlsode.f
	 cd $(MODSRCP)&& $(F77C) $(FFLAGS) $(F77FLAGS) $(SRCP)/dlsode.f

$(MODSRCP)/mod_ecfm_refr_raytrace.o:   $(SRCP)/mod_ecfm_refr_raytrace.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/dlsode.f \
	$(SRCP)/mod_ripple3d.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace_initialize.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90

$(MODSRCP)/mod_ecfm_refr_rad_transp.o:   $(SRCP)/mod_ecfm_refr_rad_transp.f90 $(STDPLIB) \
	$(SRCP)/mod_ecfm_refr_types.f90 \
	$(SRCP)/mod_ecfm_refr_interpol.f90 \
	$(SRCP)/mod_ecfm_refr_utils.f90 \
	$(SRCP)/mod_ecfm_refr_em_Hu.f90 \
	$(SRCP)/mod_ecfm_refr_abs_Al.f90 \
	$(SRCP)/mod_ecfm_refr_raytrace.f90

$(MODSRCP)/mod_ecfm_refr.o:   $(SRCP)/mod_ecfm_refr.f90 $(STDPLIB) \
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
	 rm -f modi/*

