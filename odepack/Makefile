# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
LIB=ODEPack
ifdef PREFIX
	LIBDIR = $(PREFIX)/lib
else
	LIBDIR = $(CURDIR)
endif
ifeq ($(FC),gfortran)
	FFLAGS = -O2 -mavx -fopenmp -fpic
else
	FFLAGS = -fpic -O2 -fpic -fopenmp 
endif
OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
ifdef PREFIX
all: directory $(LIBDIR)/lib$(LIB).a	
else
all: $(LIBDIR)/lib$(LIB).a
endif

directory:
	mkdir -p $(LIBDIR)

$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
$(LIBDIR)/lib$(LIB).a: $(OBJ)
	$(AR) rus $(LIBDIR)/lib$(LIB).a $?
clean:
	rm -f *.o
	rm -f $(LIBDIR)/lib$(LIB).a