# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
LIB=ODEPack
ifeq ($(COMPILER),g)
	FC = gfortran
else
	FC = ifort
endif
ifeq ($(FC),gfortran)
	FFLAGS = -O2 -mavx -fopenmp -fpic
else
	FFLAGS = -fpic -O2 -fpic -traceback -shared-intel -qopenmp 
endif
OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).a
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
clean:
	rm *.o