# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
LIB=ODEPack
FC = ifort
FFLAGS = -O2 -fp-model source -axavx -fpic -qopenmp
OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).a
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
clean:
	rm *.o






