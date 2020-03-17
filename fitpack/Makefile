# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
LIB=fit
OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).a
$(OBJ):
	ifort -c $(FFLAGS) -O2 -fp-model source  -fpe3 -qopenmp -fpic -r8 $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
#ifort -c $(FFLAGS) -O2 -r8 -vec-report -g -prof-gen -prof-dir/afs/ipp-garching.mpg.de/home/s/sdenk/F90/Ecfm_Model_new/prof_dir $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
clean:
	rm *.o
	rm lib$(LIB).a






