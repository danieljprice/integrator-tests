#----------------------------------
#
#  Makefile for simple n-body code
#  Daniel Price, May 2016
#
#----------------------------------
FC=gfortran
FFLAGS=-O3 -Wall -Wextra -pedantic -std=f2008

SRC=config.f90 force.f90 step.f90 nbody.f90
OBJ=${SRC:.f90=.o}

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

nbody: ${OBJ}
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f *.out *.ev nbody
