FC = gfortran
FCFLAGS = -O3 -Wall -Wextra
FLFLAGS =

PROG = vdw_gas
SRC = vdw_gas.f90
MOD = subroutines.f90
OBJ = ${SRC:.f90=.o}

.SUFFIXES: .f90 .o
.PHONY: all clean
.INTERMEDIATE: .o .mod

all: vdw_gas

modules: ${MOD}
	$(FC) $(FCFLAGS) -fsyntax-only -c $<

.f90.o:
	$(FC) $(FCFLAGS) -c $<

$(PROG): ${OBJ}
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm ${PROG} ${OBJ} *.mod

vdw_gas.o: modules
