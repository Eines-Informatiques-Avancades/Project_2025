FC = gfortran
FCFLAGS = -O3 -Wall -Wextra
FLFLAGS =

PROG = vdw_gas
SRC = ${PROG:=.f90}
MOD = subroutines.f90
OBJ = ${SRC:.f90=.o} ${MOD:.f90=.o}

.SUFFIXES: .f90 .o
.PHONY: all clean
.INTERMEDIATE: .o .mod

all: ${PROG}

modules: ${MOD}
	$(FC) $(FCFLAGS) -c $<

.f90.o:
	$(FC) $(FCFLAGS) -c $<

$(PROG): ${OBJ}
	$(FC) $(FCFLAGS) -o $@ ${@:=.o} ${MOD:.f90=.o}

${OBJ}: modules

clean:
	rm ${PROG} ${OBJ} *.mod

clean-output:
	rm *.xyz *.dat
