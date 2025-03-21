FC = gfortran
FCFLAGS = -O3 -Wall -Wextra
FLFLAGS =

PROG = vdw_gas
SRC = ${PROG:=.f90}
MOD = subroutines.f90
OBJ = ${SRC:.f90=.o} ${MOD:.f90=.o}

.SUFFIXES: .f90 .o
.PHONY: all clean plot
.INTERMEDIATE: .o .mod

all: ${PROG} binning jackknife

modules: ${MOD}
	$(FC) $(FCFLAGS) -c $<

.f90.o:
	$(FC) $(FCFLAGS) -c $<

$(PROG): ${OBJ}
	$(FC) $(FCFLAGS) -o $@ ${@:=.o} ${MOD:.f90=.o}

binning: binning.o
	$(FC) $(FCFLAGS) -o $@ ${@:=.o}

jackknife: jackknife.o
	$(FC) $(FCFLAGS) -o $@ ${@:=.o}

${OBJ}: modules

plot:
	sh ./plot/plot.sh

clean:
	rm ${PROG} *.o binning jackknife *.mod

clean-output:
	rm *.xyz *.dat
