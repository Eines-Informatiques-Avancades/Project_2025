.SUFFIXES:

FC = gfortran
FCFLAGS = -O3 -Wall -Wextra -J$(MODDIR)
FLFLAGS =
MODDIR = ./include

PROG = vdw_gas
SRC = ${PROG:=.f90}
MOD = ${MODDIR}/geometry.f90 ${MODDIR}/lj_forces.f90 ${MODDIR}/initial_conf.f90 \
	  ${MODDIR}/integrators.f90 ${MODDIR}/io.f90 ${MODDIR}/post_trajectory_analysis.f90 \
	  ${MODDIR}/thermodynamics.f90 ${MODDIR}/thermostat.f90
OBJ = ${MOD:.f90=.o} ${SRC:.f90=.o}
ANC = ${OBJ:.o=.anc}

.PHONY: all clean plot

all: ${PROG} binning jackknife

# Main program compilation

$(PROG): ${OBJ}
	$(FC) $(FCFLAGS) -o $@ ${@:=.o} ${MOD:.f90=.o}

%.anc: %.f90
	$(FC) $(FCFLAGS) -fsyntax-only -c -o $@ $<
	@touch $@

%.o : %.anc
	$(FC) $(FCFLAGS) -c -o $@ $(<:.anc=.f90)
	@touch $@

include/lj_forces.anc: include/geometry.anc include/lj_forces.mod
include/geometry.anc: include/geometry.mod
include/initial_conf.anc: include/initial_conf.mod
include/integrators.anc: include/lj_forces.anc include/geometry.anc include/integrators.mod
include/io.anc: include/io.mod
include/post_trajectory_analysis.anc: include/post_trajectory_analysis.mod
include/thermodynamics.anc: include/thermodynamics.mod
include/thermostat.anc: include/thermostat.mod
vdw_gas.anc: include/lj_forces.anc include/geometry.anc include/initial_conf.anc \
	include/integrators.anc include/io.anc include/post_trajectory_analysis.anc \
	include/thermodynamics.anc include/thermostat.anc

include/lj_forces.mod:
include/geometry.mod:
include/initial_conf.mod:
include/integrators.mod:
include/io.mod:
include/post_trajectory_analysis.mod:
include/thermodynamics.mod:
include/thermostat.mod:

# Statistical analysis

binning: binning.o
	$(FC) $(FCFLAGS) -o $@ ${@:=.o}

jackknife: jackknife.o
	$(FC) $(FCFLAGS) -o $@ ${@:=.o}

plot:
	sh ./plot/plot.sh

clean:
	rm -f ${PROG} binning jackknife ${OBJ} ${OBJ:.o=.mod} ${ANC} *.o

clean-output:
	rm *.xyz *.dat
