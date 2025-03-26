.SUFFIXES:

FC = gfortran
FCFLAGS = -O3 -Wall -Wextra -J$(MODDIR)
FLFLAGS =
MODDIR = ./include

PROG = vdw_gas
SRC = ${PROG:=.f90}
MOD = ${MODDIR}/geometry.f90 ${MODDIR}/lj_forces.f90 ${MODDIR}/initial_conf.f90 \
	  ${MODDIR}/integrators.f90 ${MODDIR}/io.f90 ${MODDIR}/post_trajectory_analysis.f90 \
	  ${MODDIR}/thermodynamics.f90 ${MODDIR}/thermostat.f90 ${MODDIR}/global_vars.f90 \
          ${MODDIR}/lj_potentials.f90
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

${MODDIR}/lj_forces.anc: ${MODDIR}/global_vars.anc ${MODDIR}/geometry.anc ${MODDIR}/lj_forces.mod
${MODDIR}/lj_potentials.anc: ${MODDIR}/global_vars.anc ${MODDIR}/geometry.anc ${MODDIR}/lj_potentials.mod
${MODDIR}/geometry.anc: ${MODDIR}/global_vars.anc ${MODDIR}/geometry.mod
${MODDIR}/global_vars.anc: ${MODDIR}/global_vars.mod
${MODDIR}/initial_conf.anc:${MODDIR}/global_vars.anc ${MODDIR}/initial_conf.mod
${MODDIR}/integrators.anc: ${MODDIR}/global_vars.anc ${MODDIR}/lj_forces.anc ${MODDIR}/geometry.anc ${MODDIR}/integrators.mod
${MODDIR}/io.anc: ${MODDIR}/global_vars.anc ${MODDIR}/io.mod
${MODDIR}/post_trajectory_analysis.anc: ${MODDIR}/global_vars.anc ${MODDIR}/post_trajectory_analysis.mod
${MODDIR}/thermodynamics.anc: ${MODDIR}/global_vars.anc ${MODDIR}/thermodynamics.mod
${MODDIR}/thermostat.anc: ${MODDIR}/global_vars.anc ${MODDIR}/thermostat.mod
vdw_gas.anc: ${MODDIR}/lj_forces.anc ${MODDIR}/geometry.anc ${MODDIR}/initial_conf.anc \
	${MODDIR}/integrators.anc ${MODDIR}/io.anc ${MODDIR}/post_trajectory_analysis.anc \
	${MODDIR}/thermodynamics.anc ${MODDIR}/thermostat.anc ${MODDIR}/global_vars.anc \
        ${MODDIR}/lj_potentials.anc

${MODDIR}/lj_potentials.mod:
${MODDIR}/lj_forces.mod:
${MODDIR}/geometry.mod:
${MODDIR}/global_vars.mod:
${MODDIR}/initial_conf.mod:
${MODDIR}/integrators.mod:
${MODDIR}/io.mod:
${MODDIR}/post_trajectory_analysis.mod:
${MODDIR}/thermodynamics.mod:
${MODDIR}/thermostat.mod:

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
