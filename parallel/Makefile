.SUFFIXES:

FC = mpif90
FCFLAGS = -O3 -module $(MODDIR) # ifort flags
FLFLAGS =
MODDIR = ./include
CLUSTERDIR = ./cluster

PROG = vdw_gas
SRC = ${PROG:=.f90}
MOD = ${MODDIR}/geometry.f90 ${MODDIR}/lj_forces.f90 ${MODDIR}/initial_conf.f90 \
	  ${MODDIR}/integrators.f90 ${MODDIR}/io.f90 ${MODDIR}/post_trajectory_analysis.f90 \
	  ${MODDIR}/thermodynamics.f90 ${MODDIR}/thermostat.f90 ${MODDIR}/global_vars.f90
OBJ = ${MOD:.f90=.o} ${SRC:.f90=.o}
ANC = ${OBJ:.o=.anc}

.PHONY: all clean plot

all: ${PROG} binning

# Cluster-specific commands

cluster-compile:
	qsub $(CLUSTERDIR)/submit_compile_cerqt01.sub

cluster-run-all: cluster-run-1 cluster-run-2 cluster-run-4 cluster-run-8 \
	cluster-run-16 cluster-run-32 cluster-run-40

cluster-run-1:
	rm -f vdw_gas_1_core.out
	qsub -N 'vdw_gas_1_core' -pe smp 1 \
		-o 'vdw_gas_1_core.out' -e 'vdw_gas_1_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-2:
	rm -f vdw_gas_2_core.out
	qsub -N 'vdw_gas_2_core' -pe smp 2 \
		-o 'vdw_gas_2_core.out' -e 'vdw_gas_2_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-4:
	rm -f vdw_gas_4_core.out
	qsub -N 'vdw_gas_4_core' -pe smp 4 \
		-o 'vdw_gas_4_core.out' -e 'vdw_gas_4_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-8:
	rm -f vdw_gas_8_core.out
	qsub -N 'vdw_gas_8_core' -pe smp 8 \
		-o 'vdw_gas_8_core.out' -e 'vdw_gas_8_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-16:
	rm -f vdw_gas_16_core.out
	qsub -N 'vdw_gas_16_core' -pe smp 16 \
		-o 'vdw_gas_16_core.out' -e 'vdw_gas_16_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-32:
	rm -f vdw_gas_32_core.out
	qsub -N 'vdw_gas_32_core' -pe smp 32 \
		-o 'vdw_gas_32_core.out' -e 'vdw_gas_32_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

cluster-run-40:
	rm -f vdw_gas_40_core.out
	qsub -N 'vdw_gas_40_core' -pe smp 40 \
		-o 'vdw_gas_40_core.out' -e 'vdw_gas_40_core.err' \
		$(CLUSTERDIR)/submit_run_cerqt01.sub

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
	${MODDIR}/thermodynamics.anc ${MODDIR}/thermostat.anc ${MODDIR}/global_vars.anc

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

plot:
	sh ./plot/plot.sh

clean:
	rm -f ${PROG} binning jackknife ${OBJ} ${OBJ:.o=.mod} ${ANC} *.o

clean-output:
	rm -f *.xyz *.dat
	rm -f output/*.xyz output/*.dat output/*.pdf
