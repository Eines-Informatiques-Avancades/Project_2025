#!/bin/sh
#
# run.sh
# Molecular Dynamics Simulation of a Van der Waals Gas
# Ricard Rodriguez
#

results_dir="./output"
log_dir="./log"
log_file="$(date +%Y%m%d_%H%M%S)_run.log"

# Print error message to stderr and exit the script.
err_exit() {
    printf >&2 "Error: %s\\n" "$@"
    exit 1
}

# Check if a given executable file exists and exit with an error message if it
# doesn't.
exec_exists() {
    [ -x "./$1" ] || err_exit "$1 does not exist or is not executable."
}

# Check that the necessary programs have been compiled.
exec_exists 'vdw_gas'
exec_exists 'binning'

cat > "$log_file" << EOF
# run.sh
# Molecular Dynamics Simulation of a Van der Waals Gas

Execution began at $(date '+%F %R')
EOF

./vdw_gas | tee -a "$log_file"
./binning | tee -a "$log_file"
./jackknife | tee -a "$log_file"

# Move program output and logs over to a specific directory where they can be
# stored.
[ ! -d  "$results_dir" ] && mkdir -p "$results_dir"
mv *.dat *.xyz *.out "$results_dir"

[ ! -d  "$log_dir" ] && mkdir -p "$log_dir"
mv *.log "$log_dir"
