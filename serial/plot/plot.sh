#!/bin/sh
#
# plot.sh
# Molecular Dynamics Simulation of a Van der Waals Gas
# Ricard Rodriguez
#

# Paths relative to the 'plot/' subdirectory.
plot_dir="$(dirname $(realpath "$0"))" # Full path to plot/
results_dir="../output"


# Print error message to stderr and exit the script.
err_exit() {
    printf >&2 "Error: %s\\n" "$@"
    exit 1
}

# Check if a given program is installed (if an executable with the provided
# name is in $PATH) and exit with an error message if it doesn't.
exec_exists() {
    command -v "$1" >/dev/null || \
        err_exit "$1 is not installed or is not in \$PATH."
}

exec_exists 'python'

printf "%s\\n" 'Executing plot scripts...'

cd "$plot_dir"

python "system_energy_evolution.py"
python "temperature_evolution.py"
python "rdf.py" "$results_dir/rdf.dat"
python "rmsd.py" "$results_dir/rmsd.dat"
python "binning.py" "$results_dir/binning_thermodynamics.dat"
# python "visualize_coordinates.py" "$results_dir/positions.xyz"
# python "particle_movement_animation.py" "$results_dir/positions.xyz"
