#!/bin/sh
#
# plot.sh
#
#

# Paths relative to the 'plot/' subdirectory.
results_dir="../output"
plot_dir="."


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

# python "$plot_dir/visualize_coordinates.py" "../positions.xyz"
python "$plot_dir/system_energy_evolution.py"
python "$plot_dir/temperature_evolution.py"
# python "$plot_dir/particle_movement_animation.py" "../positions.xyz"
