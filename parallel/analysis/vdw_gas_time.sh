#!/bin/sh
#
# vdw_gas_time.sh
# Molecular Dynamics Simulation of a Van der Waals Gas
# Ricard Rodriguez
#
# Retrieves the elapsed simulation time from the cluster log files and formats
# them into a plottable data file.
#
# Requirements: awk, sed, grep, paste, rm, POSIX-compliant shell
#

tmp_labels_file="tmp_labels_${run_output}"
tmp_time_file="tmp_time_${run_output}"

get_times_vdw_gas() {
    timestamp_string="$2"
    run_output="$1"

    grep -A 1 "$timestamp_string" "$run_output" | \
        awk -F ':' '{print $2}' | \
        sed 's/ s$//g' | \
        awk 'NF { if (n++ % 2) print prev "\t" $0; else prev = $0 }'
}

format_time_output() {
    run_output="$1"
    time_output="time_${run_output}"

    # Output header.
    printf "# %s, %s, %s\n" \
        'Part' 'CPU time (s)' 'Wallclock time (s)' \
        > "$time_output"

    printf "%s\n" \
        'Initial configuration      ' \
        'System evolution           ' \
        'Post-trajectory analysis   ' \
        > "$tmp_labels_file"

    # Retrieve main simulation times (vdw_gas).
    get_times_vdw_gas "$run_output" 'Cputime: ' > "$tmp_time_file"

    paste -d '\t' \
        "$tmp_labels_file" "$tmp_time_file" \
        >> "$time_output"
}

cleanup(){
    rm -f "$tmp_labels_file" "$tmp_time_file"
}

# Remove previous time outputs.
rm -f time_*.out

for i in 1 2 4 8 16 32 40; do
    run_output="vdw_gas_${i}_core.out"
    time_output="time_${run_output}"

    [ -f "vdw_gas_${i}_core.out" ] && (
        # Create a file containing the duration of each part of a simulation
        # with a certain number of cores.
        format_time_output "$run_output"

        # Create files for each time measure across all simulations with
        # different amount of cores.
        while IFS= read -r line; do
            measure="$(echo "$line" | sed 's/[[:space:]]*$//' | tr '[:upper:]' '[:lower:]' | tr ' ' '_')"

            [ ! -f "time_${measure}.out" ] && \
                printf "# %s, %s, %s\n" \
                'Cores' 'CPU time (s)' 'Wallclock time (s)' \
                > "time_${measure}.out"

            grep "$line" "$time_output" | \
                awk -F '\t' -v cores="$i" '{print cores "\t" $2 "\t" $3}' \
                >> "time_${measure}.out"
        done < "$tmp_labels_file"
    )
done

cleanup
