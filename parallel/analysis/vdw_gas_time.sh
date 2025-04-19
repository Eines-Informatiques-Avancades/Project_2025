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

get_times() {
    run_output="$1"

    grep -A 1 'Cputime: ' "$run_output" | \
        awk -F ':' '{print $2}' | \
        sed 's/ s$//g' | \
        awk 'NF { if (n++ % 2) print prev "\t" $0; else prev = $0 }'
}

format_time_output() {
    run_output="$1"
    time_output="time_${run_output}"
    tmp_labels_file="tmp_labels_${run_output}"
    tmp_time_file="tmp_time_${run_output}"

    printf "# %s, %s, %s\n" \
        'Part' 'CPU time (s)' 'Wallclock time (s)' \
        > "$time_output"

    printf "%s\n" \
        'Initial configuration      ' \
        'System evolution           ' \
        'Post-trajectory analysis   ' \
        > "$tmp_labels_file"

    get_times "$run_output" > "$tmp_time_file"

    paste -d '\t' \
        "$tmp_labels_file" "$tmp_time_file" \
        >> "$time_output"

    rm "$tmp_labels_file" "$tmp_time_file"
}


for i in 1 2 4 8 16 32 40; do
    [ -f "vdw_gas_${i}_core.out" ] && format_time_output "vdw_gas_${i}_core.out"
done
