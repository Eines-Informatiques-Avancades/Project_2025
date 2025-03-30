#!/bin/bash

out_files=("positions.xyz" "thermodynamics.dat" "rdf.dat" "rmsd.dat" "binning_thermodynamics.dat")

echo ' '
echo 'testing reproductibility...'
echo ' '

for file in ${out_files[@]}; do
    fail=FALSE
    if [ -f verified_samples/*rand_$file ]; then
        if [ -f run_out/$file ]; then
            while read line; do
                match=$(grep "$line" run_out/$file)
                if [ -z "$match" ]; then
                    fail=TRUE
                fi
            done < verified_samples/*rand_$file
            if [ "$fail" == 'TRUE' ]; then
                echo FAILED in $( echo $file | cut -d '.' -f 1)
            else
                echo VERIFIED in $( echo $file | cut -d '.' -f 1)
            fi
        else
            echo ERROR: $file not found
        fi
    else
        echo ERROR verified sample for $file not found
    fi
done


