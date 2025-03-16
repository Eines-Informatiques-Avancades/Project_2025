#!/bin/sh

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat" "rdf.dat" "rmsd.dat")

echo ' '
echo 'testing reproductibility...'
echo ' '

for file in ${out_files[@]}; do
    if ! [ -f verified_samples/*rand_$file ]; then
        echo ERROR : $file testing samples not found
    else
        while read line; do
            match=$(grep "$line" run_out/$file)
            if [ -f "$match" ]; then
                fail=TRUE
            fi
        done < verified_samples/*rand_$file
        if [ "$fail" == 'TRUE' ]; then
            echo FAILED in $( echo $file | cut -d '.' -f 1)
        else
            echo VERIFIED in $( echo $file | cut -d '.' -f 1)
        fi
    fi
    echo ' '
done


