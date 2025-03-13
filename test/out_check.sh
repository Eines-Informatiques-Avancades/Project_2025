#!/bin/bash

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

for i in "${out_files[@]}"; do
    if [ -f $i ]; then
        echo File $i exists.
    else
        echo File $i does not exist.
    fi
done

