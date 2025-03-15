#!/bin/sh

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

echo sampling data for tests

read -p "Enter number of desired sample lines: " n_samples

for file in ${out_files[@]}; do
    shuf -n $n_samples $file > $n_samples"rand_"$file
done

echo sampling accomplished

