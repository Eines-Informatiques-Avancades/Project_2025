#!/bin/bash

out_files=("positions.xyz" "thermodynamics.dat" "rdf.dat" "rmsd.dat" "binning_thermodynamics.dat")

echo sampling data for tests

read -p "Enter number of desired sample lines: " n_samples

for file in ${out_files[@]}; do
    if ! [ -f $file ]; then
        echo $file for sampling not found
    else
        shuf -n $n_samples $file > $n_samples"rand_"$file
    fi
done

echo sampling accomplished

