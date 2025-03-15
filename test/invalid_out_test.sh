#!/bin/sh

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

cd run_out || exit

echo ' '
echo 'Comparing output files with test outputs...'
echo ' '

for file in "${out_files[@]}"; do
    if [ -f "$file" ]; then
        echo "$file"

        # Check for NaN values.
        nans=$(grep NaN "$file" | wc -l)
        echo "NaN :          $nans lines"

        # Check for Infinity values.
        infinit=$(grep Infinity "$file" | wc -l)
        echo "Infinity :     $infinit lines"

        # Check if all the lines in $file display same data.
        heads=$(tail -n +2 "$file" | head -10 | awk '{print $2}')
        tails=$(tail -10 "$file" | awk '{print $2}')
        if [ "$heads" = "$tails" ]; then
            echo ERROR : print or propagation issue
        fi
    else
        echo "File $file not found."
    fi

    echo ' '
done
