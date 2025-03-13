#!/bin/sh

read -p 'Enter test name: ' test_name

exec 2>/dev/null 

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

cd $test_name

./../../vdw_gas

echo " "

mv *.dat .
mv positions.xyz . 

echo " "
echo "comparing files with test outputs"
echo " "

for file in "${out_files[@]}"; do
    if [ -f $file ]; then
        echo $file
        # check for NaN value
        nans=$(grep NaN $file | wc -c)
        echo NaN : $nans
        # check for Infinity value 
        infinit=$(grep Infinity $file | wc -c)
        echo Infinity : $infinit
        # check reproductibility issues
        differences=$(diff -r ./$file out/$file | grep "<" | wc -l)
        echo Differences : $((differences/2))
    else
        echo File $file not found.
    fi

    echo " "
done
exit 0
