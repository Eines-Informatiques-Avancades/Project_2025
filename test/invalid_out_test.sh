#!/bin/sh

#read -p 'Enter test name : ' test_name
test_name=test_template
out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

cd $test_name

./../../vdw_gas

echo ' '
mv *.dat positions.xyz ./run_out

cd run_out

echo ' ' 
echo "comparing files with test outputs"
echo ' ' 

for file in "${out_files[@]}"; do
   if [ -f $file ]; then
       echo $file
       nans=$(grep NaN $file | wc -l)   # check for NaN value
       echo NaN : '        ' $nans lines
       infinit=$(grep Infinity $file | wc -l)    # check for Infinity value 
       echo Infinity : '   ' $infinit lines
       heads=$(head -10 $file) 
       tails=$(tail -10 $file)
       if [ '$heads' == 'tails' ]; then    # check if all lines display same data
           echo print or propagation error 
       fi
   else
       echo File $file not found.
   fi
   echo ' '
done

