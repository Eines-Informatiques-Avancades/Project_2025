#!/bin/sh

out_files=("kinetic_energy.dat" "lj_potential.dat" "positions.xyz" "temperature_inst.dat" "total_energy.dat")

cd run_out

echo ' ' 
echo "comparing files with test outputs..."
echo ' ' 

for file in "${out_files[@]}"; do
   if [ -f $file ]; then
       echo $file
       nans=$(grep NaN $file | wc -l)   # check for NaN value
       echo NaN : '        ' $nans lines
       infinit=$(grep Infinity $file | wc -l)    # check for Infinity value 
       echo Infinity : '   ' $infinit lines
       heads=$(tail -n +2 $file | head -10 | awk '{print $2}')
       tails=$(tail -10 $file | awk '{print $2}')
       if [ "$heads" = "$tails" ]; then    # check if all lines display same data
           echo ERROR : print or propagation issue
       fi
   else
       echo File $file not found.
   fi
   echo ' '
done

