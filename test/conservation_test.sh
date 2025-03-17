#!/bin/sh

out_files=('temperature_inst.dat' 'total_energy.dat')

echo ' ' 
echo 'check of energy and temperature conservation'
echo ' '

echo '***************************'
for file in ${out_files[@]}; do 
    echo $file
    n_values=$(cat verified_samples/*rand_$file | awk '{print$2}' | wc -l)
    echo ' '
    echo number of samples = $n_values 
    echo ' ' 
    echo mean value : 
    cat verified_samples/*rand_$file | gnuplot -e 'stats "-" nooutput; print STATS_mean_y' 
    echo ' ' 
    echo standard deviation : 
    cat verified_samples/*rand_$file | gnuplot -e 'stats "-" nooutput; print STATS_stddev_y' 
    echo '****************************' 
done

