#!/bin/bash

read -p 'Enter test name: ' test_name

cd $test_name

./../../vdw_gas
./../../binning

echo " "

mv *.dat positions.xyz run_out

./../invalid_out_test.sh 
./../reproductibility_test.sh 
