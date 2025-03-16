#!/bin/sh

read -p 'Enter test name: ' test_name

cd $test_name

./../../vdw_gas

echo " "

mv *.dat positions.xyz ./run_out

bash ./../invalid_out_test.sh g
bash ./../reproductibility_test.sh g
