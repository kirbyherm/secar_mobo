#!/bin/bash
# Runs with $ sh quickrun.sh # # # ...
#
# Compile code.
#
# make
#
# Execute tests.
#
# echo -e "version\tn\tmax\tl2norm\ttime" 

start=100
stop=$(($start+10))

for TOSSES in $(seq $start 1 $stop) 
do
  ./draw.py ../output/output_4f_moead_FP2_FP3_225_"$TOSSES".csv
done

