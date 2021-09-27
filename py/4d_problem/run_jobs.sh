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

for TOSSES in $(seq 0 1 13) 
do
#  ./draw.py ../../output/output_4f_moead_FP2_FP3_150_"$TOSSES".csv
  ./cosy_draw.py ./4f_FP2_FP3/pygmoCosy"$TOSSES".fox
done

