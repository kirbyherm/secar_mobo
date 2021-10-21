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

for TOSSES in $(seq 0 1 15) 
do
#  ./draw.py ../../output/output_4f_moead_FP2_FP3_150_"$TOSSES".csv
  cd 4f_FP2_FP3
  cosy ./pygmoCosy"$TOSSES".fox
  mv pic001.pdf X${TOSSES}.pdf
  mv pic002.pdf Y${TOSSES}.pdf
done

