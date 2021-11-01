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

GENS=1500
BATCH=220

mkdir results_"$BATCH"
./make_db.py $GENS $BATCH
./view_db.py $BATCH
mv best"$BATCH".h5 results_"$BATCH"/
./draw.py results_"$BATCH"/best"$BATCH".h5
cd results_"$BATCH"
cd profiles

for TOSSES in $(seq 0 1 10) 
do
#  ./draw.py ../../output/output_4f_moead_FP2_FP3_150_"$TOSSES".csv
  cosy ./pygmoCosy"$TOSSES".fox > "$TOSSES".txt
  mv pic001.pdf X${TOSSES}.pdf
  mv pic002.pdf Y${TOSSES}.pdf
done
rm -r ./*.lis
