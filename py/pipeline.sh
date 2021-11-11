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

GENS=1000
BATCH=290

mkdir results_"$BATCH"
./make_db.py $GENS $BATCH
./view_db.py $BATCH
mv magnet_factors.csv results_"$BATCH"/
mv magnet_values.csv results_"$BATCH"/
mv best"$BATCH".h5 results_"$BATCH"/
./draw.py results_"$BATCH"/best"$BATCH".h5
cd results_"$BATCH"
cd profiles

for TOSSES in $(seq 0 1 5) 
do
#  ./draw.py ../../output/output_4f_moead_FP2_FP3_150_"$TOSSES".csv
  cosy ./pygmoCosy"$TOSSES".fox > "$TOSSES".txt
  mv pic001.pdf X${TOSSES}.pdf
  mv pic002.pdf Y${TOSSES}.pdf
  cosy ./pygmoCosy"$TOSSES"_DE.fox > "$TOSSES"_DE.txt
  mv pic001.pdf X${TOSSES}_DE.pdf
#  mv pic002.pdf Y${TOSSES}_DE.pdf
  rm -r ./pic*.pdf
done
rm -r ./*.lis
