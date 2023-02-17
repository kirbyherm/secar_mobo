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
#GENS=20
BATCH=510
#BATCH=8

mkdir results_"$BATCH"
./make_db.py $GENS $BATCH
./analyze_db.py $BATCH
mv magnet_factors.csv results_"$BATCH"/
mv best"$BATCH".h5 results_"$BATCH"/
./draw.py results_"$BATCH"/best"$BATCH".h5
#./draw_cluster.py results_"$BATCH"/best"$BATCH".h5 results_280/best280.h5
./draw_full.py results_"$BATCH"/best"$BATCH".h5 $BATCH
./draw_cluster_inverse.py results_"$BATCH"/best"$BATCH".h5 results_280/best280.h5
./plot_tsne.py results_"$BATCH"/best"$BATCH".h5 results_280/best280.h5
#./plot_corr.py results_"$BATCH"/best"$BATCH".h5 results_280/best280.h5
cd results_"$BATCH"
cd profiles
#
for TOSSES in $(seq 0 1 4) 
do
  cosy ./pygmoCosy"$TOSSES".fox > "$TOSSES".txt
  mv pic001.pdf X${TOSSES}.pdf
  mv pic002.pdf Y${TOSSES}.pdf
  cosy ./pygmoCosy"$TOSSES"_DE.fox > "$TOSSES"_DE.txt
  mv pic001.pdf X${TOSSES}_DE.pdf
#  cosy ./pygmoCosy"$TOSSES"_DE_FP1.fox > "$TOSSES"_DE_FP1.txt
#  mv pic001.pdf X${TOSSES}_DE_FP1.pdf
#  mv pic002.pdf Y${TOSSES}_DE.pdf
  rm -r ./pic*.pdf
done
rm -r ./*.lis
cd ../../
tar -czvf results_"$BATCH".tar.gz results_"$BATCH"/

