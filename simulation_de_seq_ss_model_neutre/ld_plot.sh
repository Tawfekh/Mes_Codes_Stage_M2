#!/bin/bash
# ABDOU R WADE
# ligne de code pour dessiner le ld decay

cd ./LD_simulation_center_r484

echo "#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs" >> new_all_simul_10000_r76.stat

for i in * ; do sed '1d' ${i} ; done >> new_all_simul_10000_r76.stat

perl /home/tawfekh/PopLDdecay-3.31/bin/Plot_OnePop.pl -inFile new_all_simul_10000_r76.stat -keepR -output new_all_simul_10000_r76
