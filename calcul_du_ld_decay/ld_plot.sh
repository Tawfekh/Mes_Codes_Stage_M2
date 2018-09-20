#!/bin/bash
# ABDOU R WADE
# ligne de code pour dessiner le ld decay

cd ./list_all_chr_cult_center

echo "#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs" >> all_frag_chr_cult_center.stat

for i in * ; do sed '1d' ${i} ; done >> all_frag_chr_cult_center.stat

perl /home/tawfekh/PopLDdecay-3.31/bin/Plot_OnePop.pl -inFile all_frag_chr_cult_center.stat -keepR -output all_frag_chr_cult_center_LDdecay
