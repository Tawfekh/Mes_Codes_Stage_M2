#!/bin/bash
# ABDOU R WADE
# Code pour calculer le ld decay en fonction de la distance entre deux loci

for i in *.vcf ; do PopLDdecay -InVCF ${i} -MaxDist 100 -MAF 0.05 -Miss 0.30 -Methold 2 -OutStat $(basename -s .vcf $i) ; done

mkdir ./LD_simulation_center_r484

mv *.gz ./LD_simulation_center_r484/

cd ./LD_simulation_center_r484

gzip -d *.gz
