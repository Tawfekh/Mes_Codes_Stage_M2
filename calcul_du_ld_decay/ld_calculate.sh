#!/bin/bash
# ABDOU R WADE
# Code pour calculer le ld decay en fonction de la distance entre deux loci

chr=$1

mkdir LD_chr_cult_center${chr}

for line in $(cut -d' '  -f3 pos_chr${chr}.txt) ; do vcftools --vcf chr${chr}_sample_190_of_allSNP_chr${chr}_1.vcf --chr chr${chr} --keep cult_c.txt --from-bp ${line} --to-bp $((${line} + 500000)) --recode --out ./LD_chr_cult_center${chr}/center_chr${chr}_${line} ; done

cd LD_chr_cult_center${chr}

for i in *.vcf ; do PopLDdecay -InVCF ${i} -MaxDist 100 -MAF 0.05 -Miss 0.30 -Methold 2 -OutStat $(basename -s .recode.vcf $i) ; done

mkdir ../list_all_chr_cult_center

cp *.gz ../list_all_chr_cult_center/

cd ../list_all_chr_cult_center

gzip -d *.gz

