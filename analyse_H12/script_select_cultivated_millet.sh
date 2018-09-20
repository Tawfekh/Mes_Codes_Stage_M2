#!/bin/bash
# ABDOU R WADE
# ligne de code pour extraire les mils du centre du Sahel

chr=$1

mkdir cult_centerVCF

vcftools --vcf chr${chr}_sample_190_of_allSNP_chr${chr}_1.vcf --keep cult_c.txt --recode --out ./cult_centerVCF/center_chr${chr}
