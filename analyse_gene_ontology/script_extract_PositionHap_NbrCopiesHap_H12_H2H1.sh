#!/bin/bash
# ABDOU R WADE
# ligne de code pour extraire Les Positions des Haplotypes le nombre de copies et leur valeur de H12 & H2H1

chr=$1

cut -d$'\t' -f2,3,4,9,10 center_chr${chr}_h12_output_400_50.txt >> chr${chr}_peaks
