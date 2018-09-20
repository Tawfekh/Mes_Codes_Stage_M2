#!/bin/bash
# ABDOU R WADE
# ligne de code pour dessiner le scan genomique H12 et pour visualiser les haplotypes

output=$1
seuil=$2
outpeak=$3
scan=$4
hapspectr=$5

python H12peakFinder.py $output -o $outpeak -t $seuil

Rscript H12_viz.R $output $outpeak $scan 50

Rscript hapSpectrum_viz.R $outpeak $hapspectr 10 190
