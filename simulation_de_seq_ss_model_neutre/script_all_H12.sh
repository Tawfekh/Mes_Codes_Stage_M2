#!/bin/bash
# ABDOU R WADE
# ligne de code pour regrouper les H12 de toute les 10000 simul

name=$1

cd 190cult

for i in *out.txt; do cat $i >> ${name}_allH12.txt; done

cp ${name}_allH12.txt ../
