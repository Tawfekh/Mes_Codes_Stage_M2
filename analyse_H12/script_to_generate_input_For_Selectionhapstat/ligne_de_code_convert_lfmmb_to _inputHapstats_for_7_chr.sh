#!/bin/bash
# ABDOU R WADE
# ligne de code pour executer convertir lfmmb to inputHapstats

for i in {1..7}; do Rscript 18-04-18_Script_convert_lfmmb_to_inputHapstats.R chr${i}_sample_190_of_allSNP_chr${i}.txt; done
