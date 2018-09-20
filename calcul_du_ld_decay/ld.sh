#!/bin/bash
# ABDOU R WADE
# ligne de code pour executer les codes ld_calculate.sh & ld_plot.sh pour les 7 chr

for i in {1..7} ; do bash ld_calculate.sh ${i} ; done

bash ld_plot.sh
