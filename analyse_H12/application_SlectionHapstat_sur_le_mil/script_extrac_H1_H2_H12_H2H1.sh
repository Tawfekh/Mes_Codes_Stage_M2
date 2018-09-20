#!/bin/bash
# ABDOU R WADE
# ligne de code pour extraire Les stats H1 H2 H12 et H2/H1

doc=$1

chr=$2

cut -d$'\t' -f7 $doc > H1_$(basename $chr)

cut -d$'\t' -f8 $doc > H2_$(basename $chr)

cut -d$'\t' -f9 $doc > H12_$(basename $chr)

cut -d$'\t' -f10 $doc > H2:H1_$(basename $chr)
