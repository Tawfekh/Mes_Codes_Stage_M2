#!/bin/bash
# @uthor: Abdou R WADE
# Script pour executer le code Script_convert_FastEPRR_InputFile_vcf_from_genotype_format.R pour les 7 chromosomes


for i in {1..7}; do Rscript Script_convert_vcf_from_genotype_format.R chr${i}_sample_190_of_allSNP_chr${i}.txt; done
