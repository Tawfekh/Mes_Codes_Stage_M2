
#############################################################################################################
## Author:Abdou Rahmane WADE                                                             Date:      2018   ##
##                                                                                                         ##
##          ##################################################################################             ##
##          ## Titre: Script pour estimer le rho sur 15 fenêtres de 100kb prises au hasard  ##             ## 
##          ##                      le long de chaque chromosomes                           ##             ##
##          ##                           Description:                                       ##             ##
##          ##################################################################################             ##
##                                                                                                         ##
##                                                                                                         ##
##     Package FastEPRR: a Fast Estimator for the Population Recombination Rate Version 1.0                ##
##              args <- numero run                                                                         ##
##              InputFiles:                                                                                ##
##                          positions_chr${i}.txt :Fichier contenant que les positions des SNPs            ##
##                          cult_c.txt : Fichier contenant que les individus du centre                     ##
##                          chr${i}_sample_190_of_allSNP_chr${i}.vcf.gz: Fichiers vcf                      ##
##                                                                                                         ##
##                                                                                                         ##
##                                                                                                         ##
##              OutputFiles:                                                                               ##
##                          Position_sample_run${j}.txt : renseigne de la position de début                ##
##                                                        des fenètres échantionnées                       ##
##                          srcOutputFilePath_run${j}_chr{i} : sortie FastEPRR step1                       ##
##                          step2_DXoutput_run${j}_chr{i} : sortie FastEPRR step2                          ##
##                          time_run${i}.txt : temps d'execution                                           ##
##                                                                                                         ##
##                                                                                                         ##
##                                                                                                         ##
#############################################################################################################

# Code FastEPRR pour estmer rhô pour le premier
    
    library(FastEPRR)

      output <- "./step3_allrun"
      dir.create(output)
      print(output)
      
      FastEPRR_VCF_step3 (srcFolderPath="./step1_allrun", DXFolderPath="./step2_allrun",finalOutputFolderPath="./step3_allrun")


