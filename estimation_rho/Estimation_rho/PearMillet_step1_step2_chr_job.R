
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

args <- commandArgs(TRUE)

t0 <- proc.time()

nb_run <- args[1]

run <- paste("run", nb_run, sep="")

print(run)

# Fichier contenant que les positions 

position_chr1 <- read.table("./position_chr1.txt", h=T)
position_chr2 <- read.table("./position_chr2.txt", h=T)
position_chr3 <- read.table("./position_chr3.txt", h=T)
position_chr4 <- read.table("./position_chr4.txt", h=T)
position_chr5 <- read.table("./position_chr5.txt", h=T)
position_chr6 <- read.table("./position_chr6.txt", h=T)
position_chr7 <- read.table("./position_chr7.txt", h=T)



    # Echantillonnage aléatoire de 15 positions de début de fenêtres
    sample_chr1 <- sample(1:nrow(position_chr1),size = 1, replace = T)
    position_sample_chr1 <- position_chr1[sample_chr1,]
    
    sample_chr2 <- sample(1:nrow(position_chr2),size = 1, replace = T)
    position_sample_chr2 <- position_chr2[sample_chr2,]
    
    sample_chr3 <- sample(1:nrow(position_chr3),size = 1, replace = T)
    position_sample_chr3 <- position_chr3[sample_chr3,]
    
    sample_chr4 <- sample(1:nrow(position_chr4),size = 1, replace = T)
    position_sample_chr4 <- position_chr4[sample_chr4,]
    
    sample_chr5 <- sample(1:nrow(position_chr5),size = 1, replace = T)
    position_sample_chr5 <- position_chr5[sample_chr5,]
    
    sample_chr6 <- sample(1:nrow(position_chr6),size = 1, replace = T)
    position_sample_chr6 <- position_chr6[sample_chr6,]
    
    sample_chr7 <- sample(1:nrow(position_chr7),size = 1, replace = T)
    position_sample_chr7 <- position_chr7[sample_chr7,]
    
      dir.create("./FastEPRR_extra")
      write.table(paste(paste("position_sample_chr1", run, sample_chr1, sep=" "),
                        paste("position_sample_chr2", run, sample_chr2, sep=" "),
                        paste("position_sample_chr3", run, sample_chr3, sep=" "),
                        paste("position_sample_chr4", run, sample_chr4, sep=" "),
                        paste("position_sample_chr5", run, sample_chr5, sep=" "),
                        paste("position_sample_chr6", run, sample_chr6, sep=" "),
                        paste("position_sample_chr7", run, sample_chr7, sep=" "), sep = "\n"),
                  file = paste("./FastEPRR_extra/Position_sample_", run, ".txt",sep=""), 
                  quote=FALSE, col.names=F, row.names=F, append=TRUE)

    # conversion des positions échantillonnées de pb en kb
    
      position_sample <- c(position_sample_chr1,
                           position_sample_chr2,
                           position_sample_chr3,
                           position_sample_chr4,
                           position_sample_chr5,
                           position_sample_chr6,
                           position_sample_chr7)
      
    position_sample <- 0.001 * position_sample
    
    # Choix fenêtre format FastEPRR
    erStart <- as.character(position_sample)
    erEnd <- as.character(99 + position_sample)
    
# Fichier contenant que les individus du centre
cult_centre <- read.table("./cult_c.txt")
cult_centre <- cult_centre$V1

    # Choix Individus format FastEPRR
    idvlConsidered <- c()
    idvlConsidered <- paste(idvlConsidered, cult_centre[1], sep="")
    for(i in 2:88){
      idvlConsidered <- paste(idvlConsidered, cult_centre[i], sep=";")
    }

# Code FastEPRR pour estmer rhô pour le premier
    
    library(FastEPRR)
    for(i in 1:7)
    {
      input <- paste("./chr", i, "_sample_190_of_allSNP_chr", i, ".vcf.gz", sep="")
      print(input)
      dir_output <- paste("./step1_", run,"_chr",i, sep="")
      dir.create(dir_output)
      output <- paste(dir_output,"/srcOutputFilePath_", run, "_chr", i, sep="")
      dir.create(output)
      print(output)
      
      FastEPRR_VCF_step1(vcfFilePath= input, erStart = erStart[i], erEnd = erEnd[i],
                       idvlChrFormat="[1:0]", winLength="100", idvlConsidered = idvlConsidered,
                       srcOutputFilePath= output)
      
      input_2 <- paste("./step1_", run,"_chr",i, sep="")
      print(input_2)
      output_2 <- paste("./step2_DXoutput_", run, "_chr", i, sep="")
      dir.create(output_2)
      print(output_2)

      FastEPRR_VCF_step2(srcFolderPath= input_2, demoParameter ="-G 2.517 -eG 0.939 0.000 -eN 0.953 3.12", 
                         DXOutputFolderPath=output_2)
      print(paste("done in", output_2, sep =" "))
    }

    
t1 <- proc.time(); t <- t1 - t0

write.table(t[3], file = paste("./FastEPRR_extra/time", run, ".txt"))
