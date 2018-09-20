Example of script using scratch space

 

#!/bin/sh



############      SGE CONFIGURATION      ###################

# Ecrit les erreur dans le fichier de sortie standard 

#$ -j y 


# Shell que l'on veut utiliser 

#$ -S /bin/bash 


# Email pour suivre l'execution 

#$ -M abdou-rahmane.wade@ird.fr


# Type de massage que l'on reçoit par mail

#    -  (b) un message au demarrage

#    -  (e) a la fin

#    -  (a)  en cas d'abandon

#$ -m bea 


# Queue que l'on veut utiliser

#$ -q bioinfo.q


# Nom du job

#$ -N zap_h12_chr6

############################################################

 

###### Creation du repertoire temporaire sur noeud

mkdir /scratch/h12/

scp -rp nas:/home/wade/data/mil_data/donnee_h12/chr6_sample_190_of_allSNP_chr6_h12.txt /scratch/h12/
scp -rp master0:/usr/local/SelectionHapStats-1.0/scripts/H12_H2H1.py /scratch/h12/


###### Execution du programme

python2.7 /scratch/h12/H12_H2H1.py /scratch/h12/chr6_sample_190_of_allSNP_chr6_h12.txt -o /scratch/h12/chr6_sample_190_of_allSNP_chr6_H12_output_400_50.txt 190 -w 400 -j 50 -d 0




##### Transfert des données du noeud vers master

scp -rp /scratch/h12/chr6_sample_190_of_allSNP_chr6_H12_output_400_50.txt nas:/home/wade/data/mil_data/donnee_h12/


#### Suppression du repertoire tmp noeud

rm -rf /scratch/h12/

