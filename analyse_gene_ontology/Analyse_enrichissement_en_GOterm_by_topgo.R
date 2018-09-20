#########################################################################################################################
## @author: Abdou Rahmane WADE                                                                                         ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ##   Titre: Script pour analyser l'enrichissement en termes de GO   ##                    ##
##                           ##                                                                  ##                    ##
##                           ######################################################################                    ##
##                                                                                                                     ## 
#########################################################################################################################
rm(list=ls())

dir <-"/media/tawfekh/zap/H12_center/peaks_center"
setwd(dir)

gene_s_sel <- read.table("gene_sous_selection.txt",h=T)

GO_ref <- read.table('GO_ref.txt')

resrefall <- read.table("resrefall",h=T)


#### Code R pour faire un test d'enrichissement sur des régions génomique sous sélection

source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
biocLite("Rgraphviz")

### Chgarger les packages suivant
library(topGO)
library(Rgraphviz)

#########################
# Lecture des donnees   #
#########################
#setwd("..........")
#Lecture du fichier d'annotations
GFF<-GO_ref

#GFF<-GFF[,c(1,2)]
names(GFF)<-c("geneid","GOterm")

#Donne le nombre de genes annotes
length(unique(GFF$geneid)) 	


###table de donnes : liste de genes identifies par SweeD
data<-gene_s_sel
names(data)<-c("geneid")
dataDF<-as.data.frame(data)


#nb de genes dans data
length(unique(dataDF$geneid))

############################
# MEF donnees format TopGO #
############################
# liste des genes du GFF - elimination de la redondance
ListGFFunique<-as.vector(unique(GFF$geneid))
geneID2GO <- list(NULL)

for (i in 1 : length(ListGFFunique)) {
  temp<-GFF[GFF$geneid==ListGFFunique[i],]
  geneID2GO[[i]]<-as.character(temp$GOterm)
}
names(geneID2GO)<-as.character(ListGFFunique)

#Construction de la liste de genes d'interets
#TopGO peut aussi prendre en compte des scores associes aux contig pour faire les tests (par exemple les pvalue des tests d'association, cf tutoriel du package) 
geneNames2GO <- names(geneID2GO)
geneListdataGO <- factor(as.integer(geneNames2GO %in% dataDF$geneid))
names(geneListdataGO) <- geneNames2GO
length(which(geneListdataGO==1))  #Donne le nombre de contig sélectionnes et annotes donc utilisables pour l'analyse
str(geneListdataGO)


##################
# Analyses TopGO #
##################
# crÃ©ation d'un objet topGO data #### 
GOdata <- new("topGOdata",
              ontology = "CC",
              allGenes = geneListdataGO,
              nodeSize =5,   #dans notre cas, grd nb de genes 5 ou 10 ne change rien : pas de catégorie trop faiblement représentée - supprime les termes faiblement reprÃ©sentÃ©s (dans le tutotiel: valeur entre 5 et 10 donne des rÃ©su + stables) 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO
)

#Affichage des caractÃ©ristiques de l'objet TopGO
GOdata

#### Il existe 5 types de test et 6 algorythm sous topGO
##### tous les algorythmes ne marchent pas avec tous les test (Cf tutoriel)
#### j'ai choisi le test exact de Fisher qui est basÃ© sur contage de gÃ¨nes, 
#### et le test de Kolmogorov-Smirnov basÃ© sur les scores de gÃ¨nes
###J'ai choisi les algorythmes Classic et elim pour les test F et KS

##########################################################
#### Test de Fisher classique effectue terme a terme  ####
##########################################################


#### Test d'enrichissement classique de la surreprÃ©sentation des termes GO au sein du groupe de gÃ¨nes exprimÃ©s 
#### ceci de maniÃ¨re diffÃ©rentielle. donc chaque catÃ©gorie de GO est testÃ© de maniÃ¨re indÃ©pendante

Fisherclassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
Fisherclassic.table<-GenTable(GOdata, classicFisher = Fisherclassic, topNodes=20)
Fisherclassic.table
write.table(Fisherclassic.table, "resultFisherclassic-pmillet_center-CC.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(Fisherclassic), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, Fisherclassic, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE)


#### Tests de Fisher en tenant compte des liens de hierarchie des termes entre eux  ####

FisherWeight01<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
FisherWeight01.table<-GenTable(GOdata, Weight01 = FisherWeight01, topNodes=20)
FisherWeight01.table
write.table(FisherWeight01.table, "FisherWeight01-pmillet_center-CC.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(FisherWeight01), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, FisherWeight01, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE)

allRes<-GenTable(GOdata, classicFisher = Fisherclassic, weight01=FisherWeight01, classic=Fisherclassic, orderBy="weight01", ranksOf="classic", topNodes=10)
allRes   
write.table(allRes  , "pmillet_center_allRes.CC.txt" , sep=";", quote=FALSE)

##################
# Analyses TopGO # MF
##################
# crÃ©ation d'un objet topGO data #### 
GOdata <- new("topGOdata",
              ontology = "MF",
              allGenes = geneListdataGO,
              nodeSize =5,   #dans notre cas, grd nb de genes 5 ou 10 ne change rien : pas de catégorie trop faiblement représentée - supprime les termes faiblement reprÃ©sentÃ©s (dans le tutotiel: valeur entre 5 et 10 donne des rÃ©su + stables) 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO
)

#Affichage des caractÃ©ristiques de l'objet TopGO
GOdata

#### Il existe 5 types de test et 6 algorythm sous topGO
##### tous les algorythmes ne marchent pas avec tous les test (Cf tutoriel)
#### j'ai choisi le test exact de Fisher qui est basÃ© sur contage de gÃ¨nes, 
#### et le test de Kolmogorov-Smirnov basÃ© sur les scores de gÃ¨nes
###J'ai choisi les algorythmes Classic et elim pour les test F et KS

##########################################################
#### Test de Fisher classique effectue terme a terme  ####
##########################################################


#### Test d'enrichissement classique de la surreprÃ©sentation des termes GO au sein du groupe de gÃ¨nes exprimÃ©s 
#### ceci de maniÃ¨re diffÃ©rentielle. donc chaque catÃ©gorie de GO est testÃ© de maniÃ¨re indÃ©pendante

Fisherclassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
Fisherclassic.table<-GenTable(GOdata, classicFisher = Fisherclassic, topNodes=20)
Fisherclassic.table
write.table(Fisherclassic.table, "resultFisherclassic-pmillet_center-MF.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(Fisherclassic), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, Fisherclassic, firstSigNodes = 5, useInfo = "all_MF", pdfSW = TRUE)


#### Tests de Fisher en tenant compte des liens de hierarchie des termes entre eux  ####

FisherWeight01<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
FisherWeight01.table<-GenTable(GOdata, Weight01 = FisherWeight01, topNodes=20)
FisherWeight01.table
write.table(FisherWeight01.table, "FisherWeight01-pmillet_center-MF.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(FisherWeight01), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, FisherWeight01, firstSigNodes = 5, useInfo = "all-MF", pdfSW = TRUE)

allRes<-GenTable(GOdata, classicFisher = Fisherclassic, weight01=FisherWeight01, classic=Fisherclassic, orderBy="weight01", ranksOf="classic", topNodes=10)
allRes   
write.table(allRes  , "pmillet_center_allRes.MF.txt" , sep=";", quote=FALSE)


##################
# Analyses TopGO # BP
##################
# crÃ©ation d'un objet topGO data #### 
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneListdataGO,
              nodeSize =5,   #dans notre cas, grd nb de genes 5 ou 10 ne change rien : pas de catégorie trop faiblement représentée - supprime les termes faiblement reprÃ©sentÃ©s (dans le tutotiel: valeur entre 5 et 10 donne des rÃ©su + stables) 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO
)

#Affichage des caractÃ©ristiques de l'objet TopGO
GOdata

#### Il existe 5 types de test et 6 algorythm sous topGO
##### tous les algorythmes ne marchent pas avec tous les test (Cf tutoriel)
#### j'ai choisi le test exact de Fisher qui est basÃ© sur contage de gÃ¨nes, 
#### et le test de Kolmogorov-Smirnov basÃ© sur les scores de gÃ¨nes
###J'ai choisi les algorythmes Classic et elim pour les test F et KS

##########################################################
#### Test de Fisher classique effectue terme a terme  ####
##########################################################


#### Test d'enrichissement classique de la surreprÃ©sentation des termes GO au sein du groupe de gÃ¨nes exprimÃ©s 
#### ceci de maniÃ¨re diffÃ©rentielle. donc chaque catÃ©gorie de GO est testÃ© de maniÃ¨re indÃ©pendante

Fisherclassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
Fisherclassic.table<-GenTable(GOdata, classicFisher = Fisherclassic, topNodes=20)
Fisherclassic.table
write.table(Fisherclassic.table, "resultFisherclassic-pmillet_center-BP.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(Fisherclassic), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, Fisherclassic, firstSigNodes = 5, useInfo = "all_BP", pdfSW = TRUE)


#### Tests de Fisher en tenant compte des liens de hierarchie des termes entre eux  ####

FisherWeight01<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
FisherWeight01.table<-GenTable(GOdata, Weight01 = FisherWeight01, topNodes=20)
FisherWeight01.table
write.table(FisherWeight01.table, "FisherWeight01-pmillet_center-BP.table.txt" , sep=";", quote=FALSE)
showSigOfNodes(GOdata, score(FisherWeight01), firstSigNodes = 4, useInfo ='all')
printGraph(GOdata, FisherWeight01, firstSigNodes = 5, useInfo = "all-BP", pdfSW = TRUE)

allRes<-GenTable(GOdata, classicFisher = Fisherclassic, weight01=FisherWeight01, classic=Fisherclassic, orderBy="weight01", ranksOf="classic", topNodes=10)
allRes   
write.table(allRes  , "pmillet_center_allRes.BP.txt" , sep=";", quote=FALSE)
