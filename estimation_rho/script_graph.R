#############################################################################################################
## Author:Abdou Rahmane WADE                                                             Date:      2018   ##
##                                                                                                         ##
##          ##################################################################################             ##
##          ## Titre: Script pour générer des graph de distribution de Rho estimé           ##             ## 
##          ##                           Description:                                       ##             ##
##          ##################################################################################             ##
##                                                                                                         ##
#############################################################################################################
rm(list = ls())
setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data/vcfFile/step3_allrun")

data <- read.csv2("rho_allchr.csv", h=T)

par(mfrow = c(1,1))

plot(density(data$RHO))

pdf("sorties_Distribution_rho_estimés_genome.pdf")
plot(density(data$RHO), cex=5,
     xlab = "Rho", main = "Distribution des rho estimés du génome")
legend(50, 0.016, legend= paste("All Chromosomes", "mean =", round(mean(data$RHO),2), "sd =", round(sd(data$RHO),2), sep=" "), bty = "n", cex = 0.75)
dev.off()

pdf("sorties_estimationsDistribution_rho.pdf")
i=1
plot(density(data$RHO[(((i-1)*15)+1):(i*15)]), col=i+1, cex=5, ylim =c(0,0.027),
     xlab = "Rho", main = "Distribution des rho estimés pour chaque chromosome")
legend(50, (0.027-(i*0.0012)), legend= paste(data$CHR[i*15], "mean =", round(mean(data$RHO[(((i-1)*15)+1):(i*15)]),2), 
                                              "sd =", round(sd(data$RHO[(((i-1)*15)+1):(i*15)]),2), sep=" "), 
       bty = "n", cex = 0.75, text.col=i+1)
for(i in 2:7){
  lines(density(data$RHO[(((i-1)*15)+1):(i*15)]), col=i+1, cex=5)
  legend(50, (0.027-(i*0.0012)), legend= paste(data$CHR[i*15], "mean =", round(mean(data$RHO[(((i-1)*15)+1):(i*15)]),2), 
                                                "sd =", round(sd(data$RHO[(((i-1)*15)+1):(i*15)]),2), sep=" "), 
         bty = "n", cex = 0.75, text.col=i+1)
}
dev.off()


dchr1 <- data[1:15,]
dchr2 <- data[16:30,]
dchr3 <- data[31:45,]
dchr4 <- data[46:60,]
dchr5 <- data[61:75,]
dchr6 <- data[76:90,]
dchr7 <- data[91:105,]

pdf("sorties_Distribution_rho_sur_chr.pdf")
par(mfrow = c(4,2))
attach(dchr1)
plot(data=dchr1[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,275407), ylim=c(0,576),
     main = "Distribution rho sur chr1", col="red", pch=20, xlab = "Position (kb)")
detach(dchr1)
attach(dchr2)
plot(data=dchr2[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,249546), ylim=c(0,576),
     main = "Distribution rho sur chr2", col="red", pch=20, xlab = "Position (kb)")
detach(dchr2)
attach(dchr3)
plot(data=dchr3[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,297401), ylim=c(0,576),
     main = "Distribution rho sur chr3", col="red", pch=20, xlab = "Position (kb)")
detach(dchr3)
attach(dchr4)
plot(data=dchr4[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,203198), ylim=c(0,576),
     main = "Distribution rho sur chr4", col="red", pch=20, xlab = "Position (kb)")
detach(dchr4)
attach(dchr5)
plot(data=dchr5[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,164707), ylim=c(0,576),
     main = "Distribution rho sur chr5", col="red", pch=20, xlab = "Position (kb)")
detach(dchr5)
attach(dchr6)
plot(data=dchr6[order(PosDebut),], RHO~PosDebut, type='b', xlim =c(1,270064), ylim=c(0,576),
     main = "Distribution rho sur chr6", col="red", pch=20, xlab = "Position (kb)")
detach(dchr5)
attach(dchr7)
plot(data=dchr7[order(PosDebut),], RHO~PosDebut, type='b', xlim=c(1,226187), ylim=c(0,576),
     main = "Distribution rho sur chr7", col="red", pch=20, xlab = "Position (kb)")
detach(dchr7)
dev.off()

Pos <- data.frame(table(data$PosDebut))

Pos$Var1 <- as.numeric(Pos$Var1)

pdf("sorties_estimations_rho_1.pdf")
par(mfrow = c(4,2))
plot(1:275406.818, rep(2,275406 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[1:15], Pos$Freq[1:15], size =1.5, pch=20)
lines(data$PosDebut[1:15], Pos$Freq[1:15], size =2, pch=20, col='blue')
points(data$PosDebut[6], Pos$Freq[6], size =1.5, pch=20, col="red")
title(main = "Estimation rhô \n Positions échantillonage sur chr1")

plot(1:249546, rep(2,249546 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[16:30], Pos$Freq[16:30], size =1.5, pch=20)
lines(data$PosDebut[16:30], Pos$Freq[16:30], size =2, pch=20, col='blue')
title(main = "Estimation rhô \n Positions échantillonage sur chr2")

plot(1:297401, rep(2,297401 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[31:45], Pos$Freq[31:45], size =1.5, pch=20)
lines(data$PosDebut[31:45], Pos$Freq[31:45], size =2, pch=20, col='blue')
title(main = "Estimation rhô \n Positions échantillonage sur chr3")

plot(1:203198, rep(2,203198 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[46:60], Pos$Freq[46:60], size =1.5, pch=20)
lines(data$PosDebut[46:60], Pos$Freq[46:60], size =2, pch=20, col='blue')
title(main = "Estimation rhô \n Positions échantillonage sur chr4")

plot(1:164707, rep(2,164707 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[61:75], Pos$Freq[61:75], size =1.5, pch=20)
lines(data$PosDebut[61:75], Pos$Freq[61:75], size =2, pch=20, col='blue')
points(data$PosDebut[65], Pos$Freq[65], size =1.5, pch=20, col="red")
title(main = "Estimation rhô \n Positions échantillonage sur chr5")

plot(1:270064, rep(2,270064 ),ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[76:90], Pos$Freq[76:90], size =1.5, pch=20)
lines(data$PosDebut[76:90], Pos$Freq[76:90], size =2, pch=20, col='blue')
points(data$PosDebut[88], Pos$Freq[88], size =1.5, pch=20, col="red")
title(main = "Estimation rhô \n Positions échantillonage sur chr6")

plot(1:226187, rep(2,226187 ), ylim = c(0.5,1.5), yaxt="n" , xlab = "Position (kb)", ylab = "" )
points(data$PosDebut[91:105], Pos$Freq[91:105], size =1.5, pch=20)
lines(data$PosDebut[91:105], Pos$Freq[91:105], size =2, pch=20, col='blue')
title(main = "Estimation rhô \n Positions échantillonage sur chr7")

dev.off()
