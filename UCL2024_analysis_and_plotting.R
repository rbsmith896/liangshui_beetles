#RB Smith
#MSc Thesis 2024: Liangshui Beetles
#Analysis and Plotting

library(dplyr)
library(rarestR)
library(ggplot2)
library(abdiv)
library(tidyverse)
library(viridis)
library(fossil)
library(vegan)
library(scales)
library(gridExtra)

#Figure 1: abundance of each species across all plots
tbeets <- as.data.frame(t(as.matrix(beets)))
tbeets$abdunance <- rowSums(tbeets)
tbeets <- tbeets[order(tbeets$abdunance), ]
rownames(tbeets) <- paste("ms",rownames(tbeets),sep="")
orderedsizevec <- c("darkblue","lightblue","lightblue","lightblue",
                              "lightblue","blue","darkblue",
                              "darkblue","darkblue","blue","darkblue",
                              "darkblue","darkblue","blue",
                              "lightblue","darkblue","darkblue",
                              "lightblue","blue","blue","blue") #Sorry I know this is clunky, I'm just adding this last minute, I promise it's right
                              
ggplot(tbeets) +
  geom_bar( aes(x=fct_inorder(rownames(tbeets)), y=(tbeets$abdunance)), fill=orderedsizevec, stat="identity", alpha=0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Sample Abundance by Morphospecies") + 
  ylab("Abudndance (absolute)") + 
  xlab("Morphospecies ID")

ggplot(tbeets) +
  geom_bar( aes(x=fct_inorder(rownames(tbeets)), y=log(tbeets$abdunance+1)), fill=orderedsizevec,stat="identity", alpha=0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Sample Abundance by Morphospecies") + 
  ylab("Abdundance (ln(x+1) transformed)") + 
  xlab("Morphospecies ID")

#Figure 2: Total abundance per plot
beets.samps <- beets
beets.samps$samplesize <- rowSums(beets.samps)
beets.samps$foresttype <-  c(rep("Korean pine plantation",4),rep("Larch plantation",4),
                             rep("Mixed coniferous",4), rep("Mixed broadleaf",4),
                             rep("Korean pine - broadleaf mixed",4),rep("Birch forest",4),
                             "Larch plantation","Korean pine - broadleaf mixed",
                             "Mixed broadleaf","Mixed coniferous")
beets.samps <- beets.samps[c(1:8,25,9:12,28,13:16,27,17:20,26,21:24),]
beets.reord <- beets[c(1:8,25,9:12,28,13:16,27,17:20,26,21:24),]

ggplot(beets.samps) +
  geom_bar( aes(x=fct_inorder(rownames(beets.samps)), y=beets.samps$samplesize, fill=foresttype), stat="identity", alpha=0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Sample Abundance by Plot") + 
  xlab("Plot name") + 
  ylab("Abundance") +
  scale_fill_discrete(breaks=c("Korean pine plantation", "Larch plantation",
                               "Mixed coniferous", "Mixed broadleaf",
                               "Korean pine - broadleaf mixed","Birch forest"),
                      name = "Forest type")



#Figure 3: Boxplots of species richness, Simpsons index, Chao1
beets.samps$simpson <- 0
for(i in 1:nrow(beets.samps)){
  beets.samps$simpson[i] <- abdiv::simpson(beets[i,])
}
beets.samps.4s <- beets.samps[-c(9,14,19,24),]

par(mfrow = c(1, 1))

ggplot(beets.samps.4s, aes(x=fct_inorder(foresttype), y=simpson, fill=foresttype)) + 
  geom_boxplot(show.legend=FALSE) + 
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Simpson's Index") + 
  xlab("Forest type") + 
  ylab("Simpson's Index") +
  scale_fill_discrete(breaks=c("Korean pine plantation", "Larch plantation",
                               "Mixed coniferous", "Mixed broadleaf",
                               "Korean pine - broadleaf mixed","Birch forest"),
                      name = "Forest type")

simpsons.anova <- aov(simpson ~ foresttype, data=beets.samps.4s)
summary(simpsons.anova) #not significant

beets.samps$richness <- rowSums((beets.reord != 0))
beets.samps.4s <- beets.samps[-c(9,14,19,24),]

ggplot(beets.samps.4s, aes(x=fct_inorder(foresttype), y=richness, fill=foresttype)) + 
  geom_boxplot(show.legend=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(c(0,13)) + 
  ggtitle("Species Richness") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Forest type") + 
  ylab("Species Richness") +
  scale_fill_discrete(breaks=c("Korean pine plantation", "Larch plantation",
                               "Mixed coniferous", "Mixed broadleaf",
                               "Korean pine - broadleaf mixed","Birch forest"),
                      name = "Forest type")


richness.anova <- aov(richness ~ foresttype, data=beets.samps.4s)
summary(richness.anova) #significant to .05

beets.samps$chao1 <- getChaos(as.matrix(beets.reord))
beets.samps.4s <- beets.samps[-c(9,14,19,24),]

ggplot(beets.samps.4s, aes(x=fct_inorder(foresttype), y=chao1, fill=foresttype)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("bc-Chao1 Estimator") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Forest type") + 
  ylab("bc-Chao1 Estimator") +
  scale_fill_discrete(breaks=c("Korean pine plantation", "Larch plantation",
                               "Mixed coniferous", "Mixed broadleaf",
                               "Korean pine - broadleaf mixed","Birch forest"),
                      name = "Forest type")

chao1.anova <- aov(chao1 ~ foresttype, data=beets.samps.4s)
summary(chao1.anova) #not significant

#Figure 4: TESab for each forest type combined
all_PKP <- colSums(beets[1:4,])
all_LGP <- colSums(beets[c(5:8,25),])
all_MC <- colSums(beets[c(9:12,28),])
all_MB <- colSums(beets[c(13:16,27),])
all_KPBM <- colSums(beets[c(17:20,26),])
all_BF <- colSums(beets[c(21:24),])

combtypes <- rbind(all_PKP,all_LGP,all_MC,
                   all_MB,all_KPBM,all_BF)

TESabs <- TESify(combtypes)
tesifybeets <- TESify(beets)
TESabs$foresttype <- c("Korean pine plantation", "Larch plantation",
                       "Mixed coniferous", "Mixed broadleaf",
                       "Korean pine - \nbroadleaf mixed","Birch forest")

ggplot(TESabs, aes(x=fct_inorder(foresttype), y=TESab.est, color=foresttype)) + 
  geom_point(stat="identity", 
             position=position_dodge(width=1), size=3) +
  geom_errorbar(aes(ymin=TESab.est-TESab.sd, ymax=TESab.est+TESab.sd), width=.2,
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Total Expected Species (TESab) Estimates") + 
  xlab("Forest type") + 
  ylab("TESab Estimate") + 
  scale_color_discrete(breaks=TESabs$foresttype,
                       name = "Forest type") +
  ylim(c(0,35.5))


#Figure 5: Actual species shared NMDS
beets.4s <- beets[-c(25:28),]
spec <- c("#00BA38","#00BFC4","#F564E3",
                   "#619CFF","#B79F00","#F8766D")
cols <- c(rep(spec[1],4),rep(spec[2],4),rep(spec[3],4),
        rep(spec[4],4),rep(spec[5],4),rep(spec[6],4))
ftype <- c(rep("PKP",4),rep("LGP",4),rep("MC",4),
        rep("MB",4),rep("KPBM",4),rep("BF",4))

beetles.spp_NMS <- metaMDS(beets.4s,
                           distance = "jaccard",
                           k = 3,
                           maxit = 999, 
                           trymax = 500,
                           wascores = TRUE)

ordiplot(beetles.spp_NMS,type="n", main="NMDS: Jaccard Dissimilarity (actual shared species)",
         xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
shapes <- c(17, 8, 0, 15, 18, 1)
group_shapes <- shapes[as.factor(ftype)]
points(beetles.spp_NMS, display="sites", col=cols, pch=group_shapes, cex=1.5)
ordihull(beetles.spp_NMS, groups=ftype, draw="polygon", col=c("#F8766D","#B79F00","#00BFC4",
                                                                       "#619CFF","#F564E3","#00BA38"), label=FALSE)
                                                                       
legend(x=1.5,y=1, legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                           "Mixed broadleaf","Korean pine - mixed broadleaf",
                           "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=1)

#Figure 6: CNESS, m=1 shared species
essprod <- ess(beets.4s, m = 1, index = "CNESS")

beetles.cnessm1_NMS <- metaMDS(essprod,
                               distance = "",
                               k = 2,
                               maxit = 999, 
                               trymax = 500,
                               wascores = TRUE)
ordiplot(beetles.cnessm1_NMS,type="n", main=expression(paste("CNESS, ", italic("m"), " = 1")))
shapes <- c(17, 8, 0, 15, 18, 1)
group_shapes <- shapes[as.factor(ftype)]
points(beetles.cnessm1_NMS, display="sites", col=cols, pch=group_shapes, cex=1.5)
ordihull(beetles.cnessm1_NMS, groups=ftype, draw="polygon", col=c("#F8766D","#B79F00","#00BFC4",
                                                                           "#619CFF","#F564E3","#00BA38"), label=FALSE)
                                                                           
legend("topright", legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                            "Mixed broadleaf","Korean pine - mixed broadleaf",
                            "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=1)


#Figure 6: CNESS, m=37 shared species
essprod <- ess(beets.4s, m = 37, index = "CNESS")

beetles.cnessm37_NMS <- metaMDS(essprod,
                                distance = "",
                                k = 2,
                                maxit = 999, 
                                trymax = 500,
                                wascores = TRUE)
ordiplot(beetles.cnessm37_NMS,type="n", main=expression(paste("NMDS: CNESS, ", italic("m"), " = 37")),
         xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
shapes <- c(17, 8, 0, 15, 18, 1)
group_shapes <- shapes[as.factor(ftype)]
points(beetles.cnessm37_NMS, display="sites", col=cols, pch=group_shapes, cex=1.5)
ordihull(beetles.cnessm37_NMS, groups=ftype, draw="polygon", col=c("#F8766D","#B79F00","#00BFC4",
                                                                            "#619CFF","#F564E3","#00BA38"), label=FALSE)
                                                                            
legend("topright", legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                            "Mixed broadleaf","Korean pine - mixed broadleaf",
                            "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=.8)


#Figure 7: TESS ordination
goodTESSmatrix <- getWorkingTESSdist(beets.4s)
goodBeets <- beets.4s[rownames(beets.4s) %in% rownames(goodTESSmatrix), ]
goodBeetsSamps <- TESify(goodBeets)
TESSdissimilarities <- goodTESSmatrix
for(i in 2:nrow(TESSdissimilarities)){
  for(j in 1:(i-1)){
    plot1 <- rownames(goodTESSmatrix)[i]
    plot2 <- rownames(goodTESSmatrix)[j]
    plot1TESab <- goodBeetsSamps$TESab.est[rownames(goodBeetsSamps)==plot1]
    plot2TESab <- goodBeetsSamps$TESab.est[rownames(goodBeetsSamps)==plot2]
    TESS1and2 <- goodTESSmatrix[i,j]
    TESSdissimilarities[i,j] <- TESS1and2 / (plot1TESab + plot2TESab - TESS1and2)
  }
}
#This doesn't make any sense - the TESS is greater than TESab for many plots,
#which is impossible. Highlights how the data is not suited to the method.

#let's see anyway - take out the 2 that cause problems here
fTd <- TESSdissimilarities[c(2,4:5,7:8),c(2,4:5,7:8)]
goodBeets2 <- beets.4s[rownames(beets.4s) %in% rownames(fTd), ]
templatemat <- ess(goodBeets2)
TESStoplot <- c(fTd[2:5,1],fTd[3:5,2],fTd[4:5,3],fTd[5,4])
for(i in 1:length(TESStoplot)){
  templatemat[i] <- TESStoplot[i]
}

TESS.NMDS <- metaMDS(templatemat,
                     k = 2,
                     maxit = 999, 
                     trymax = 500,
                     wascores = TRUE)
ordiplot(TESS.NMDS,type="n", main=expression(paste("NMDS: TESS")),
         xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
#text(TESS.NMDS, display="sites")
shapes2 <- c(0, 18, 15, 8, 17)
cols2 <- c("#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D")
points(TESS.NMDS, display="sites", col=cols2, pch=shapes2, cex=1.5)

#not really meaningful

#Figure 8: plants pca
setwd("~/Desktop/UCL Thesis/plant data")
all_herbs <- read.csv("all_herbs.csv")
all_herbs[is.na(all_herbs)] <- 0
all_shrubs <- read.csv("all_shrubs.csv")
all_shrubs[is.na(all_shrubs)] <- 0
all_trees <- read.csv("all_trees.csv")
all_trees[is.na(all_trees)] <- 0
all_plants <- cbind(all_trees, all_shrubs, all_herbs)
plotnames <- all_plants[,1]
all_plants <- all_plants[,-c(1,30,49)] #take out the columns called "Plot"
rownames(all_plants) <- plotnames

all_herbs <- all_herbs[,-1]
all_shrubs <- all_shrubs[,-1]
all_trees <- all_trees[,-1]


for(i in 1:nrow(all_plants)){
  for(j in 1:ncol(all_plants)){
    all_plants[i,j] <- as.numeric(all_plants[i,j])
  }
}


plant.dists <- vegdist(all_plants, method = "jaccard")
herbs.dists <- vegdist(all_herbs, method= "jaccard")
shrubs.dists <- vegdist(all_shrubs, method= "jaccard")
trees.dists <- vegdist(all_trees, method= "jaccard")


plants.pca <- wcmdscale(plant.dists,
                        #distance = "jaccard",
                        eig=TRUE
                        #maxit = 999, 
                        #trymax = 500,
                        #wascores = TRUE,
                        #w=1)
)

cols.plants <- c(rep(spec[1],4),rep(spec[2],5),rep(spec[3],5),
                 rep(spec[4],5),rep(spec[5],5),rep(spec[6],4))
ftype.plants=c(rep("PKP",4),rep("LGP",5),rep("MC",5),
               rep("MB",5),rep("KPBM",5),rep("BF",4))

ordiplot(plants.pca,type="n", main="PCoA: Plant Jaccard dissimilarity \n(actual shared plant species)",
         xlab = "PCoA Axis 1", ylab = "PCoA Axis 2")
shapes <- c(17, 8, 0, 15, 18, 1)
group_shapes.plants <- shapes[as.factor(ftype.plants)]
points(plants.pca$points[,1:2], col=cols.plants, pch=group_shapes.plants, cex=1.5)
#points(plants.pca, display="species", col="red", pch=4, cex=1)
ordihull(plants.pca, groups=ftype.plants, draw="polygon", col=c("#F8766D","#B79F00","#00BFC4",
                                                                         "#619CFF","#F564E3","#00BA38"), label=FALSE)
                                                                         
legend("topright", legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                            "Mixed broadleaf","Korean pine - mixed broadleaf",
                            "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=1)

#Figure out how many components explain 
sum(plants.pca$eig[1:4])/sum((plants.pca$eig))
#trees: first 7 explain 74.3% of variance; 5 explain 61.65%
#shrubs: first 6 explain 76.44% of variance; 4 explain 62.63%
#herbs: first 7 explain 76.63% of variance; 5 explain 65.69%; 4 explain 59.43%

#without the abs
sum(plants.pca$eig[1:6])/sum((plants.pca$eig))
#trees: first 4 explain 62.80
#shrubs: first 3 explain 62.83, first 4 explain 73.26
#herbs: first 4 explain 60.46


#all plants
#first 4 explain 48.9
#first 5 explain 55.93
#first 6 explain 61.97

#Figure 9: Plants and beetles CCA
PCvals <- plants.pca$points
PCvals.4s <- PCvals[-c(9,14,19,24),]
beets.4s.log <- log(beets.4s+1)

beets.cca <- cca(beets.4s.log ~ (PCvals.4s[,1:4]))
ordiplot(beets.cca, type = "n",main="CCA Ordination\nFour all-plant PCoA axes as constraints",
         xlab = "CCA Axis 1", ylab = "CCA Axis 2")
constraints <- scores(beets.cca, display = "bp")
#constraints[4,1] <- -1*constraints[4,1] #flipping this one around for legibility
arrows(0, 0, constraints[, 1], constraints[, 2],
       col = "blue", length = 0.1, angle = 15)
text(constraints[,1]*1.2, constraints[,2]*1.2,col="blue",labels = c("PCoA 1",
                                                                    "PCoA 2",
                                                                    "PCoA 3",
                                                                    "PCoA4"),cex=0.75)
points(beets.cca, display="sites", col=cols, pch=group_shapes, cex=1)
#points(beets.cca, display="species",col="red",pch=4,cex=1)
text(beets.cca, display="species", col="red", labels = c("ms1","ms2","ms3","ms4",
                                                         "ms5","ms6","ms7","ms8",
                                                         "ms9","ms10","ms11","ms12",
                                                         "ms13","ms14","ms15","ms16",
                                                         "ms17","ms18","ms19","ms20","ms21"), cex=0.5)

ordihull(beets.cca, groups=ftype, draw="polygon", col=c("#F8766D","#B79F00","#00BFC4",
                                                                 "#619CFF","#F564E3","#00BA38"), label=FALSE)
legend("topright", legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
"Mixed broadleaf","Korean pine - mixed broadleaf",
"Birch forest"), 
pch=c(1, 0, 18, 15, 8, 17),
col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
title="Forest type", cex=1)

#trees pcoa explains 39.65% of beetle variation with 7 components
#with 4 it's 21.21%
#with 5 it's 25.70%
#shrubs explain 27.35% of beetle variation with 6 components
#with 3 it's 14.99%
#with 4 it's 18.28%
#with 5 it's 22.19%
#herbs explain 39.51% of beetle variation with 7 components
#with 4 it's 25.84% 
#with 5 it's 28.95%

#all plants: with 6 it's 35.68
#with 5 it's 32.33%
#with 4 it's 28.41%

#Figure 10: Correlations: richness and TESb over plant richness
allTESresults <- TESify(beets)
allTESresults <- allTESresults[c(1:8,25,9:12,28,13:16,27,17:20,26,21:24),]
allTESresults.4s <- allTESresults[-c(9,14,19,24),]


plantsAlpha <- rowSums(all_plants)
plantsAlpha.4s <- rowSums(all_plants)[-c(9,14,19,24)]

lm.r= lm(formula = allTESresults$alpha ~ plantsAlpha)
summary(lm.r) #no significance

par(mar = c(5, 5, 5, 5))
plot(plantsAlpha.4s,allTESresults.4s$alpha,main="A) Species Richness",
     col = cols,pch=group_shapes,
     ylab="Beetle species richness",
     xlab="Plant species richness")#no pattern
legend("topright",legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                           "Mixed broadleaf","Korean pine - mixed broadleaf",
                           "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=1, ncol=1,
       xpd = TRUE, inset = c(-.6, .3))


lm.r= lm(formula = allTESresults$TESb.est ~ plantsAlpha)
summary(lm.r) #no significance

par(mar = c(5, 5, 5, 10))
plot(plantsAlpha.4s,allTESresults.4s$TESb.est,main="B) TESb estimate",
     col = cols,pch=group_shapes,
     ylab="Beetle TESb estimate",
     xlab="Plant species richness")#no pattern
legend("topright",legend=c("Korean pine plantation","Larch plantation","Mixed coniferous",
                           "Mixed broadleaf","Korean pine - mixed broadleaf",
                           "Birch forest"), 
       pch=c(1, 0, 18, 15, 8, 17),
       col=c("#00BA38", "#00BFC4","#F564E3","#619CFF","#B79F00","#F8766D"), 
       title="Forest type", cex=1, ncol=1,
       xpd = TRUE, inset = c(-.6, .3))



#Figure 11: t-tests of similarity values


dists <- vegdist(beets.4s, method = "jaccard")
dists <- ess(beets.4s,m=37,index="CNESS")

#alt: dists <- ess(beets.4s, m=37, index="CNESS")
distmat <- makeTESSMatrix(dists,rownames(beets.4s))

PKPins <- c(distmat[substr(rownames(distmat),start=1,stop=3)=="PKP",
                    substr(colnames(distmat),start=1,stop=3)=="PKP"])
PKPins <- PKPins[PKPins!=0]
PKPouts <- c(distmat[substr(rownames(distmat),start=1,stop=3)!="PKP",
                     substr(colnames(distmat),start=1,stop=3)=="PKP"],
             distmat[substr(rownames(distmat),start=1,stop=3)=="PKP",
                     substr(colnames(distmat),start=1,stop=3)!="PKP"])
PKPouts <- PKPouts[PKPouts!=0]

t.test(PKPins,PKPouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=PKPins),fill = "#00BA38") + 
  geom_boxplot(aes(x= "Without-type", y=PKPouts),fill = "#00BA38") +
  ggtitle("A) Korean pine plantation") + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") + 
  ylim(c(0,1))

LGPins <- c(distmat[substr(rownames(distmat),start=1,stop=3)=="LGP",
                    substr(colnames(distmat),start=1,stop=3)=="LGP"])
LGPins <- LGPins[LGPins!=0]
LGPouts <- c(distmat[substr(rownames(distmat),start=1,stop=3)!="LGP",
                     substr(colnames(distmat),start=1,stop=3)=="LGP"],
             distmat[substr(rownames(distmat),start=1,stop=3)=="LGP",
                     substr(colnames(distmat),start=1,stop=3)!="LGP"])
LGPouts <- LGPouts[LGPouts!=0]

t.test(LGPins,LGPouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=LGPins),fill = "#00BFC4") + 
  geom_boxplot(aes(x= "Without-type", y=LGPouts),fill = "#00BFC4") +
  ggtitle(expression(bold("B) Larch plantation ***"))) + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") +
  ylim(c(0,1))


MCins <- c(distmat[substr(rownames(distmat),start=1,stop=2)=="MC",
                   substr(colnames(distmat),start=1,stop=2)=="MC"])
MCins <- MCins[MCins!=0]
MCouts <- c(distmat[substr(rownames(distmat),start=1,stop=2)!="MC",
                    substr(colnames(distmat),start=1,stop=2)=="MC"],
            distmat[substr(rownames(distmat),start=1,stop=2)=="MC",
                    substr(colnames(distmat),start=1,stop=2)!="MC"])
MCouts <- MCouts[MCouts!=0]

t.test(MCins,MCouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=MCins),fill = "#F564E3") + 
  geom_boxplot(aes(x= "Without-type", y=MCouts),fill = "#F564E3") +
  ggtitle("C) Mixed coniferous") + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") +
  ylim(c(0,1))


MBins <- c(distmat[substr(rownames(distmat),start=1,stop=2)=="MB",
                   substr(colnames(distmat),start=1,stop=2)=="MB"])
MBins <- MBins[MBins!=0]
MBouts <- c(distmat[substr(rownames(distmat),start=1,stop=2)!="MB",
                    substr(colnames(distmat),start=1,stop=2)=="MB"],
            distmat[substr(rownames(distmat),start=1,stop=2)=="MB",
                    substr(colnames(distmat),start=1,stop=2)!="MB"])
MBouts <- MBouts[MBouts!=0]

t.test(MBins,MBouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=MBins),fill = "#619CFF") + 
  geom_boxplot(aes(x= "Without-type", y=MBouts),fill = "#619CFF") +
  ggtitle("D) Mixed broadleaf") + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") +
  ylim(c(0,1))

KPBMins <- c(distmat[substr(rownames(distmat),start=1,stop=4)=="KPBM",
                     substr(colnames(distmat),start=1,stop=4)=="KPBM"])
KPBMins <- KPBMins[KPBMins!=0]
KPBMouts <- c(distmat[substr(rownames(distmat),start=1,stop=4)!="KPBM",
                      substr(colnames(distmat),start=1,stop=4)=="KPBM"],
              distmat[substr(rownames(distmat),start=1,stop=4)=="KPBM",
                      substr(colnames(distmat),start=1,stop=4)!="KPBM"])
KPBMouts <- KPBMouts[KPBMouts!=0]

t.test(KPBMins,KPBMouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=KPBMins),fill = "#B79F00") + 
  geom_boxplot(aes(x= "Without-type", y=KPBMouts),fill = "#B79F00") +
  ggtitle("E) Korean pine - broadleaf mixed") + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") +
  ylim(c(0,1))


BFins <- c(distmat[substr(rownames(distmat),start=1,stop=2)=="BF",
                   substr(colnames(distmat),start=1,stop=2)=="BF"])
BFins <- BFins[BFins!=0]
BFouts <- c(distmat[substr(rownames(distmat),start=1,stop=2)!="BF",
                    substr(colnames(distmat),start=1,stop=2)=="BF"],
            distmat[substr(rownames(distmat),start=1,stop=2)=="BF",
                    substr(colnames(distmat),start=1,stop=2)!="BF"])
BFouts <- BFouts[BFouts!=0]

t.test(BFins,BFouts)

ggplot() + 
  geom_boxplot(aes(x="Within-type", y=BFins),fill = "#F8766D") + 
  geom_boxplot(aes(x= "Without-type", y=BFouts),fill = "#F8766D") +
  ggtitle("F) Birch forest") + 
  ylab(expression(paste("CNESS, ",italic("m")," = 37"))) + 
  xlab("") + 
  ylim(c(0,1))

#moral of the story: significant differences only in the Larch plantation,
#and only when CNESS is used rather than Jaccard


#Figure 12: ANOVA of 7 dissimilarities
Betweens<- unique(c(PKPouts,LGPouts,MCouts,MBouts,KPBMouts,BFouts)) #this works because there are no repeated values
anovInput <- data.frame(
  Diss=c(PKPins, LGPins, MCins, MBins, KPBMins, BFins, Betweens),
  Ftype =factor(rep(c("PKP", "LGP", "MC", "MB","KPBM","BF","Between-type"),
                    times=c(length(PKPins), length(LGPins), length(MCins), 
                            length(MBins), length(KPBMins), length(BFins),
                            length(Betweens))))
)

anova.ftype <- aov(Diss~Ftype, data=anovInput)
anova(anova.ftype)

ggplot(anova.ftype, aes(x=fct_inorder(Ftype), y=Diss)) + 
  geom_boxplot() 

par(mar = c(5, 10, 5, 5))
tuks <- TukeyHSD(anova.ftype, conf.level=.95)
plot(tuks, las=2)


#no differences with Jaccard
#CNESS shows LGP more similar than the rest, difference between LGP and KPBM




