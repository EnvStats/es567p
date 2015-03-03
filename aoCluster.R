# Bern - Looking for patterns in morpological distribution

apo<-read.csv("apoch.csv")
##Data exploration##
br <-apo[,c(2:8)] #br is body ratios

#### cluster analysis ####
library(ade4)
library (vegan)  # should be loaded after ade4 to avoid some conflicts!
library (gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)

# Apo morphology data: compute matrix of chord distance among 
# sites followed by single linkage agglomerative clustering
# **************************************************************
bod.norm <- decostand(br, "normalize") #normalized data
bod.ch <- vegdist(bod.norm, "euc") #calculate 'chord' distance which is normalized Euclidean distance
bod.ch.single <- hclust(bod.ch, method="single") #single linkage cluster

# Plot a dendrogram using the default options
par(mfrow=c(1,1))
plot(bod.ch.single)

# Compute complete-linkage agglomerative clustering
# **************************************************
bod.ch.complete <- hclust(bod.ch, method="complete")
plot(bod.ch.complete, label=apo$sex)

# Compute UPGMA agglomerative clustering
# **************************************
par(mfrow=c(2,1))
bod.ch.UPGMA <- hclust(bod.ch, method="average")
plot(bod.ch.UPGMA)
plot(bod.ch.UPGMA,labels=apo$sex, fg=gray(0.7))

# Compute centroid clustering 
# ********************************************
bod.ch.centroid <- hclust(bod.ch, method="centroid")
plot(bod.ch.centroid, labels=apo$sex)

# Compute Ward's minimum variance clustering
# ******************************************
bod.ch.ward <- hclust(bod.ch, method="ward")
plot(bod.ch.ward, labels=apo$form)

### Cophenetic correlation
# **********************

# Single linkage clustering
bod.ch.single.coph <- cophenetic(bod.ch.single)
cor(bod.ch, bod.ch.single.coph)   ##output: .4603724
# Complete linkage clustering
bod.ch.comp.coph <- cophenetic(bod.ch.complete)
cor(bod.ch, bod.ch.comp.coph)   ##output: .4826873
# Average clustering
bod.ch.UPGMA.coph <- cophenetic(bod.ch.UPGMA)
cor(bod.ch, bod.ch.UPGMA.coph)   ##output:.5928187
# Ward clustering
bod.ch.ward.coph <- cophenetic(bod.ch.ward)
cor(bod.ch, bod.ch.ward.coph)   ##output: .4342331

#***********************************
### Shepard-like diagrams
# *********************
par(mfrow=c(2,2))
plot(bod.ch, bod.ch, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)), 
     main=c("Single linkage",paste("Cophenetic correlation ",
     round(cor(bod.ch, bod.ch.comp.coph),3))))
abline(0,1)
lines(lowess(bod.ch, bod.ch.single.coph), col="red")

plot(bod.ch, bod.ch.comp.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Complete linkage", paste("Cophenetic correlation ",
     round(cor(bod.ch, bod.ch.comp.coph),3))))
abline(0,1)
lines(lowess(bod.ch, bod.ch.comp.coph), col="red")

plot(bod.ch, bod.ch.UPGMA.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)), 
     main=c("UPGMA", paste("Cophenetic correlation ",
     round(cor(bod.ch, bod.ch.UPGMA.coph),3))))
abline(0,1)
lines(lowess(bod.ch, bod.ch.UPGMA.coph), col="red")

plot(bod.ch, bod.ch.ward.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), 
     ylim=c(0,max(bod.ch.ward$height)),
     main=c("Ward clustering", paste("Cophenetic correlation ",
     round(cor(bod.ch, bod.ch.ward.coph),3))))
abline(0,1)
lines(lowess(bod.ch, bod.ch.ward.coph), col="red")

