#Bern Romey, 04Feb15 ~ ESM567 Term Project
#PCA - MANOVA

#Goal: is there a pattern in genetic variability with morphology?

#Data- Ratio and total body length
dta <- read.csv("apoch.csv")
am <- dta[c(2:8)] # scale, tl not a ratio
env <-dta[c(9:12)]
str(am)

#Assumptions
boxplot(am, main = "Not scaled")
boxplot(scale(am), main="Scaled (centered) with Z-score")

cor.matrix(am)# source cor.matrix function

cov(am) # calculate cov matrix with the standardized data

#PCA Analysis
require(MASS) # loads the PCA package
pca <- princomp(am, cor=T) # creates a PC matrix using the cor matrix
biplot(pca, expand = 1.05,main = "Biplot", xlab = "Comp.1 (42.0%)", ylab = "Comp.2 (33.3%)")
# Scale for sites(PC matrix-pca$scores) on top, scale for variables (vectors-loadings) along bottom
summary(pca) # proportion of variance is eigenvalues for each PC

library(vegan)
screeplot(pca, bstick = TRUE, main="PCA") # inertia= variance n PCA

round(loadings(pca),2) # Check eigenvectors: length of vector is relative variance and how much it contributes to the PC
# Principal component loading (pg 50).  The further from zero, the greater the contribution.
round(loadings(pca)[,c(1:2)],2) # Loading for PC1 & 2 only

round((pca$scores),2) # PC matrix showing site scores for all PCs. How far each is(SD) from the the grand centroid
# This is the distribution of PC1 and PC2 site scores (top scale).  Each variable for each component. 
# In this case due to broken stick, PC1 and PC2

# Create shepard diagram
apo<-dist(am) # Calculate Euclidian distance among sites (multidimentional spcae). Check transformation and matrix used.
apo.1<-dist(pca$scores[,c(1,2)]) # Calculate Euclidian distance among sites in PCA space using only first 2 PCs (reduced space).
plot(apo,apo.1,main="PC=2", xlab="Distance in Multidimensional space", ylab="Distance in Reduced space") #x=euc, y=euc.1  

