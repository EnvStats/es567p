### Bern - 03Feb15 - ESM567 - pg. 66 in textbook, 512 in Ecology

###---- Assumptions and what to report ----
### Nonmetric Multidimensional Scaling (NMDS) with Stable Solution from Random Starts, 
### Axis Scaling and Species Scores . 

### Adequacy is dependant on: 
### Stress value < .10, check Shepard plot for outliers, Check distortion: select top 10% of values 
### in Bray-Curtis matrix and compare to the NMDS ordination distance
### What to report:
### Distance measure used (BC), Algorithm/Software used, Number of runs, Number of dimensions (2) in the final solution
### Interpretation of each axis

### ---- Data ----
dta <- read.csv("apoch.csv")
#apoch <- droplevels(subset(dta,eco4!="WCLV" & eco4!="WKF" & eco4!="PVB" & form!="R"))
apoch <- droplevels(subset(dta,form!="R"))
apo <-apoch[-c(70),] #outlier

ap1 <- apo[,c(2:8)]# multivariate response variable matrix (morphology)
ap <-as.data.frame(scale(ap1, center=F, scale=T)) # total lenght is not a ratio like other variables
env <- apo[,c(9:14)] # Environemtal groups
summary(apo)
boxplot(ap)
cor(ap)
cor.matrix(ap)

### -----nMDS Model-----
library(vegan)
library(MASS)
## Run multiple times to reduce the risk of getting stuck in local min, trace=T so you can see
## Euclidean distances of a matrix where columns are centred, have unit variance, and are uncorrelated
mod <- metaMDS(ap,k=3,trace=T,distance="euclidian",autotransform=F,trymax = 20) #nMDS with k=2 dimensions
names(mod)# Type of informatin produced from nMDS.
round(mod$species,2) # Variables contribution to nMDS axis

par(mfrow=c(1,2)) #Stress plot and nMDS plot
br<-vegdist(ap, method="bray") # calculate Bray-Curtis distance, or other dependant on data type
stressplot(mod,br,p.col="lemonchiffon3",l.col="darkolivegreen",lwd=2, main="nMDS/Euc.")
legend("bottomright",c(paste("Stress =",round(mod$stress*100,2))))
plot(mod,type="t", main=paste("nMDS Dim 1 vs Dim 2"))

#Variable Correlation with MDS1 & 2 Axis.
cor(mod$points[,1],apo[,c(2:8)]) # Variable correlation with axis MDS1
cor(mod$points[,2],apo[,c(2:8)]) # Variable correlation with axis MDS2

par(mfrow=c(2,2)) # Correlation plots
plot(mod$points[,1],ap1[,3],xlab="nMDS 1",ylab="tib", col="blue") # Plots original variable against MDS.
legend("topright",c("r = -0.89"))
plot(mod$points[,2],ap1[,5],xlab="nMDS 2",ylab="abd", col="red")
legend("topright",c("r = -0.80"))
plot(mod$points[,2],ap1[,4],xlab="nMDS 2",ylab="carp", col="green")
legend("topright",c("r = -0.61"))


### Source cor.matrix, for NMDS 1 & 2 correlated with original variables
MDS <- mod$points
apMDS <-cbind(MDS,ap)
cor.matrix(apMDS)

####~~~~ Sydney's Plots ~~~~####
par(mfrow=c(2,2))
plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$sex),col=c("red","darkgreen")[apo$sex], xlab="nMDS 1", ylab="nMDS 2")
legend("topright",c("F","M"),pch=c(1,2),col=c("red","darkgreen"),inset=.02) 
title(main="Sex - Dim 1 vs Dim 2")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
text(mod, "species", col="grey", cex=0.8)

plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$eco4),col=c(1,2,3,"grey",5)[apo$eco4], xlab="nMDS 1", ylab="nMDS 2") 
legend("topright",c("CPL","NFRF","OCF","V","WH"),pch=c(1,2,3,4,5),col=c(1,2,3,"grey",5),inset=.02) 
title(main="Eco4 - Dim 1 vs Dim 2")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$form),col=c(2,"lightgrey",4)[apo$form], xlab="nMDS 1", ylab="nMDS 2") 
legend("topright",c("A","C","T"),pch=c(1,2,3),col=c(2,"lightgrey",4),inset=.02) 
title(main="Form - Dim 1 vs Dim 2")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

### distance plot
plot(apo$LAT, apo$LONG,type='n', xlab="Latitude", ylab="Longitude")
symbols(apo$LAT, apo$LONG,circles=apo$eco4, inches=0.2, add=T, lwd=2)  
title(main="Sample site dispersion")

### ~~ Hypothesis Tests ~~ ####

### ANOSIM-test statistically whether there is a significant difference between 
### two or more groups of sampling units.

# Mahalanobis distances are Euclidean distances of a matrix where columns are centred, 
## have unit variance, and are uncorrelated. The index is not commonly used for community data, 
## but it is sometimes used for environmental variables. The calculation is based on transforming 
## data matrix and then using Euclidean distances following Mardia et al. (1979).

## If two groups of sampling units are really different in their species composition, 
## then compositional dissimilarities between the groups ought to be greater than those within the groups.
## The anosim statistic R is based on the difference of mean ranks between groups (r_B) and within groups (r_W):
ap.dist <-vegdist(ap, "euclidean")
attach(env)
apoS.ano <-anosim(ap.dist,sex)
apoF.ano <-anosim(ap.dist,form)
apoE.ano <-anosim(ap.dist,eco4)
summary(apoS.ano)
summary(apoF.ano)
summary(apoE.ano)
par(mfrow=c(2,2)) # For there to be a sig dif, the between disimilarity needs to be greater than the within groups
plot(apoS.ano, title(main="Gender-Dim 1"), col = "lightgray") # Plot width is disimalarity within group
plot(apoF.ano, title(main="Land form-Dim 2"), col = "blue")
plot(apoE.ano, title(main="Eco4-Dim2"), col = "darkgreen")
detach(env)

### Pairwise t-test and aov for differences between groups
pairwise.t.test(mod$points[,1], env$sex, p.adj = "bonferroni") 
pairwise.t.test(mod$points[,2], env$form, p.adj = "bonferroni")
pairwise.t.test(mod$points[,2], env$eco4, p.adj = "bonferroni")

TukeyHSD(aov(mod$points[,1] ~ env$sex, ap)) 
TukeyHSD(aov(mod$points[,2] ~ env$form, ap))
TukeyHSD(aov(mod$points[,2] ~ env$eco4, ap))

### Optional adonis ANOVA
adonis(ap ~ sex*form*eco4, data=env, permutations=999) #Original data

### ---- Optional Spider Plots ----
par(mfrow=c(1,2))
plot(mod, type="t") # gender(sex)
### Hulls show treatment
#with(env, ordihull(mod, group=sex, show="m"))
#with(env, ordihull(mod, group=sex, show="f", col="green")) 
### Spider shows fields
with(env, ordispider(mod, group=sex, lty=5, col="blue"))
title(main="Gender - Axis 1")

plot(mod, type="t") # land Form
### Hulls show treatment
#with(env, ordihull(mod, group=form, show="T"))
#with(env, ordihull(mod, group=form, show="A", col="green"))
#with(env, ordihull(mod, group=form, show="C", col="red")) 
### Spider shows fields
with(env, ordispider(mod, group=form, lty=4, col="darkgreen"))
title(main="Land Form - Axis 2")

### ---- example of strata, nested (block) design ----
plot(mod, type="t") # Eco4, marginal p-value
### Hulls show treatment, within field
with(env, ordihull(mod, group=form, show="A"))
with(env, ordihull(mod, group=form, show="C", col="green")) #group = strata
with(env, ordihull(mod, group=form, show="T", col="blue"))
### Spider shows fields
with(env, ordispider(mod, group=sex, lty=5, col="blue"))

