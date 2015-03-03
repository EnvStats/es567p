### Bern - NMDS - 03Feb15 - ESM567

###---- Assumptions and what to report ----
### Nonmetric Multidimensional Scaling with Stable Solution from Random Starts, 
### Axis Scaling and Species Scores (NMDS). 

### Adequacy is dependant on: 
### Stress value < 10, check Shepard plot for outliers, Check distortion: select top 10% of values 
### in Bray-Curtis matrix and compare to the MDS ordination distance
### What to report:
### Distance measure used (BC), Algorithm/Software used, Number of runs, Number of dimensions (2) in the final solution
### Interpretation of each axis

### ---- Data ----
apo <-read.csv("apoch.csv")
### ap <- apo[-c(32:33),] # only two sites for land form R
head(apo)
Y<-as.matrix(apo[,c(2:8)]) # multivariate response variable matrix (morphology)
env <-apo[,c(9:14)] # Environemtal group

### -----NMDS Model-----
library(vegan)
library(MASS)
mod <- metaMDS(Y,k=2,trace=T, distance="bray",autotransform=F, trymax = 20) #NMDS with 2 dimension
names(mod)# what type of informatin is produced from NMDS.

br<-vegdist(Y, method="bray") # calculate Bray-Curtis distance, or other dependant on data type
stressplot(mod,br)  # Generate a Sharperd diagram by plotting BC between sites against Eucledian distance between sites in NMDS 
mod$stress # Check for NMDS' stress value
round(mod$species,2) # Variables contribution to NMDS axis

cor(mod$points[,1],apo[,c(2:8)]) # Variable correlation with axis NMDS1
cor(mod$points[,2],apo[,c(2:8)]) # Variable correlation with axis NMDS2

par(mfrow=c(2,2))
plot(mod$points[,1],apo[,3], xlab="NMDS 1", ylab="chela") # Plots original variable against NMDS.
plot(mod$points[,1],apo[,2], xlab="NMDS 1", ylab="fem") # Only variables with r > .50
plot(mod$points[,2],apo[,2], xlab="NMDS 2", ylab="fem")
plot(mod$points[,2],apo[,3], xlab="NMDS 2", ylab="chela")

### Source cor.matrix, for NMDS 1 & 2 correlated with original variables
mrph <- apo[,c(2:8)]
NMDS1 <- mod$points[,1]
mrphMDS1 <-cbind(NMDS1,mrph)
cor.matrix(mrphMDS1)

NMDS2 <- mod$points[,2]
mrphMDS2 <-cbind(NMDS2,mrph)
cor.matrix(mrphMDS2)

####~~~~ Sydney Plots ~~~~####
par(mfrow=c(1,2))
plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$sex),col=c("red","blue")[apo$sex], xlab="NMDS 1", ylab="NMDS 2")
legend("topright",c("F","M"),pch=c(1,2),col=c("red","blue"),inset=.02) 
title(main="Chela gradient - Axis 1")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$form),col=c(1,2,3,4)[apo$form], xlab="NMDS 1", ylab="NMDS 2") 
legend("topright",c("A","C","R","T"),pch=c(1,2,3,4),col=c(1,2,3,4),inset=.02) 
title(main="Femor gradient - Axis 2")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

plot(apo$LAT, apo$LONG,type='n', xlab="Latitude", ylab="Longitude")
symbols(apo$LAT, apo$LONG,circles=apo$form, inches=0.2, add=T, lwd=2)  
title(main="Sample site dispersion")


### hypothesis test (with strata)
adonis(Y ~ form*sex*eco4, data=env, perm=999)


### ---- Optional Spider Plots ----
par(mfrow=c(1,2))
plot(mod, type="t") # gender(sex)
### Hulls show treatment
#with(env, ordihull(mod, group=sex, show="m"))
#with(env, ordihull(mod, group=sex, show="f", col="green")) 
### Spider shows fields
with(env, ordispider(mod, group=sex, lty=5, col="blue"),
     spiders = c("centroid", "median"),  show.groups, 
     label = FALSE)
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
with(env, ordihull(mod, group=sex, show="m"))
with(env, ordihull(mod, group=sex, show="f", col="green")) #group = strata
### Spider shows fields
with(env, ordispider(mod, group=eco4, lty=5, col="blue"))

