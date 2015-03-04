### Bern - 03Feb15 - ESM567 - pg. 66 in textbook, 512 in Ecology

###---- Assumptions and what to report ----
### Nonmetric Multidimensional Scaling (NMDS) with Stable Solution from Random Starts, 
### Axis Scaling and Species Scores . 

### Adequacy is dependant on: 
### Stress value < 10, check Shepard plot for outliers, Check distortion: select top 10% of values 
### in Bray-Curtis matrix and compare to the NMDS ordination distance
### What to report:
### Distance measure used (BC), Algorithm/Software used, Number of runs, Number of dimensions (2) in the final solution
### Interpretation of each axis

### ---- Data ----
dta <-read.csv("apoch.csv")
apo <- dta[-c(32:33),] # only two sites for land form R
head(apo)
Y<-as.matrix(apo[,c(2:8)]) # multivariate response variable matrix (morphology)
env <-apo[,c(9:14)] # Environemtal group

### -----NMDS Model-----
library(vegan)
library(MASS)
mod <- metaMDS(Y,k=2,trace=T,distance="bray",autotransform=F,trymax = 20) #MDS with 2 dimension
names(mod)# what type of informatin is produced from NMDS.

br<-vegdist(Y, method="bray") # calculate Bray-Curtis distance, or other dependant on data type
stressplot(mod,br,p.col="lemonchiffon3",l.col="darkolivegreen",lwd=2) 
legend("bottomright",c("Stress = 0.191"))
stressplot(mod,br,pch=21,bg="beige",p.col="lemonchiffon3",l.col="darkolivegreen")  # Generate a Sharperd diagram by plotting BC between sites against Eucledian distance between sites in NMDS 
legend("bottomright",c("Stress = 0.191"))

mod$stress # Check for MDS' stress value
round(mod$species,2) # Variables contribution to MDS axis

cor(mod$points[,1],apo[,c(2:8)]) # Variable correlation with axis MDS1
cor(mod$points[,2],apo[,c(2:8)]) # Variable correlation with axis MDS2

par(mfrow=c(2,2))
plot(mod$points[,1],apo[,3],xlab="MDS 1",ylab="chela") # Plots original variable against MDS.
plot(mod$points[,1],apo[,2],xlab="MDS 1",ylab="fem") # Only variables with r > .50
plot(mod$points[,2],apo[,2],xlab="MDS 2",ylab="fem")
plot(mod$points[,2],apo[,3],xlab="MDS 2",ylab="chela")

### Source cor.matrix, for NMDS 1 & 2 correlated with original variables
mrph <- apo[,c(2:8)]
MDS1 <- mod$points[,1]
mrphMDS1 <-cbind(MDS1,mrph)
cor.matrix(mrphMDS1)

MDS2 <- mod$points[,2]
mrphMDS2 <-cbind(MDS2,mrph)
cor.matrix(mrphMDS2)

####~~~~ Sydney's Plots ~~~~####
par(mfrow=c(1,2))
plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$sex),col=c("red","blue")[apo$sex], xlab="NMDS 1", ylab="NMDS 2")
legend("topright",c("F","M"),pch=c(1,2),col=c("red","blue"),inset=.02) 
title(main="Chela gradient - Axis 1")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

plot(mod$points[,1], mod$points[,2], pch=as.numeric(apo$form),col=c(1,2,3,4)[apo$form], xlab="NMDS 1", ylab="NMDS 2") 
legend("topright",c("A","C","T"),pch=c(1,2,4),col=c(1,2,4),inset=.02) 
title(main="Femor gradient - Axis 2")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

plot(apo$LAT, apo$LONG,type='n', xlab="Latitude", ylab="Longitude")
symbols(apo$LAT, apo$LONG,circles=apo$eco4, inches=0.2, add=T, lwd=2)  
title(main="Sample site dispersion")


### hypothesis test (with strata)
adonis(Y ~ form*sex*eco4, data=env, perm=999) # Eco4 is similar to land form

apo.dist <-vegdist(Y)
attach(env)
apoS.ano <-anosim(apo.dist,sex)
apoF.ano <-anosim(apo.dist,form)
summary(apoS.ano)
summary(apoF.ano)
plot(apoS.ano)
plot(apoF.ano)

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




