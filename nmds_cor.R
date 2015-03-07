apo <-read.csv("apoch.csv")
ap <- apo[-c(32:33),] #removed becuase only two land form for R
head(ap)
Y<-as.matrix(ap[,c(2:8)]) #multivariate response variable matrix (morphology)
env <-ap[,c(9:12)] #Environemtal group
x <-ap[,c(12)]

names(mod)# what type of informatin is produced from NMDS.
mod$points

plot(mod$points[,1],mod$points[,2]) # NMDS1 vs NMDS2

cor(mod$points[,1],ap[,c(2:8)]) #Is there a correlation with NMDS axis one?
cor(mod$points[,2],ap[,c(3:8)]) # Is there a correlation with NMDS axis two?

plot(mod$points[,1],ap[,3]) #Plots original variable  Chela against NMDS.
symbols(mod$points[,1],mod$points[,2],circles=ap[,3],inches=0.2)

plot(mod$points[,1],mod$points[,2],type='n') #Type=n leaves the plot empty for symbols
symbols(mod$points[,1],mod$points[,2],circles=ap[,3],inches=0.1)

par(mfrow=c(1,2))
plot(mod$points[,1], mod$points[,2], pch=as.numeric(ap$sex),col=c("red","blue")[x], xlab="NMDS 1", ylab="NMDS 2") #eigenvalues for each axis
legend("topright",c("F","M"),pch=c(1,2),col=c("red","blue"),inset=.02) 
title(main="PC1-Chela gradient")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")


