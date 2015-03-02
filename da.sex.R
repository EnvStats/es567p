#Bern 22Feb15: MANOVA and Linear Discriminant Function Analysis - pg 673 in Legendre, 141 in text.
#The hypothesis: is there statistical sig-difference in morphology among the different groups?
#(sex, land type)
#Does the genetic information support this hypothesis?

#Difference from PCA: Lecture 8, slide 52

ap <-read.csv("apoch.csv")
head(ap)
x<-as.factor(ap[,12]) #independent variable (factor)
y<-as.matrix(ap[,c(2:8)]) #multivariate response variable matrix
env <-ap[,c(9:12)] #Environemtal group
morph <-ap[,c(2:8)]
cor.matrix(scale(morph)) #linear relationship among variables poor (<0.40) not a good PCA candidate

##Run PCA first to see if skewed site distribution (horshoe)
#----
pca<-princomp(morph,cor=T)#total length not same unit of measure, thus scale all to Z-score
summary(pca) #PC1 is a chela gradient, while pc2 is a fem gradient
biplot(pca, xlab="PC 1(41%)", ylab="PC 2(32%)")
library(vegan)
screeplot(pca, bstick = TRUE, main="PCA") #inertia= variance n PCA
round(loadings(pca)[,c(1:2)],2)#contribution to PC1 & 2

par(mfrow=c(1,2))
plot(pca$scores[,1], pca$scores[,2],pch=as.numeric(ap$sex),col=c("red","blue")[x], xlab="PC 1(41%)", ylab="PC 2(32%)") #eigenvalues for each axis
legend("topright",c("F","M"),pch=c(1,2),col=c("red","blue"),inset=.02) 
title(main="PC1-Chela gradient")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

# PC2 is a C T gradient
x3 <-as.factor(ap[,9])
plot(pca$scores[,1], pca$scores[,2],pch=as.numeric(ap$form),col=c(1,2,3,4)[x3],xlab="PC 1 (41%)", ylab="PC 2(32%)") #eigenvalues for each axis
legend("topright",c("A","C","R","T"),pch=c(1,2,3,4),col=c(1,2,3,4),inset=.02) 
title(main="PC2-Femor gradient")
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")

##Run one-way MANOVA
#----
mod<-manova(y~x)# Is there a sig-dif in variance?
summary.manova(mod) #MANOVA table, compare it with univariate ANOVA table 
summary.manova(mod, test = "Pillai") #4 MANOVA statistics are calculated and you can choose
summary.manova(mod, test = "Wilks")
summary.manova(mod, test = "Hotelling-Lawley")
summary.manova(mod, test = "Roy") #with a large sample size, all 4 statistics should converge well
#Rejected the null, there is a sig-dif in variances between morphologies

summary.aov(mod)# Individual variable significance test
# If there is a sig-dif, where is it? Chela, carp, abd, tl

##### oOo~Assumptions~oOo ####

### check for the two most important assumptions for MANOVA
#1. Multi-Normality assumption
#Graphic assessment of multi-normality: chi-square plot written by Everitt
#it is similar to a q-q plot in univariate ANOVA
chisplot <- function(x) {
    if (!is.matrix(x)) stop("x is not a matrix")

    ### determine dimensions
    n <- nrow(x)
    p <- ncol(x)
    #
    xbar <- apply(x, 2, mean)
    S <- var(x)
    S <- solve(S)
    index <- (1:n)/(n+1)
    #
    xcent <- t(t(x) - xbar)
    di <- apply(xcent, 1, function(x,S) x %*% S %*% x,S)
    #
    quant <- qchisq(index,p)
    plot(quant, sort(di), ylab = "Ordered distances",
         xlab = "Chi-square quantile", lwd=2,pch=1)
         return (di)

}
chisplot(resid(mod))  # the dataset must be in matrix format, original data

#run a normality test "mshapiro.test"(similar to "Shapiro.test") 
#but you need to install a package (mvnormtest) first
library(mvnormtest)
mshapiro.test(t(resid(mod))) #data has to be a matrix and has to be transposed

#check for individual variable's normality using Q-Q plot
#carp has the largest deviation from norm.
par(mfrow=c(3,3))
names<-colnames(y)
for (i in 1:ncol(y)){
  qqnorm(y[,i],main=paste("Q-Q plot of",names[i]))
  qqline(y[,i])
}

#####  oOo*~Refit~*oOo #####
library(mvnormtest)
#refit manova excluding outlier datapoints
ap.fit<-ap[-c(11,18,35,72),] #removed sites to improve multinormality-skewed (cooks.r)
x.fit<-as.factor(ap.fit[,12]) #independent variable (factor)
y.fit<-as.matrix(ap.fit[,c(2:7)])
mod.fit<-manova(y.fit~x.fit)

##test the normality assumption again
chisplot(resid(mod.fit))
mshapiro.test(t(resid(mod.fit)))  #small sample size issue!

### hypothesis test, alternative to ANOVA
library(vegan)
adonis(y ~ sex, data=env, perm=999)


#check normality after refit using Q-Q plot, again
par(mfrow=c(3,2))
names<-colnames(y.fit)
for (i in 1:ncol(y.fit)){
  qqnorm(y.fit[,i],main=paste("Q-Q plot of",names[i]))
  qqline(y.fit[,i])
}

par(mfrow=c(1,2))
chisplot(resid(mod)) #no transformation
chisplot(resid(mod.fit))#data points/skew removed

### check for the second most important assumptions for MANOVA
#2. homogeneity of covariance matrices
#boxplot to check if variance is equal among the groups
par(mfrow=c(2,4))
names<-colnames(y.fit)
for (i in 1:ncol(y)){
  boxplot(y.fit[,i]~x.fit,main=paste("Boxplot of",names[i]),xlab="Groups", ylab=paste(names[i]))
}

#log-transform
par(mfrow=c(2,4))
names<-colnames(y.fit)
for (i in 1:ncol(y)){
  boxplot(log(y.fit[,i]+1)~x.fit,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]))
}

par(mfrow=c(1,2))
##test homogeneity of covariance matrices 
##the test is developed by Anderson and is similar to Levene's univariate test for equal variance 
library(vegan) #Required functions for testing are in the 'vegan' Package
m<-betadisper(dist(y.fit),ap.fit$sex, type="centroid") #calculate average distance to its group centroid for each group
m.HSD<-TukeyHSD(m) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m.HSD) #graphic display of pair-wise comparisons of difference in dispersion between groups

#repeat the same test with log-transformed data
m.Log<-betadisper(dist(log(y.fit)+1),ap.fit$sex, type="centroid") #calculate average distance to its group centroid for each group
m.Log#Ave dist to centeroid should be small, if not small, assumptions not met
m.Log.HSD<-TukeyHSD(m.Log) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m.Log.HSD) #graphic display of pair-wise comparisons of difference in dispersion between groups


##### ~~~ o*Sex*o ~~~ ####

###Perform linear Discriminant Function Analysis (DFA)- Tries to seperate the groups
library(MASS) #"lda" function is in the MASS package
dfa<-lda(x.fit~y.fit) #run linear Discriminant Function Analysis
dfa #outputs, similar to PCA's output
plot(dfa, cex=1.5) #DFA plot (similar to a PCA plot); if more than two grouping catagories
tb<-table(Predicted=predict(dfa,ap.fit[,-c(1,8:12)],type='class')$class,Observed=ap.fit[,"sex"])
tb #show a confusion table (classification results using the DFA model)
((1-sum(diag(tb)/sum(tb)))*100) #mis-classification rate

library(klaR)
g<-greedy.wilks(scale(y),x)
g

