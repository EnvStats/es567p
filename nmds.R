# Bern - NMDS
# ~~ Nonmetric Multidimensional Scaling with Stable Solution from Random Starts, 
# Axis Scaling and Species Scores. Alternative to MANOVA ~~

# Adequacy is dependant on: 
# Stress value < 10, check Shepard plot for outliers, Check distortion: select top 10% of values 
# in Bray-Curtis matrix and compare to the MDS ordination distance

# What to report:
# Distance measure used, Algorithm/Software used, Number of runs, Number of dimensions in the final solution
# Interpretation of each axis

apo <-read.csv("apoch.csv")
ap <- apo[-c(32:33),] #removed becuase only two land form for R
head(ap)
Y<-as.matrix(ap[,c(2:8)]) #multivariate response variable matrix (morphology)
env <-ap[,c(9:12)] #Environemtal group

library(vegan)
library(MASS)
mod<-metaMDS(Y,k=2,trace=T,autotransform =F) #MDS with 2 dimension
br<-vegdist(Y)  #calculate Bray-Curtis distance among all sites
stressplot(mod,br)  #generate a Sharperd diagram by plotting BC between sites against Eucledian distance between sites in NMDS 
mod$stress #check for NMDS' stress value
round(mod$species,2)

par(mfrow=c(1,2))
plot(mod, type="t") # land Form
### Hulls show treatment
with(env, ordihull(mod, group=form, show="T"))
with(env, ordihull(mod, group=form, show="A", col="green"))
with(env, ordihull(mod, group=form, show="C", col="red")) 
### Spider shows fields
with(env, ordispider(mod, group=form, lty=4, col="darkgreen"))

plot(mod, type="t") # Sex
### Hulls show treatment
with(env, ordihull(mod, group=sex, show="m"))
with(env, ordihull(mod, group=sex, show="f", col="green")) 
### Spider shows fields
with(env, ordispider(mod, group=sex, lty=5, col="blue"))

### hypothesis test (with strata)
adonis(Y ~ form*sex*eco4, data=env, perm=999)

#-------------------
##### example of strata, nested (block) design.
plot(mod, type="t") # Eco4
### Hulls show treatment, within field
with(env, ordihull(mod, group=sex, show="m"))
with(env, ordihull(mod, group=sex, show="f", col="green")) #group = strata
### Spider shows fields
with(env, ordispider(mod, group=eco4, lty=5, col="blue"))
#-------------------



#example from R
data(dune)
data(dune.env)
adonis(dune ~ Management*A1, data=dune.env, permutations=99)

### Example of use with strata, for nested (e.g., block) designs.

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment,within field
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))
### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)

