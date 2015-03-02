#Cooks Distance

apo <-read.csv("apoch.csv")
head(apo)
ap<-apo[,c(2:8)]# only morphology variables

lshap <- lapply(ap, shapiro.test) #shapiro test for normality; used to select best norm for lm 
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres) #transposed

##### --linear model and cooks distance analysis-- ####

par(mfrow = c(3,2))
lm.ap <- lm(tl ~ chela + carp + cheli + tib + fem + abd, data = ap) 
plot(lm.ap, which = 1:6)
summary(lm.ap)
anova(lm.ap)

inflm.ap <- influence.measures(lm.ap)
which(apply(inflm.ap$is.inf, 1, any)) # which observations are influential

summary(inflm.ap) # only the influential ones
inflm.ap          # all, with influance column


##### ~reduced~ ####
ap.r<-ap[-c(11,18,35,72),] #removed

lshap <- lapply(ap.r, shapiro.test) #shapiro test for normality; used to select best norm for lm 
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres) #transposed

##### --linear model and cooks distance analysis-- ####

par(mfrow = c(3,2))
lm.ap <- lm(tl ~ chela + carp + cheli + tib + fem + abd, data = ap.r) 
plot(lm.ap, which = 1:6)


