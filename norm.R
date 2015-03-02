#Bern Romey; 04Feb15
# Test normality on entire data frame

dta <- read.csv("apoch.csv")
am <- na.omit(dta)
rm(dta)
am <- am[c(6:12)]
am.lg <-log(am+1)

lshap <- lapply(am, shapiro.test) #shapiro test on normal data
lshap[[1]] ## look at the first column results
# You will need to extract the things you want from these objects, which all have the structure:
str(lshap[[1]])
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres) #transposed

lshap <- lapply(am.lg, shapiro.test) #shapiro test on log transformed data
lres <- sapply(lshap, `[`, c("statistic","p.value"))
lres <- t(lres)

norm <-table(lres)


#Theoretical norm, exp, log plots with qqplots
norm <- par(mfrow=c(3,2))
nm <- rnorm(53) #Normal distribution
hist(nm)
qqnorm(nm)
qqline(nm)

exp <- rexp(30, rate=1) #Exponential distribution
hist(exp)
qqnorm(exp)
qqline(exp)

lg <-log(exp)
hist(lg)
qqnorm(lg)
qqline(lg)
par(norm)



