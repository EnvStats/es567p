names(mod)# what type of informatin is produced from NMDS.
mod$points

plot(mod$points[,1],mod$points[,2]) 
head(apo)
cor(mod$points[,1],ap[,c(2:8)]) #Is there a correlation with NMDS axis one?
cor(mod$points[,2],ap[,c(3:8)]) # Is there a correlation with NMDS axis two?

plot(mod$points[,1],ap[,3]) #Plots original variable  Chela against NMDS.
symbols(mod$points[,1],mod$points[,2],circles=ap[,3],inches=0.2)
plot(mod$points[,1],mod$points[,3],type='n')

plot(mod$points[,1],mod$points[,2],type='n') #Type=n leaves the plot empty for symbols
symbols(mod$points[,1],mod$points[,2],circles=ap[,3],inches=0.1)




