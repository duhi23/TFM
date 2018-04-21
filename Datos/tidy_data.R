# Universidad de Alcalá
# Master en Ciencias Actuariales y Financieras
# Diego Paul Huaraca Shagñay
# TFM: Riesgo de Longevidad en la Población Española

# Datos - Human Mortality Database
options(tz="Europe/Madrid")
libs <- c('readxl', 'sme', 'tidyverse', 'forecast', 'foreign', 'MASS', 'demography', 'lifecontingencies',
          'tseries','pspline','psych','xtable', 'dplyr', 'StMoMo')
lapply(libs, require, character.only = T)
rm(list=c("libs"))

# Poblacion al 1 de Enero (Expuesto inicial al riesgo)
Pxt <- read_excel("./Datos/Data_Spain.xlsx", sheet=1)
Pxt_m <- Pxt %>% dplyr::select(1:3) %>% spread(key="Year", value="Female")
Pxt_h <- Pxt %>% dplyr::select(1,2,4) %>% spread(key="Year", value="Male")
Pxt_t <- Pxt %>% dplyr::select(1,2,5) %>% spread(key="Year", value="Total")
rm(list=c("Pxt"))

# Gráfico tridimensional - data original
z <- as.matrix(Pxt_t[,2:66])
edad <- seq(0, 110, by=1)
tiempo <- seq(1950, 2014, by=1)
color <- colorRampPalette(c("red", "yellow", "green"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='Población', ticktype="detailed", zlim=c(0,max(z)))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))

# Número de defunciones correspondiente al año t
Dxt <- read_excel("./Datos/Data_Spain.xlsx", sheet=2)
Dxt_m <- Dxt %>% dplyr::select(1:3) %>% spread(key="Year", value="Female")
Dxt_h <- Dxt %>% dplyr::select(1,2,4) %>% spread(key="Year", value="Male")
Dxt_t <- Dxt %>% dplyr::select(1,2,5) %>% spread(key="Year", value="Total")
rm(list=c("Dxt"))

# Poblacion Expuesta central al riesgo 
Ext <- read_excel("./Datos/Data_Spain.xlsx", sheet=3)
Ext_m <- Ext %>% dplyr::select(1:3) %>% spread(key="Year", value="Female")
Ext_h <- Ext %>% dplyr::select(1,2,4) %>% spread(key="Year", value="Male")
Ext_t <- Ext %>% dplyr::select(1,2,5) %>% spread(key="Year", value="Total")
rm(list=c("Ext"))

# Función suavizado superficies - P-Splines
suavizar <- function(data){
      n_dim <- dim(data)
      data <- as.matrix(data)
      data_s <- mat.or.vec(n_dim[1],(n_dim[2]-1))
      for(j in 2:n_dim[2]){
            sp.spline <- smooth.Pspline(x = data[,1], y = data[,j], method = 4, norder = 3)
            ajuste <- predict(sp.spline, x = data[,1])
            data_s[,j-1] <- abs(round(ajuste[1:n_dim[1]],0))
      }
      colnames(data_s) <- colnames(data)[-1]
      return(data_s)
}

# Acotar data
acotar <- function(data, age){
      ndim <- dim(data)
      ndata <- matrix(0, ncol=ndim[2], nrow=(age+1))
      ndata[1:age,] <- as.matrix(data)[1:age,]
      ndata[(age+1),] <- matrix(apply(as.matrix(data)[(age+1):ndim[1],], 2, sum),nrow=1)
      #ndata[(age+1),1] <- age
      colnames(ndata) <- colnames(data)
      return(ndata)
}

# Gráfico tridimensional - data suavizada & acotada
z <- as.matrix(acotar(suavizar(Pxt_t), 105))
edad <- seq(0, 105, by=1)
tiempo <- seq(1950, 2014, by=1)
color <- colorRampPalette(c("red", "yellow", "green"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='Población', ticktype="detailed", zlim=c(0,max(z)))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))


Ext_m <- acotar(suavizar(Ext_m), 105)
Ext_h <- acotar(suavizar(Ext_h), 105)
Ext_t <- acotar(suavizar(Ext_t), 105)

Pxt_m <- acotar(suavizar(Pxt_m), 105)
Pxt_h <- acotar(suavizar(Pxt_h), 105)
Pxt_t <- acotar(suavizar(Pxt_t), 105)

Dxt_m <- acotar(suavizar(Dxt_m), 105)
Dxt_h <- acotar(suavizar(Dxt_h), 105)
Dxt_t <- acotar(suavizar(Dxt_t), 105)


# Creación de tasas y probabilidades de muerte
mxt_m <- Dxt_m[,-c(57:66)]/Ext_m[,-c(57:66)]
mxt_h <- Dxt_h[,-c(57:66)]/Ext_h[,-c(57:66)]
mxt_t <- Dxt_t[,-c(57:66)]/Ext_t[,-c(57:66)]
qxt_m <- Dxt_m[,-c(57:66)]/Pxt_m[,-c(57:66)]
qxt_h <- Dxt_h[,-c(57:66)]/Pxt_h[,-c(57:66)]
qxt_t <- Dxt_t[,-c(57:66)]/Pxt_t[,-c(57:66)]


# Creación objeto demogdata
data0_m <- StMoMoData(demogdata(data=qxt_m, pop=Pxt_m[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="female"), type="initial")
data0_h <- StMoMoData(demogdata(data=qxt_h, pop=Pxt_h[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="male"), type="initial")
data0_t <- StMoMoData(demogdata(data=qxt_t, pop=Pxt_t[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="total"), type="initial")

dataC_m <- StMoMoData(demogdata(data=mxt_m, pop=Ext_m[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="female"), type="central")
dataC_h <- StMoMoData(demogdata(data=mxt_h, pop=Ext_h[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="male"), type="central")
dataC_t <- StMoMoData(demogdata(data=mxt_t, pop=Ext_t[,-c(57:66)], ages=c(0:105), years=c(1950:2005), 
                     type="mortality", label="Spain", name="total"), type="central")

rm(list=c("mxt_m", "mxt_h", "mxt_t", "qxt_m", "qxt_h", "qxt_t", "Pxt_m", "Pxt_h", "Pxt_t",
          "Ext_m", "Ext_h", "Ext_t", "Dxt_m", "Dxt_h", "Dxt_t", "suavizar"))

#############################################
####   Ajuste modelo Lee Carter Poisson   ###
#############################################

LC <- lc(link="log", const = "sum")
LCfit_m <- fit(LC, data = dataC_m, ages.fit = 0:105, years=1950:2005)
LCfit_h <- fit(LC, data = dataC_h, ages.fit = 0:105, years=1950:2005)
LCfit_t <- fit(LC, data = dataC_t, ages.fit = 0:105, years=1950:2005)

# Gráfico parametro ax
plot(LCfit_m$ax, type='o', col=2, pch=20, ylim=c(-10,0), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$ax, type='o', col=3, pch=20, ylim=c(-10,0), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_t$ax, type='o', col=4, pch=18, ylim=c(-10,0), xlab='Edad', ylab='ax')
legend(65, -6, legend=c("Mujeres", "Hombres", "Total"), col=c("red", "green", "blue"), cex=0.5,
       box.lty=0, text.font=10, lwd=2)

# Gráfico parametro bx
plot(LCfit_m$bx, type='o', col=2, pch=20, ylim=c(-0.025,0.045), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$bx, type='o', col=3, pch=20, ylim=c(-0.025,0.045), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_t$bx, type='o', col=4, pch=18, ylim=c(-0.025,0.045), xlab='Edad', ylab='bx')
abline(h=0,lty=3)
legend(45, -0.002, legend=c("Mujeres", "Hombres", "Total"), col=c("red", "green", "blue"), cex=0.5,
       box.lty=0, text.font=10, lwd=2)

# Gráfico parametro kt
plot(seq(1950,2005),as.numeric(LCfit_m$kt), type='o', col=2, pch=20, ylim=c(-75,80), xlab='', ylab='')
par(new=TRUE)
plot(seq(1950,2005),as.numeric(LCfit_h$kt), type='o', col=3, pch=20, ylim=c(-75,80), xlab='', ylab='')
par(new=TRUE)
plot(seq(1950,2005),as.numeric(LCfit_t$kt), type='o', col=4, pch=18, ylim=c(-75,80), xlab='Periodo', ylab='kt')
abline(h=0,lty=3)
legend(1980, 70, legend=c("Hombres", "Mujeres", "Total"), col=c("red", "green", "blue"), cex=0.5,
       box.lty=0, text.font=10, lwd=2)

# Proyección parámetro kt
est_m <- forecast(LCfit_m, h=40, level=c(90,95))
est_h <- forecast(LCfit_h, h=40, level=c(90,95))
est_t <- forecast(LCfit_t, h=40, level=c(90,95))
plot(est_m, only.kt=TRUE, col='lightgreen')
plot(est_h, only.kt=TRUE, col='lightpink')
plot(est_t, only.kt=TRUE, col='lightblue')


#######################################
######      Ajuste modelo CBD     #####
#######################################

CBD <- cbd()
CBDfit_m <- fit(CBD, data = data0_m, ages.fit = 60:105, years=1950:2005)
CBDfit_h <- fit(CBD, data = data0_h, ages.fit = 60:105, years=1950:2005)
CBDfit_t <- fit(CBD, data = data0_t, ages.fit = 60:105, years=1950:2005)

# Gráfico parametro kt(1)
plot(CBDfit_m$kt[1,], type='o', col=2, pch=20, ylim=c(-3,-1.5), xlab='', ylab='')
par(new=TRUE)
plot(CBDfit_h$kt[1,], type='o', col=3, pch=20, ylim=c(-3,-1.5), xlab='', ylab='')
par(new=TRUE)
plot(CBDfit_t$kt[1,], type='o', col=4, pch=18, ylim=c(-3,-1.5), xlab='Edad', ylab='kt(1)')
legend(43, -1.5, legend=c("Mujeres", "Hombres", "Total"), col=c("red", "green", "blue"), cex=0.5,
       box.lty=0, text.font=10, lwd=2)

# Gráfico parametro kt(2)
plot(CBDfit_m$kt[2,], type='o', col=2, pch=20, ylim=c(0.07,0.15), xlab='', ylab='')
par(new=TRUE)
plot(CBDfit_h$kt[2,], type='o', col=3, pch=20, ylim=c(0.07,0.15), xlab='', ylab='')
par(new=TRUE)
plot(CBDfit_t$kt[2,], type='o', col=4, pch=18, ylim=c(0.07,0.15), xlab='Edad', ylab='kt(2)')
legend(35, 0.10, legend=c("Mujeres", "Hombres", "Total"), col=c("red", "green", "blue"), cex=0.5,
       box.lty=0, text.font=10, lwd=2)

# Proyección parámetro kt
pro_m <- forecast(CBDfit_m, h=40, level=c(90,95))
pro_h <- forecast(CBDfit_h, h=40, level=c(90,95))
pro_t <- forecast(CBDfit_t, h=40, level=c(90,95))
plot(pro_m, only.kt=TRUE, col='lightgreen')
plot(pro_h, only.kt=TRUE, col='lightpink')
plot(pro_t, only.kt=TRUE, col='lightblue')


########################################
#####     Ajuste modelo li-lee     #####
########################################

