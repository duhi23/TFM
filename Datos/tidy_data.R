# Universidad de Alcalá
# Master en Ciencias Actuariales y Financieras
# Diego Paul Huaraca Shagñay
# TFM: Riesgo de Longevidad en la Población Española

start <- Sys.time()
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
z <- as.matrix(Pxt_t[,17:66])
edad <- seq(0, 110, by=1)
tiempo <- seq(1965, 2014, by=1)
color <- colorRampPalette(c("red", "yellow", "green"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-25, phi=25, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='', ticktype="detailed", zlim=c(0,max(z)))
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
# suavizar <- function(data){
#       n_dim <- dim(data)
#       data <- as.matrix(data)
#       data_s <- mat.or.vec(n_dim[1],(n_dim[2]-1))
#       for(j in 2:n_dim[2]){
#             sp.spline <- smooth.Pspline(x = data[,1], y = data[,j], method = 4, norder = 3)
#             ajuste <- predict(sp.spline, x = data[,1])
#             data_s[,j-1] <- abs(round(ajuste[1:n_dim[1]],0))
#       }
#       colnames(data_s) <- colnames(data)[-1]
#       return(data_s)
# }


# Acotar data
acotar <- function(data, age){
      data <- as.matrix(data)[,-c(1)]
      ndim <- dim(data)
      ndata <- matrix(0, ncol=ndim[2], nrow=(age+1))
      ndata[1:age,] <- as.matrix(data)[1:age,]
      ndata[(age+1),] <- matrix(apply(as.matrix(data)[(age+1):ndim[1],], 2, sum),nrow=1)
      #ndata[(age+1),1] <- age
      colnames(ndata) <- colnames(data)
      return(ndata)
}


Ext_m <- acotar(Ext_m, 105)
Ext_h <- acotar(Ext_h, 105)
Ext_t <- acotar(Ext_t, 105)

Pxt_m <- acotar(Pxt_m, 105)
Pxt_h <- acotar(Pxt_h, 105)
Pxt_t <- acotar(Pxt_t, 105)

Dxt_m <- acotar(Dxt_m, 105)
Dxt_h <- acotar(Dxt_h, 105)
Dxt_t <- acotar(Dxt_t, 105)
rm(list=c("acotar"))

# Creación de tasas y probabilidades de muerte
ctasas <- function(dxt,ext,ini,fin){
      res <- dxt[,c(ini:fin)]/ext[,c(ini:fin)]
      return(res)
}
# data.frame(Asc=seq(1,65), Des=seq(65,1), Year=colnames(Dxt_h))
# 16=1965 - 55=2004
mxt_m <- ctasas(Dxt_m,Ext_m,16,55)
mxt_h <- ctasas(Dxt_h,Ext_h,16,55)
mxt_t <- ctasas(Dxt_t,Ext_t,16,55)
qxt_m <- ctasas(Dxt_m,Pxt_m,16,55)
qxt_h <- ctasas(Dxt_h,Pxt_h,16,55)
qxt_t <- ctasas(Dxt_t,Pxt_t,16,55)

# Gráfico tridimensional - Log tasas de mortalidad
z <- log(qxt_m)
edad <- seq(0, 105, by=1)
tiempo <- seq(1965, 2004, by=1)
color <- colorRampPalette(c("green", "yellow", "orange", "red"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='', ticktype="detailed", zlim=c(min(z),0))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))

# Creación objeto demogdata
data0_m <- StMoMoData(demogdata(data=qxt_m, pop=Pxt_m[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="female"), type="initial")
data0_h <- StMoMoData(demogdata(data=qxt_h, pop=Pxt_h[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="male"), type="initial")
data0_t <- StMoMoData(demogdata(data=qxt_t, pop=Pxt_t[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="total"), type="initial")

dataC_m <- StMoMoData(demogdata(data=mxt_m, pop=Ext_m[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="female"), type="central")
dataC_h <- StMoMoData(demogdata(data=mxt_h, pop=Ext_h[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="male"), type="central")
dataC_t <- StMoMoData(demogdata(data=mxt_t, pop=Ext_t[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                     type="mortality", label="Spain", name="total"), type="central")



################################################
#####   Ajuste modelo Lee Carter Poisson   #####
################################################

LC <- lc(link="log", const = "sum")
LCfit_m <- fit(LC, data = dataC_m)
LCfit_h <- fit(LC, data = dataC_h)
LCfit_t <- fit(LC, data = dataC_t)
rm(list = c("LC"))

# Gráfico parámetro ax
plot(LCfit_m$ax, type='b', col=2, pch=18, ylim=c(-10,0), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$ax, type='b', col=3, pch=20, ylim=c(-10,0), xlab='', ylab='')
legend(60, -5, legend=c("Hombres","Mujeres"), col=c("green","red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))


# Gráfico parametro bx
plot(LCfit_m$bx, type='b', col=2, pch=18, ylim=c(-0.045,0.05), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$bx, type='b', col=3, pch=20, ylim=c(-0.045,0.05), xlab='', ylab='')
abline(h=0,lty=3)
legend(45, 0.001, legend=c("Hombres","Mujeres"), col=c("green","red"), bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))


# Gráfico parametro kt
plot(seq(1965,2004),as.numeric(LCfit_m$kt), type='b', col=2, pch=18, ylim=c(-50,50), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),as.numeric(LCfit_h$kt), type='b', col=3, pch=20, ylim=c(-50,50), xlab='', ylab='')
abline(h=0,lty=3)
legend(1990, 50, legend=c("Hombres","Mujeres"), col=c("green","red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))


# Proyección parámetro kt
est_m <- forecast(LCfit_m, h=50, level=c(90,95))
est_h <- forecast(LCfit_h, h=50, level=c(90,95))
est_t <- forecast(LCfit_t, h=50, level=c(90,95))
plot(est_h, only.kt=TRUE, col='lightgreen')
plot(est_m, only.kt=TRUE, col='lightpink')
plot(est_t, only.kt=TRUE, col='lightblue')

# Comprobación modelos arima
tsm <- ts(as.numeric(LCfit_m$kt), start=c(1965, 1), end=c(2004, 1), frequency=1)
tsh <- ts(as.numeric(LCfit_h$kt), start=c(1965, 1), end=c(2004, 1), frequency=1)
summary(Arima(tsm, order=c(0,1,0), include.drift=TRUE))
summary(Arima(tsh, order=c(0,1,0), include.drift=TRUE))
rm(list = c("tsm", "tsh"))

# Obtención qxt estimados
for_qxt <- function(est){
      m <- exp(est$model$ax + est$model$bx %*% est$kt.f$mean)
      q <- 2*m/(2+m)
      colnames(q) <- colnames(m)
      return(m)
}

# Acota las probabilidades a 1

q1 <- function(mat){
      res <- apply(mat, 2, function(x){ifelse(x>1,1,x)})
      return(res)
}


# Gráfico tridimensional - Log tasas de mortalidad proyectadas
z <- log(q1(for_qxt(est_m)))
edad <- seq(0, 105, by=1)
tiempo <- seq(2005, 2054, by=1)
color <- colorRampPalette(c("green", "yellow", "orange", "red"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='', ticktype="detailed", zlim=c(min(z),0.1))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))



#######################################
######      Ajuste modelo CBD     #####
#######################################

fun_logit <- function(qxt){
      res<- log(1-qxt)-log(qxt)
      return(res)
}

# Gráfico logit qxt
q <- qxt_m
plot(fun_logit(q)[,1], col=2, type='p', ylim=c(-1.2,9.7), pch=18, ylab='')
par(new=TRUE)
plot(fun_logit(q)[,11], col=3, type='p', ylim=c(-1.2,9.7), pch=20, ylab='')
par(new=TRUE)
plot(fun_logit(q)[,21], col=4, type='p', ylim=c(-1.2,9.7), pch=18, ylab='')
par(new=TRUE)
plot(fun_logit(q)[,31], col=5, type='p', ylim=c(-1.2,9.7), pch=20, ylab='')
par(new=TRUE)
plot(fun_logit(q)[,40], col=6, type='p', ylim=c(-1.2,9.7), pch=18, ylab='logit(qxt)')
legend(6, 7, legend=c("1965", "1975", "1985", "1995", "2004"), col=c(2,3,4,5,6),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(18,20,18,20,18))
rm(list=c("q"))

# Modelo
CBD <- cbd("logit")
CBDfit_m <- fit(CBD, data = data0_m, ages.fit = 0:105, years=1965:2004)
CBDfit_h <- fit(CBD, data = data0_h, ages.fit = 0:105, years=1965:2004)
CBDfit_t <- fit(CBD, data = data0_t, ages.fit = 0:105, years=1965:2004)
rm(list=c("CBD"))

# Gráfico parámetro kt(1)
plot(seq(1965,2004),CBDfit_m$kt[1,], type='o', col=2, pch=18, ylim=c(-7,-4), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),CBDfit_h$kt[1,], type='o', col=3, pch=20, ylim=c(-7,-4), xlab='', ylab='')
legend(1980, -5.8, legend=c("Hombres", "Mujeres"), col=c("green", "red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))

# Gráfico parámetro kt(2)
plot(seq(1965,2004),CBDfit_m$kt[2,], type='o', col=2, pch=20, ylim=c(0.04,0.13), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),CBDfit_h$kt[2,], type='o', col=3, pch=20, ylim=c(0.04,0.13), xlab='', ylab='')
legend(1980, 0.077, legend=c("Hombres", "Mujeres"), col=c("green", "red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))

# Proyección parámetro kt
pro_m <- forecast(CBDfit_m, h=50, level=c(90,95))
pro_h <- forecast(CBDfit_h, h=50, level=c(90,95))
pro_t <- forecast(CBDfit_t, h=50, level=c(90,95))
plot(pro_m, only.kt=TRUE, col='lightpink')
plot(pro_h, only.kt=TRUE, col='lightgreen')
plot(pro_t, only.kt=TRUE, col='lightblue')


tsm <- ts(as.numeric(CBDfit_m$kt[2,]), start=c(1965, 1), end=c(2004, 1), frequency=1)
tsh <- ts(as.numeric(CBDfit_h$kt[2,]), start=c(1965, 1), end=c(2004, 1), frequency=1)
summary(Arima(tsm, order=c(0,1,0), include.drift=TRUE))
summary(Arima(tsh, order=c(0,1,0), include.drift=TRUE))
rm(list=c("tsm", "tsh"))


# Transforma vectores en matrices fila o columna
fm <- function(vector, rc){
      lon <- length(vector)
      if(rc==1){
            res <- matrix(unlist(vector), nrow=1, ncol=lon)
      } else {
            res <- matrix(unlist(vector), nrow=lon, ncol=1)
      }
      return(res)
} 


lg_qxt <- function(pro){
      # generar secuencia de 1
      x <- fm(rep(1,106),2)%*%fm(pro$kt.f$mean[1,],1) + fm(pro$ages-mean(pro$ages),2)%*% fm(pro$kt.f$mean[2,],1)
      q <- exp(x)/(1+exp(x))
      return(q)
}

lg_qxt(pro_h) - pro_h$rates

# Gráfico tridimensional - Log tasas de mortalidad
z <- log(pro_m$rates)
edad <- seq(0, 105, by=1)
tiempo <- seq(2005, 2054, by=1)
color <- colorRampPalette(c("green", "yellow", "orange", "red"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='', ticktype="detailed", zlim=c(min(z),0))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))


plot(log(pro_h$rates)[,1], type='b', col=2, ylim=c(-12,0))
par(new=TRUE)
plot(log(pro_h$rates)[,10], type='b', col=3, ylim=c(-12,0))
par(new=TRUE)
plot(log(pro_h$rates)[,20], type='b', col=4, ylim=c(-12,0))


########################################
#####     Ajuste modelo li-lee     #####
########################################

dataL_m <- demogdata(data=mxt_m, pop=Ext_m[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                                type="mortality", label="Spain", name="female")
dataL_h <- demogdata(data=mxt_h, pop=Ext_h[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                                type="mortality", label="Spain", name="male")
dataL_t <- demogdata(data=mxt_t, pop=Ext_t[,c(16:55)], ages=c(0:105), years=c(1965:2004), 
                                type="mortality", label="Spain", name="total")

LCT_m <- lca(dataL_m, max.age = 105)
LCT_h <- lca(dataL_h, max.age = 105)
LCT_t <- lca(dataL_t, max.age = 105)

# Gráfico parámetro ax - Hombres & Mujeres
plot(LCT_m$ax, type='b', col=2, pch=18, ylim=c(-10,0), xlab='', ylab='')
par(new=TRUE)
plot(LCT_h$ax, type='b', col=3, pch=20, ylim=c(-10,0), xlab='', ylab='')
legend(60, -5, legend=c("Hombres","Mujeres"), col=c("green","red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))

# Gráfico parámetro Bx - Total
plot(LCT_t$bx, type='o', col=4, pch=15, ylim=c(-0.026,0.036), xlab='', ylab='')
abline(h=0,lty=3)
legend(50, -0.003, legend=c("Total"), col=c("blue"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(15))

# Gráfico parámetro Kt - Total
plot(LCT_t$kt, type='o', col=4, pch=15, ylim=c(-38,32), xlab='', ylab='')
abline(h=0,lty=3)
legend(1980, -15, legend=c("Total"), col=c("blue"),   bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(15))


res_m <- log(mxt_m)-LCT_m$ax-fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1)
pres_m <- log(mxt_m)-LCT_m$ax
res_h <- log(mxt_h)-LCT_h$ax-fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1)
pres_h <- log(mxt_h)-LCT_h$ax


# Ratio ajuste factor común
Rc_m <- 1-(sum(res_m^2))/(sum(pres_m^2))
Rc_m
Rc_h <- 1-(sum(res_h^2))/(sum(pres_h^2))
Rc_h
rm(list = c("Rc_m", "Rc_h"))

# SVD residuos
svd_res <- function(residuos){
      res <- residuos #- apply(residuos, 1, mean)
      const <- sum(svd(res)$u[,1])
      b <- svd(res)$u[,1]/const
      k <- const*max(svd(res)$d)*svd(res)$v[,1]
      return(list(bx=b, kt=k))
}

 
# Gráfico parámetro bx - Hombres & Mujeres
plot(svd_res(res_m)$bx, type='o', col=2, pch=18, ylim=c(-0.015,0.06), xlab='', ylab='')
par(new=TRUE)
plot(svd_res(res_h)$bx, type='o', col=3, pch=20, ylim=c(-0.015,0.06), xlab='', ylab='')
abline(h=0,lty=3)
legend(50, 0.055, legend=c("Hombres","Mujeres"), col=c("green","red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))

# Gráfico parámetro kt - Hombres & Mujeres
plot(svd_res(res_m)$kt, type='o', col=2, pch=18, ylim=c(-15,15), xlab='', ylab='')
par(new=TRUE)
plot(svd_res(res_h)$kt, type='o', col=3, pch=20, ylim=c(-15,15), xlab='', ylab='')
abline(h=0,lty=3)
legend(15, -5, legend=c("Hombres","Mujeres"), col=c("green","red"),  bty = 'n', cex=0.9, 
       pt.cex=2, inset=0.2, pch=c(20,18))

# Ratio ajuste incluye residuos SVD
nres_m <- log(mxt_m) - LCT_m$ax - fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1) - fm(svd_res(res_m)$bx,2)%*%fm(svd_res(res_m)$kt,1)
nres_h <- log(mxt_h) - LCT_h$ax - fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1) - fm(svd_res(res_h)$bx,2)%*%fm(svd_res(res_h)$kt,1)
Rca_m <- 1-(sum(nres_m^2))/(sum(pres_m^2))
Rca_m
Rca_h <- 1-(sum(nres_h^2))/(sum(pres_h^2))
Rca_h
rm(list = c("Rca_m", "Rca_h", "nres_m", "nres_h"))


# Prueba de hipótesis kt: (H0: beta=1, entonces es un RW, caso contrario es un AR(1))
hreg <- function(res){
      alfa <- 0.05
      nkt <- length(svd_res(res)$kt)
      xm <- svd_res(res)$kt[-nkt]
      ym <- svd_res(res)$kt[-1]
      btest <- abs(summary(lm(ym ~ xm))$coef[2,1]-1)/summary(lm(ym ~ xm))$coef[2,2]
      if(pt(btest,summary(lm(ym ~ xm))$df[2], lower.tail=F)>alfa/2){
            res <- "Es un RW"
      } else {
            res <- "Es un AR(1)"
      }
      btest <- abs(summary(lm(ym ~ xm+0))$coef[1,1]-1)/summary(lm(ym ~ xm+0))$coef[1,2]
      if(pt(btest,summary(lm(ym ~ xm+0))$df[2], lower.tail=F)>alfa/2){
            res1 <- "Es un RW"
      } else {
            res1 <- "Es un AR(1)"
      }
      return(list(res, res1))
}

hreg(res_m)
hreg(res_h)



# Proyección parámetro Kt y kt
pry_t <- forecast(auto.arima(LCT_t$kt), h=50, level=c(90,95))
plot(pry_t,col='blue')

tpm <- ts(svd_res(res_m)$kt, start=c(1965, 1), end=c(2004, 1), frequency=1)
tph <- ts(svd_res(res_h)$kt, start=c(1965, 1), end=c(2004, 1), frequency=1)

pry_m <- forecast(auto.arima(tpm), h=50, level=c(90,95))
plot(pry_m, col=1, fcol=1, shadecols="lightpink")
pry_h<- forecast(auto.arima(tph), h=50, level=c(90,95))
plot(pry_h, col=1, fcol=1, shadecols="lightgreen")


li_qxt <- function(LCa, LCTb, res, pryt, pry){
      estm <- exp(LCa$ax + fm(LCTb$bx,2)%*%fm(as.numeric(pryt$mean),1) + fm(svd_res(res_m)$bx,2)%*%fm(as.numeric(pry$mean),1))
      q <- 2*estm/(2+estm)
      colnames(q) <- seq(2005, 2054)
      return(q)
}


# Gráfico tridimensional - Log probabilidades de muerte
z <- log(q1(li_qxt(LCT_m, LCT_t, res_m, pry_t, pry_m)))
edad <- seq(0, 105, by=1)
tiempo <- seq(2005, 2054, by=1)
color <- colorRampPalette(c("green", "yellow", "orange", "red"))(50)
zfacet <- z[-1, -1] + z[-1, -length(tiempo)] + z[-length(edad), -1] + z[-length(edad), -length(tiempo)]
facetcol <- cut(zfacet, 50)
persp(edad, tiempo, z, theta=-30, phi=30, expand=0.75, xlab='Edad', col=color[facetcol],
      ylab='Periodo', zlab='', ticktype="detailed", zlim=c(min(z),0))
rm(list=c("z", "edad", "tiempo", "color", "zfacet", "facetcol"))


########################################
#####     Análisis de residuos     #####
########################################

# Error cuadrático medio #

ecm <- function(datar, datae){
      error <- apply(datar - datae, 2, function(x){sqrt(sum(x^2/length(x)))})
      return(error)
}

# Error porcentual absoluto medio #

epam <- function(datar, datae){
      dife <- abs(datar-datae)/datar
      error <- apply(dife, 2, function(x){sum(x)/length(x)})
      return(error)
}

# Almacenamiento de resultados
LCPm <- q1(for_qxt(est_m))
LCPh <- q1(for_qxt(est_h))
CBDm <- q1(lg_qxt(pro_m))
CBDh <- q1(lg_qxt(pro_h))
LLm <- q1(li_qxt(LCT_m, LCT_t, res_m, pry_t, pry_m))
LLh <- q1(li_qxt(LCT_h, LCT_t, res_h, pry_t, pry_h))
save(list = c("LCPm", "LCPh", "LLm", "LLh", "CBDm", "CBDh"), file = "Estimaciones.RData", envir = .GlobalEnv)

# LCP
em1 <- ecm(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], for_qxt(est_m)[15:106,c(1:10)])
eh1 <- ecm(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], for_qxt(est_h)[15:106,c(1:10)])

em4 <- epam(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], for_qxt(est_m)[15:106,c(1:10)])
eh4 <- epam(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], for_qxt(est_h)[15:106,c(1:10)])

# LL
em2 <- ecm(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], li_qxt(LCT_m, LCT_t, res_m, pry_t, pry_m)[15:106,c(1:10)])
eh2 <- ecm(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], li_qxt(LCT_h, LCT_t, res_h, pry_t, pry_h)[15:106,c(1:10)])

em5 <- epam(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], li_qxt(LCT_m, LCT_t, res_m, pry_t, pry_m)[15:106,c(1:10)])
eh5 <- epam(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], li_qxt(LCT_h, LCT_t, res_h, pry_t, pry_h)[15:106,c(1:10)])

# CBD
em3 <- ecm(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], lg_qxt(pro_m)[15:106,c(1:10)])
eh3 <- ecm(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], lg_qxt(pro_h)[15:106,c(1:10)])

em6 <- epam(ctasas(Dxt_m,Pxt_m,56,65)[15:106,], lg_qxt(pro_m)[15:106,c(1:10)])
eh6 <- epam(ctasas(Dxt_h,Pxt_h,56,65)[15:106,], lg_qxt(pro_h)[15:106,c(1:10)])

xtable(data.frame(LCP=eh1, LL=eh2, CBD=eh3, LCP1=eh4, LL1=eh5, CBD1=eh6), digits=6)

xtable(data.frame(LCP=em1, LL=em2, CBD=em3, LCP1=em4, LL1=em5, CBD1=em6), digits=6)


# Tiempo ejecución & Limpieza Workspace
(Sys.time() - start)
rm(list = ls())

