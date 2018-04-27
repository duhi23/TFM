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


# Creación de tasas y probabilidades de muerte
ctasas <- function(dxt,ext,ini,fin){
      res <- dxt[,c(ini:fin)]/ext[,c(ini:fin)]
      return(res)
}
# data.frame(Asc=seq(1,65), Des=seq(65,1), Year=colnames(Dxt_h))

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

#rm(list=c("mxt_m", "mxt_h", "mxt_t", "qxt_m", "qxt_h", "qxt_t", "Pxt_m", "Pxt_h", "Pxt_t",
#          "Ext_m", "Ext_h", "Ext_t", "Dxt_m", "Dxt_h", "Dxt_t", "suavizar"))


################################################
#####   Ajuste modelo Lee Carter Poisson   #####
################################################

LC <- lc(link="log", const = "sum")
LCfit_m <- fit(LC, data = dataC_m)
LCfit_h <- fit(LC, data = dataC_h)
LCfit_t <- fit(LC, data = dataC_t)

# Gráfico parámetro ax
plot(LCfit_m$ax, type='o', col=2, pch=20, ylim=c(-10,0), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$ax, type='o', col=3, pch=20, ylim=c(-10,0), xlab='', ylab='')
#par(new=TRUE)
#plot(LCfit_t$ax, type='o', col=4, pch=18, ylim=c(-10,0), xlab='Edad', ylab='ax')
#op <- par(cex = 2.8)
legend(60, -5, legend=c("Mujeres", "Hombres"), col=c("red", "green"), cex=0.9, lwd=3, bty = 'n')



# Gráfico parametro bx
plot(LCfit_m$bx, type='o', col=2, pch=20, ylim=c(-0.045,0.05), xlab='', ylab='')
par(new=TRUE)
plot(LCfit_h$bx, type='o', col=3, pch=20, ylim=c(-0.045,0.05), xlab='', ylab='')
#par(new=TRUE)
#plot(LCfit_t$bx, type='o', col=4, pch=18, ylim=c(-0.035,0.05), xlab='Edad', ylab='bx')
abline(h=0,lty=3)
legend(45, 0.005, legend=c("Mujeres", "Hombres"), col=c("red", "green"), cex=0.9, lwd=3, bty = 'n')

# Gráfico parametro kt
plot(seq(1965,2004),as.numeric(LCfit_m$kt), type='o', col=2, pch=20, ylim=c(-55,55), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),as.numeric(LCfit_h$kt), type='o', col=3, pch=20, ylim=c(-55,55), xlab='', ylab='')
#par(new=TRUE)
#plot(seq(1965,2004),as.numeric(LCfit_t$kt), type='o', col=4, pch=18, ylim=c(-55,55), xlab='Periodo', ylab='kt')
abline(h=0,lty=3)
legend(1985, 60, legend=c("Hombres", "Mujeres"), col=c("red", "green"), cex=0.9, lwd=3, bty = 'n')

# Proyección parámetro kt
est_m <- forecast(LCfit_m, h=50, level=c(90,95))
est_h <- forecast(LCfit_h, h=50, level=c(90,95))
est_t <- forecast(LCfit_t, h=50, level=c(90,95))
plot(est_m, only.kt=TRUE, col='lightgreen')
plot(est_h, only.kt=TRUE, col='lightpink')
plot(est_t, only.kt=TRUE, col='lightblue')


# Obtención qxt estimados
for_qxt <- function(est){
      m <- exp(est$model$ax + est$model$bx %*% est$kt.f$mean)
      q <- 2*m/(2+m)
      colnames(q) <- colnames(m)
      return(m)
}

for_qxt(est_m)

# Error cuadrático medio
error_m <- apply(ctasas(Dxt_m,Pxt_m,56,65) - for_qxt(est_m)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
error_m
error_h <- apply(ctasas(Dxt_h,Pxt_h,56,65) - for_qxt(est_h)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
error_h
error_t <- apply(ctasas(Dxt_t,Pxt_t,56,65) - for_qxt(est_t)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
error_t


plot(log(ctasas(Dxt_m,Pxt_m,56,65))[,1], col=2, type='o', ylim=c(-10,0))
par(new=TRUE)
plot(log(for_qxt(est_m))[,1], col=3, type='o', ylim=c(-10,0))


plot(log(ctasas(Dxt_h,Pxt_h,56,65))[,1], col=2, type='o', ylim=c(-10,0))
par(new=TRUE)
plot(log(for_qxt(est_h))[,1], col=3, type='o', ylim=c(-10,0))


#######################################
######      Ajuste modelo CBD     #####
#######################################

fun_logit <- function(qxt){
      res<- log(1-qxt)-log(qxt)
      return(res)
}

plot(fun_logit(qxt_h)[,20])

CBD <- cbd("logit")
CBDfit_m <- fit(CBD, data = data0_m, ages.fit = 0:105, years=1965:2004)
CBDfit_h <- fit(CBD, data = data0_h, ages.fit = 0:105, years=1965:2004)
CBDfit_t <- fit(CBD, data = data0_t, ages.fit = 0:105, years=1965:2004)

# Gráfico parámetro kt(1)
plot(seq(1965,2004),CBDfit_m$kt[1,], type='o', col=2, pch=20, ylim=c(-7,-4), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),CBDfit_h$kt[1,], type='o', col=3, pch=20, ylim=c(-7,-4), xlab='', ylab='')
#par(new=TRUE)
#plot(seq(1965,2004),CBDfit_t$kt[1,], type='o', col=4, pch=18, ylim=c(-3,-1.5), xlab='Edad', ylab='kt(1)')
legend(1980, -5.8, legend=c("Hombres", "Mujeres"), col=c("red", "green"), cex=0.9, lwd=3, bty = 'n')

# Gráfico parámetro kt(2)
plot(seq(1965,2004),CBDfit_m$kt[2,], type='o', col=2, pch=20, ylim=c(0.04,0.13), xlab='', ylab='')
par(new=TRUE)
plot(seq(1965,2004),CBDfit_h$kt[2,], type='o', col=3, pch=20, ylim=c(0.04,0.13), xlab='', ylab='')
#par(new=TRUE)
#plot(seq(1965,2004),CBDfit_t$kt[2,], type='o', col=4, pch=18, ylim=c(0.07,0.15), xlab='Edad', ylab='kt(2)')
legend(1980, 0.07, legend=c("Hombres", "Mujeres"), col=c("red", "green"), cex=0.9, lwd=3, bty = 'n')

# Proyección parámetro kt
pro_m <- forecast(CBDfit_m, h=50, level=c(90,95))
pro_h <- forecast(CBDfit_h, h=50, level=c(90,95))
pro_t <- forecast(CBDfit_t, h=50, level=c(90,95))
plot(pro_m, only.kt=TRUE, col='lightgreen')
plot(pro_h, only.kt=TRUE, col='lightpink')
plot(pro_t, only.kt=TRUE, col='lightblue')

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
      x <- fm(rep(1,46),2)%*%fm(pro$kt.f$mean[1,],1) + fm(pro$ages-mean(pro$ages),2)%*% fm(pro$kt.f$mean[2,],1)
      q <- exp(x)/(1+exp(x))
      return(q)
}

lg_qxt(pro_m)

# Error cuadrático medio
err_m <- apply(ctasas(Dxt_m,Pxt_m,56,65)[61:106,] - lg_qxt(pro_m)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
err_m
err_h <- apply(ctasas(Dxt_h,Pxt_h,56,65)[61:106,] - lg_qxt(pro_h)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
err_h
err_t <- apply(ctasas(Dxt_t,Pxt_t,56,65)[61:106,] - lg_qxt(pro_t)[,c(1:10)], 2, function(x){sqrt(sum(x^2/length(x)))})
err_t



########################################
#####     Ajuste modelo li-lee     #####
########################################


LCT <- lc(link="logit", const = "sum")
LCT_m <- fit(LCT, data = data0_m, ages.fit = 0:105, years=1950:2005)
LCT_h <- fit(LCT, data = data0_h, ages.fit = 0:105, years=1950:2005)
LCT_t <- fit(LCT, data = data0_t, ages.fit = 0:105, years=1950:2005)

res_m <- log(qxt_m)-LCT_m$ax-fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1)
pres_m <- log(qxt_m)-LCT_m$ax
res_h <- log(qxt_h)-LCT_h$ax-fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1)
pres_h <- log(qxt_h)-LCT_h$ax

# Ratio ajuste factor común
Rc_m <- 1-(sum(res_m^2))/(sum(pres_m^2))
Rc_h <- 1-(sum(res_h^2))/(sum(pres_h^2))

# SVD residuos
svd_res <- function(residuos){
      res <- residuos - apply(residuos, 1, mean)
      const <- sum(svd(res)$u[,1])
      b <- svd(res)$u[,1]/const
      k <- const*svd(res)$d[1]*svd(res)$v[,1]
      return(list(bx=b, kt=k))
}
 
svd_res(res_m)
 
sum(svd_res(res_m)$bx)
sum(svd_res(res_m)$kt)

# Ratio ajuste incluye residuos SVD
nres_m <- log(qxt_m) - LCT_m$ax - fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1) - fm(svd_res(res_m)$bx,2)%*%fm(svd_res(res_m)$kt,1)
nres_h <- log(qxt_h) - LCT_h$ax - fm(LCT_t$bx,2)%*%fm(LCT_t$kt,1) - fm(svd_res(res_h)$bx,2)%*%fm(svd_res(res_h)$kt,1)
Rca_m <- 1-(sum(nres_m^2))/(sum(pres_m^2))
Rca_h <- 1-(sum(nres_h^2))/(sum(pres_h^2))



## Obtención de la ratio Rac

res1M=resM-bxM%*%t(ktM)
num=sum(res1M^2)
den=sum((log(dM)-axM)^2)
RacM=1-num/den
RacM

res1F=resF-bxF%*%t(ktF)
num=sum(res1F^2)
den=sum((log(dF)-axF)^2)
RacF=1-num/den
RacF

## Comprobación de la existencia de un AR(1) en los kti (H0: beta=1, entonces es un RW)

alfa=0.05
N=length(ktM)
XM=ktM[-N]
YM=ktM[-1]

MX=lm(YM~XM)
summary(MX)

i=length(MX$coef)
test=abs((MX$coef[i]-1)/summary(MX)$coef[i,2])
p.val=pt(test,summary(MX)$df[2],lower.tail=F)
if(p.val>alfa/2) "Es un RW" else "Es un AR(1)"
MP[1,1]=summary(MX)$coef[1,1]
MP[1,2]=summary(MX)$coef[1,4]
MP[1,3]=summary(MX)$coef[2,1]
MP[1,4]=summary(MX)$coef[2,4]
MP[1,5]=if(p.val>alfa/2) print("RW",quote=F) else print("AR(1)",quote=F)

MX2=lm(YM~XM+0)
summary(MX2)

i=length(MX2$coef)
test=abs((MX2$coef[i]-1)/summary(MX2)$coef[i,2])
p.val=pt(test,summary(MX2)$df[2],lower.tail=F)
if(p.val>alfa/2) "Es un RW" else "Es un AR(1)"
MP2[1,3]=summary(MX2)$coef[1,1]
MP2[1,4]=summary(MX2)$coef[1,4]
MP2[1,5]=if(p.val>alfa/2) "RW" else "AR(1)"

XF=ktF[-N]
YF=ktF[-1]

MX=lm(YF~XF)
summary(MX)

i=length(MX$coef)
test=abs((MX$coef[i]-1)/summary(MX)$coef[i,2])
p.val=pt(test,summary(MX)$df[2],lower.tail=F)
if(p.val>alfa/2) "Es un RW" else "Es un AR(1)"
MP[2,1]=summary(MX)$coef[1,1]
MP[2,2]=summary(MX)$coef[1,4]
MP[2,3]=summary(MX)$coef[2,1]
MP[2,4]=summary(MX)$coef[2,4]
MP[2,5]=if(p.val>alfa/2) "RW" else "AR(1)"

MX2=lm(YF~XF+0)
summary(MX2)

i=length(MX2$coef)
test=abs((MX2$coef[i]-1)/summary(MX2)$coef[i,2])
p.val=pt(test,summary(MX2)$df[2],lower.tail=F)
if(p.val>alfa/2) "Es un RW" else "Es un AR(1)"
MP2[2,3]=summary(MX2)$coef[1,1]
MP2[2,4]=summary(MX2)$coef[1,4]
MP2[2,5]=if(p.val>alfa/2) "RW" else "AR(1)"


#### Almacenamiento de los kt (global y por sexos)

kappa=cbind(K,ktM,ktF)
write.table(kappa,"kappaT.txt",row.names=T)


final=Sys.time()

(duración=final-inicio)

