# Universidad de Alcalá
# Master en Ciencias Actuariales y Financieras
# Diego Paul Huaraca Shagñay
# TFM: Riesgo de Longevidad en la Población Española

start <- Sys.time()
options(tz="Europe/Madrid")
libs <- c('readxl', 'sme', 'tidyverse', 'forecast', 'foreign', 'MASS', 'demography', 'lifecontingencies',
          'tseries','pspline','psych','xtable', 'dplyr', 'StMoMo')
lapply(libs, require, character.only = T)
rm(list=c("libs"))

# Acota las probabilidades a 1
q1 <- function(mat){
      res <- apply(mat, 2, function(x){ifelse(x>1,1,x)})
      res[nrow(res),] <- rep(1,ncol(res))
      return(res)
}

# Obtención de probabilidades proyectadas
perm <- read_excel("Datos/Tablas_PERMF2000.xlsx", sheet = 1) %>% dplyr::select(-Nacimiento, -Edad)
qxt_perm <- function(data, proy, genero){
      qxt <- matrix(0, nrow = nrow(data), ncol=proy)
      if(genero=="h"){
            base <- as.data.frame(data)[,1]
            factor <- as.data.frame(data)[,3]
      } else if (genero=="m"){
            base <- as.data.frame(data)[,2]
            factor <- as.data.frame(data)[,4]
      } else {
            cat("Género ingresado no válido")
      }
      for(j in 1:proy){
            qxt[,j] <- base * exp(-factor *j)
      }
      colnames(qxt) <- seq(2001, (2000+proy))
      return(qxt)
}

# Generación tablas PERMF (2018--2054)
TPh <- q1(qxt_perm(perm, proy=97, genero="h"))[1:106,-seq(1,4)]
TPm <- q1(qxt_perm(perm, proy=97, genero="m"))[1:106,-seq(1,4)]
rm(list=c("perm", "qxt_perm"))

# Carga de resultados modelos LCP, LL, CBD
load("Estimaciones.RData")


###   Seguro Capital Diferido ###
nominal <- function(tasa, num){
      res <- numeric(num)
      for(i in 0:(num-1)){
            res[i+1] <- (1+tasa)^i
      }
      return(res)
}

descuento <- function(tasa, num, tipo){
      vec <- numeric(num+1)
      for(i in 0:num){
            vec[i+1] <- 1/(1+tasa)^i
      }
      if(tipo=="pre"){
            res <- vec[-length(vec)]
      } else {
            res <- vec[-1]
      }
      return(res)
}

seg_dif <- function(edad, dur, capi, tipo, cre, tabla){
      tabla <- tabla[,-seq(1,13)] # filtramos periodo (2018-2097)
      prb <- cumprod(1 - diag(tabla[(edad+1):(edad+dur+1),]))
      psv <- prb[length(prb)]
      pxt <- c(1,prb[-length(prb)]) 
      res <- matrix(0, nrow=length(pxt), ncol=5)
      row.names(res) <- seq(edad, length.out = length(pxt))
      res[,1] <- pxt
      res[,2] <- nominal(cre,nrow(res))
      res[,3] <- pxt * nominal(cre,nrow(res)) * descuento(tipo, nrow(res),"pre")
      vap <- (capi * psv) / (1 + tipo)^dur
      prima <- vap / sum(res[,3])
      res[,4] <- round(prima * res[,2],2)
      res[,5] <- round(prima * res[,3],2)
      return(list(res,vap,prima))
}

seg_dif(edad=30, dur=20, capi=100000, tipo=0.02, cre=0.01, TPh)


### Seguro Mixto
# Cálculo de primas
smixto <- function(mqxt, edad, dur, inte, capi=1){
      mqxt <- mqxt[,-seq(1,13)] # se excluye la mortalidad anterior a 2018
      dif <- 105 - edad
      if(dif < dur){
            cat("Duración no adecuada, se realiza el cálculo con duración de ", dif, " años\n")
            qxt <- diag(mqxt[(edad+1):106,])
      } else {
            qxt <- diag(mqxt[(edad+1):(edad+dur),])
      }
      n <- length(qxt)
      pxta <- cumprod(1-qxt)
      st <- cumprod(rep(1/(1+inte), times=n))
      smuerte <- sum(st * c(1,pxta[-n]) * qxt)
      ssuper <- st[n] * pxta[n]
      return(round(capi*(smuerte+ssuper),2))
}

#smixto(LLh, edad=30, dur=5, inte=0.02, capi=100000)

iter_mixto <- function(tabla, iedad, idur, tipo, capital){
      res <- matrix(0, ncol=length(idur), nrow=length(iedad))
      colnames(res) <- idur
      rownames(res) <- iedad
      for(i in seq_along(iedad)){
            for(j in seq_along(idur)){
                  res[i,j] <- smixto(tabla, iedad[i], idur[j], tipo, capital)
            }
      }
      return(res)
}

iter_mixto(LLm, seq(30,60, by=5), seq(5,30, by=5), tipo=0.02, capital=100000)

# Plan de Pensiones (renta anual vitalicia)
aportacion <- function(mqxt, edad, inte, alfa=1, beta=1, prest=1 ,salida=65){
      mqxt <- mqxt[,-seq(1,13)] # se excluye la mortalidad anterior a 2018
      if(edad >= salida){
            cat("Edad de ingreso debe ser inferior a", salida, " años\n")
            break
      } else {
            qxt <- diag(mqxt[(edad+1):106,])
            pxta <- cumprod(1-qxt) # probabilidades de supervivencia acumuladas
            st <- cumprod(rep(1/(1+inte), times=length(pxta)))
      }
      dur <- salida - edad
      # flujos aportaciones
      p1 <- c(1,pxta[1:(dur-1)])
      st1 <- c(1,st[1:(dur-1)])
      axr <- sum(p1 * st1)
      # factor de actualizacion
      p2 <- pxta[dur]
      st2 <- st[dur]
      rEx <- p2 * st2
      # flujos prestaciones
      p3 <- pxta[(dur+1):length(pxta)]
      st3 <- st[(dur+1):length(st)]
      ax <- sum(p3 * st3)
      # calculo del aporte
      val <- (rEx * ax)/axr
      return(round(prest*val,2))
}

#aportacion(TPm, edad=60, inte=0.04)
  
iter_apor <- function(tabla1, tabla2, iedad, inte, alfa, beta, capital=1){
      res <- matrix(0, nrow = length(iedad), ncol=2)
      for(i in seq_along(iedad)){
            res[i,1] <- aportacion(tabla1, edad=iedad[i], inte, alfa, beta, capital)
            res[i,2] <- aportacion(tabla2, edad=iedad[i], inte, alfa, beta, capital)
      }
      rownames(res) <- iedad
      return(res)
}

iter_apor(TPh, LLh, seq(30,60,by=5), inte=0.03, alfa=1, beta=1, capital=10000)


# Cálculo de la reserva matemática (método retrospectivo)
reserva <- function(mqxt, edad, inte, alfa=1, beta=1, salida=65){
      
}

(Sys.time() - start)




