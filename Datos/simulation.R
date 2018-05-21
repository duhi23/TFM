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

qxt_perm(perm, 30, "h")


# Cálculo de primas
load("Estimaciones.RData")

smixto <- function(mqxt, edad, dur, inte){
      dif <- 105 - edad
      if(dif < dur){
            cat("Duración no adecuada, se realiza el cálculo con duración de ", dif, " años\n")
            val <- diag(mqxt[(edad+1):106,])
            qxt <- ifelse(val>1,1,val)
      } else {
            cat("Bien\n")
            val <- diag(mqxt[(edad+1):(edad+dur),])
            qxt <- ifelse(val>1,1,val)
      }
      n <- length(qxt)
      pxta <- cumprod(1-qxt)
      st <- cumprod(rep(1/(1+inte), times=n))
      smuerte <- sum(st * c(1,pxta[-n]) * qxt)
      ssuper <- st[n] * pxta[n]
      #return(list(qxt,pxta,st,smuerte,ssuper))
      return(smuerte+ssuper)
}


smixto(LLh, 50, 40, 0.04)






(Sys.time() - start)




