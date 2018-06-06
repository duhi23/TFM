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

pasem <- read_excel("Datos/Tablas_PERMF2000.xlsx", sheet = 3) %>% dplyr::select(-Edad)
rep_col <- function(vec, col){
      res <- matrix(0, nrow=length(vec), ncol=col)
      for(j in 1:col){
            res[,j] <- vec
      }
      return(res)
}

PSh <- q1(rep_col(pasem[[1]], col=93))
PSm <- q1(rep_col(pasem[[2]], col=93))

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
      prb <- cumprod(1 - diag(tabla[(edad+1):(edad+dur),]))
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
      return(list(res,vap, sum(res[,3]),prima))
}


tab1 <- seg_dif(edad=45, dur=20, capi=100000, tipo=0.02, cre=0.01, PSh)
tab2 <- seg_dif(edad=45, dur=20, capi=100000, tipo=0.02, cre=0.01, LLh)

dif_tabla <- function(tabla1, tabla2){
      res <- matrix(0, nrow=nrow(tabla1[[1]]), ncol=6)
      row.names(res) <- row.names(tabla1[[1]])
      res[,1] <- tabla1[[1]][,4] # cobrado
      res[,2] <- tabla1[[1]][,5] # cobrado descontado
      res[,3] <- tabla2[[1]][,4] # simulado
      res[,4] <- tabla2[[1]][,5] # simulado descontado
      res[,5] <- res[,1] - res[,3]
      res[,6] <- res[,2] - res[,4]
      rlong <- sum(res[,6])
      return(list(res,rlong))
}

dif_tabla(tab1, tab2)

# Simulación riesgo de longevidad
sim_long <- function(vec1, vec2, tab1, tab2){
      res <- matrix(0, nrow=length(vec1), ncol=length(vec2))
      for(i in 1:length(vec1)){
            for(j in 1:length(vec2)){
                  if(vec1[i] + vec2[j] <= 65){
                        tab01 <- seg_dif(edad=vec1[i], dur=vec2[j], capi=100000, tipo=0.06, cre=0.01, tab1)
                        tab02 <- seg_dif(edad=vec1[i], dur=vec2[j], capi=100000, tipo=0.06, cre=0.01, tab2)
                        res[i,j] <- dif_tabla(tab01, tab02)[[2]]
                  } else {
                        res[i,j] <- 0
                  }
            }
      }
      return(res)
}

sim_long(seq(30,60,by=5), seq(5,35, by=5), PSm, LLm)
sum(sim_long(seq(30,60,by=5), seq(5,35, by=5), PSm, LLm))



# Plan de Pensiones (renta anual vitalicia)
aportacion <- function(mqxt1, mqxt2, edad, inte, alfa=0.01, beta=0.01,salida=65, prest=1){
      mqxt1 <- mqxt1[,-seq(1,13)] # se excluye la mortalidad anterior a 2018
      mqxt2 <- mqxt2[,-seq(1,13)] # se excluye la mortalidad anterior a 2018
      pxt1 <- cumprod(1 - diag(mqxt1[(edad+1):106,]))
      pxt2 <- cumprod(1 - diag(mqxt2[(edad+1):106,]))
      des <- descuento(inte, length(pxt1)+1, "pre")
      # Tabla aportaciones
      pini <- c(1,pxt1[1:(salida-edad-1)])
      tabini <- matrix(0, nrow=length(pini), ncol=5)
      colnames(tabini) <- c("tpx", "nom", "va", "vaa", "vaax")
      tabini[,1] <- pini
      tabini[,2] <- prest * nominal(alfa, length(pini))
      tabini[,3] <- tabini[,2] * des[1:length(pini)]
      tabini[,4] <- tabini[,1] * tabini[,3]
      tabini[,5] <- c(1,pxt2[1:(salida-edad-1)]) * tabini[,3]
      # Tabla prestaciones
      pfin <- pxt1[-c(1:(salida-edad-1))]
      tabfin <- matrix(0, nrow=length(pfin), ncol=5)
      tabfin[,1] <- pfin
      tabfin[,2] <- nominal(beta, length(pfin))
      tabfin[,3] <- tabfin[,2] * des[-c(1:length(pini))]
      tabfin[,4] <- tabfin[,1] * tabfin[,3]
      tabfin[,5] <- pxt2[-c(1:(salida-edad-1))] * tabfin[,3]
      # Prima
      prima <- sum(tabini[,4])/sum(tabfin[,4])
      # Tabla final con todos los flujos +/-
      res <- rbind(tabini[,-c(1)], -prima*tabfin[,-c(1)])
      rlong <- round(apply(res[,c(3,4)],2,sum),2)
      return(list(tabini,tabfin,prima,res,rlong))
}

aportacion(PSh, LLh, edad=40, inte=0.01, prest=1000)


(Sys.time() - start)




