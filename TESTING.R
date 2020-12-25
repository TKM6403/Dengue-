rm(list=ls())
setwd("~/Desktop/stats")
library(reshape2)
library(ggplot2)
library(spatialpred)
library(dengueThailand)
library(foreach)
library(ggplot2)
library(dplyr)
library(grid)
library(readr)
library(gridExtra)
library(cowplot)
library(reshape2)
library(rjags)
library(wesanderson)
library(cowplot)
library(scales)
library(pROC)
library(ResourceSelection)
library(readxl)
library(relaimpo)
library(ridge)
library(car)
library(mgcv)
library(glmnet)
library(caret)
library(tidyverse)

#tauO(1), tauOJAGS(2-6), tauP(7), tauPJAGS(8-12)
#tauS(13), tauSJAGS(14-18), theta(19), thetaJAGS(20-24), rho(25), rho(26-30), gamma(31), gamma(32-36)
#Sp0(37), Sp0(38-42), CASES(43-94), CASES(95-354), I(355-406), I(407-666), Sp(667-718), Sp(719-978)
#lambda(979-1030), lambda(1031-1290), C(1291-1342)

ntrials <- 1
M <- matrix(0, nrow=ntrials, ncol=1446)
#before for loop, go through one case and make it work. Then generate 100 trials of different combinations
#need to run from sim Dengue

#R0=thetainit[1:100]/gammainit[1:100]

#full model
run_mod <- function(Sp0inp, I1inp, tauOinp, tauPinp, tauSinp, thetainp, gammainp, rhoinp){

population<- c(seq(4167000,4266000
                   ,(4266000-4167000)/52))
#births per population obtained from website
births <- 0.088*population
year<-c(1:52)

#Number of time observations

J <- 1 # Number of Serotypes
L<-length(seq(from=1, to=52, by=1)) #Number of time observations

numInf<-1 # Number of times you can get infected

if(J==2){
  serotypes<-c("DEN-1", "DEN-2")# ##Names of the serotypes
} else if(J==1){
  serotypes <- c("DEN-1")
} else if(J==4){ 
  serotypes<-c("DEN-1", "DEN-2", "DEN-3", "DEN-4") ##Names of the serotypes
}else{print("THIS DOESNT WORK")}

#initialize everything
I<-array(NA, c(numInf,L,J), dimnames=list(c("Primary"), c(year), c(serotypes))) # Infecteds 
Sp<-list(array(NA, c(1,L,1), dimnames=list(c("Primary"), c(year), c("DEN"))))

C<-array(NA, c(1,L,J), dimnames=list(c("Cases"), c(year), c(serotypes))) # Cases
K<-array(NA, c(1,L,J), dimnames=list(c("Cases"), c(year), c(serotypes))) # Subsample of serotyped cases
CASES<-array(NA, c(1,L,1), dimnames=list(c("Cases"), c(year))) # Cases
## Set the observed varibles that we know
B<-array(births, c(1,L,1), dimnames = list(c("Births"), c(year))) # Births
N<-array(population, c(1,L,1), dimnames = list(c("Pop"), c(year))) # Population
rho<-array(c(rhoinp), c(numInf, 1, 1), dimnames = list(c("Primary"), c("rho"))) # Reporting rate
CVar<-array(NA, c(1,L,J), dimnames=list(c("CasesCV"), c(year), c(serotypes))) ## Array of coefficients of variation for cases
IVar<-array(NA, c(numInf,L,J), dimnames=list(c("Primary"), c(year), c(serotypes))) ## Array of coefficients of variation for Infections
SpVar<-list(array(NA, c(1,L,1), dimnames=list(c("Primary"), c(year), c("DEN"))))


#rho=0.005
## Initialize the transition parameters
lambda<-array(NA, c(1,L,J), dimnames = list(c("FOI"), c(year), c(serotypes))) # FOI
tauO=tauOinp
tauP= tauPinp
tauS= tauSinp
theta= thetainp
gamma= gammainp




I[1] <-I1inp	 #figure out in JAGS
lambda[1] <- 1-exp(-theta*(I[1]/N[1]))
Sp0 <- Sp0inp	 #figure out in JAGS
#Sp0 <-runif(1,1600000,2000000)
IVar[1] <- (lambda[1] * Sp0)*tauP

SpVar[[1]][1] = (Sp0-I[1]+B[1])*tauS
Sp[[1]][1] <- rnorm(1,Sp0-I[1]+B[1], SpVar[[1]][1] )

CVar[1] <- (rho*I[1])*tauO
C[1] <- rnorm(1,rho*I[1], CVar[1])
CASES[1] <-C[1]




################################## LOOP OVER THE REST OF THE VALUES ################################## 
for(l in 2:L){
  for(j in 1:J){
    lambda[1,l,j] <- 1- exp(-theta*(I[1,l-1,1])/N[l-1])
  }
  
  
  ####### INFECTIONS ########
  ## Calculate the IVar
  for(j in 1:J){
    IVar[1,l,j] = (lambda[1,l,j]*Sp[[1]][1,l-1,1]-gamma*I[1,l-1,j])*tauP
    # IVar[2,l,j] = (sum(lambda[1,l,j]*S[[2]][1,l-1,])-lambda[1,l,j]*S[[2]][1,l-1,j])*tauP
  }
  
  
  ## Calculate infections at time t=1
  for(j in 1:J){
    I[1,l,j] = rnorm(1, lambda[1,l,j]*Sp[[1]][1,l-1,1]-gamma*I[1,l-1,j], IVar[1,l,j])
    # I[2,l,j] = rnorm(1, sum(lambda[1,l,j]*S[[2]][1,l-1,])-lambda[1,l,j]*S[[2]][1,l-1,j], IVar[1,l,j])
  }
  
  
  ####### SUSCEPTIBLES ########
  ## Calculate the SVar
  SpVar[[1]][1,l,1] = (Sp[[1]][1,l-1,1] -  sum(I[1,l,]) + B[1,l,1])*tauS
  
  ## Calculate susceptibles at time t=1
  Sp[[1]][1,l,1] <- rnorm(1, Sp[[1]][1,l-1,1] -  sum(I[1,l,]) + B[1,l,1], SpVar[[1]][1,l,1])
  
  
  ## Calculate susceptibles and SVar at time t=1 for secondary infections
  temp<-matrix(NA, nrow=J, ncol=J)
  ####### CASES ########
  ## Calculate the CVar
  for (j in 1:J){
    CVar[1,l,j] = (rho[,1,1]*I[,l,j])*tauO
  }
  
  
  ## Calculate cases at t=1
  
  for (j in 1:J){
    C[1,l,j] <- rnorm(1, (rho[,1,1]*I[,l,j]), CVar[1,l,j]) ## Cases of serotype J
  }
 
   CASES[l] <- sum(C[1,l,], na.rm=TRUE) ## Total number of cases
  }

return(list(lambda,I,Sp,C,CASES,B,N))
}


#Sp0inp, I1inp, tau0inp, tauPinp, tauSinp, thetainp, gammainp, rho1inp
i <- 0
while(i <= ntrials){
Sp0inp= runif(1, min=0, max=4167000) #Sp0
I1inp= runif(1, min=0, max=100000) #I[1]
tauOinp= runif(1, min=0, max=1) #tauO
tauPinp=runif(1, min=0, max=1) #tauP
tauSinp= runif(1, min=0, max=1) #tauS
thetainp= runif(1, min=1, max=5) #theta
gammainp= runif(1, min=0.7, max=1) #gamma
rhoinp=runif(1, min=0, max=1) #rho

output <- run_mod(Sp0inp, I1inp,tauOinp,tauPinp,tauSinp,thetainp,gammainp,rhoinp)

#make sure the trial was fully successful. There would be NAs if something went wrong
if(sum(is.na(output[[1]]>0)+is.na(output[[2]]>0)+is.na(output[[3]][[1]]>0)+is.na(output[[4]]>0)+is.na(output[[5]]>0)+is.na(output[[6]]>0)+is.na(output[[7]]>0))==FALSE){
M[i,37] <- Sp0inp
M[i,355] <- I1inp
M[i,1] <- tauOinp
M[i,7] <- tauPinp
M[i,13] <- tauSinp
M[i,19] <- thetainp
M[i,31] <- gammainp
M[i,25] <- rhoinp
M[i,979:1030] <- output[[1]][1:52]# generated lambda 1-52 for each trial
M[i,355:406] <- output[[2]][1:52] #generated I 1-52 for each trial
M[i,667:718] <- output[[3]][[1]][1:52] #generated Sp 1-52 for each trial
M[i,1291:1342] <- output[[4]][1:52]  #generated C 1-52 for each trial
M[i,43:94] <- output[[5]][1:52]
M[i, 1343:1394] <- output[[6]][1:52]
M[i, 1395:1446] <- output[[7]][1:52]#generated CASES 1-52 for each trial
i <- i+1
}
#lambda,I,Sp,C,CASES
#tauO(1), tauOJAGS(2-6), tauP(7), tauPJAGS(8-12)
#tauS(13), tauSJAGS(14-18), theta(19), thetaJAGS(20-24), rho(25), rho(26-30), gamma(31), gamma(32-36)
#Sp0(37), Sp0(38-42), CASES(43-94), CASES(95-354), I(355-406), IJAGS(407-666), Sp(667-718), Sp(719-978)
#lambda(979-1030), lambda(1031-1290), C(1291-1342), Births (1343-1394), Pop(1395-1446)


}
for (i in 1:10){
  print(M[i,31])
}

###JAGS
jagsmodel <- function(births,population, cases){

JAGS.sim.data<-NULL
JAGS.sim.data$T<-length(seq(from=1, to=52, by=1))
JAGS.sim.data$B<-c(births)
JAGS.sim.data$N<-c(population)
JAGS.sim.data$C <- cases

params <- c("CASES","Sp", "tauS", "tauO", "tauP","theta","rho1", "gamma", "lambda", "Sp0", "I" ,"theta")

jmod <- jags.model("simdengue.jags", 
                   data=JAGS.sim.data,
                   n.chains=1,
                   inits = list(tauS=0.1),
                   n.adapt=100000)

jsamp <- coda.samples(jmod,
                      variable.names=params,
                      n.iter=100000,
                      thin=10)

fit <- summary(as.mcmc.list(jsamp))

fitquant <- as.data.frame(fit$quantiles)
return(fitquant)
}

for (i in 1:ntrials){
  for (j in 1:52){
births <- list()
births[[i]]= M[i, 1343:1394]

population<- list()
population[[i]]= M[i, 1395:1446] 

cases <- list()
cases[[i]]=M[i,1291:1342]

output <- list()
output[[i]] <- jagsmodel(births[[i]][1:52],population[[i]][1:52],cases[[i]][1:52])
#put in all CASES that JAGS makes
M[i,(5*j+90):(5*j+94)] = c(output[[i]][j,1:5]$`2.5%`,output[[i]][j,1:5]$`25%`,output[[i]][j,1:5]$`50%`,output[[i]][j,1:5]$`75%`,output[[i]][j,1:5]$`97.5%`)
#put in all I that JAGS makes
M[i,(5*j+402):(5*j+406)]= c(output[[i]][j+52,1:5]$`2.5%`,output[[i]][j+52,1:5]$`25%`,output[[i]][j+52,1:5]$`50%`,output[[i]][j+52,1:5]$`75%`,output[[i]][j+52,1:5]$`97.5%`)
#put in all Sp
M[i,(5*j+714):(5*j+718)]= c(output[[i]][j+104,1:5]$`2.5%`,output[[i]][j+104,1:5]$`25%`,output[[i]][j+104,1:5]$`50%`,output[[i]][j+104,1:5]$`75%`,output[[i]][j+104,1:5]$`97.5%`)
#put in Sp0
M[i,38:42]= c(output[[i]][157,1:5]$`2.5%`,output[[i]][157,1:5]$`25%`,output[[i]][157,1:5]$`50%`,output[[i]][157,1:5]$`75%`,output[[i]][157,1:5]$`97.5%`)
#put in gamma
M[i,32:36]= c(output[[i]][158,1:5]$`2.5%`,output[[i]][158,1:5]$`25%`,output[[i]][158,1:5]$`50%`,output[[i]][158,1:5]$`75%`,output[[i]][158,1:5]$`97.5%`)
#put in lambda
M[i,(5*j+1026):(5*j+1030)]= c(output[[i]][j+158,1:5]$`2.5%`,output[[i]][j+158,1:5]$`25%`,output[[i]][j+158,1:5]$`50%`,output[[i]][j+158,1:5]$`75%`,output[[i]][j+158,1:5]$`97.5%`)
#put in rho
M[i,26:30]= c(output[[i]][211,1:5]$`2.5%`,output[[i]][211,1:5]$`25%`,output[[i]][211,1:5]$`50%`,output[[i]][211,1:5]$`75%`,output[[i]][211,1:5]$`97.5%`)
#put in tauO
M[i,2:6]= c(output[[i]][212,1:5]$`2.5%`,output[[i]][212,1:5]$`25%`,output[[i]][212,1:5]$`50%`,output[[i]][212,1:5]$`75%`,output[[i]][212,1:5]$`97.5%`)
#put in tauP
M[i,8:12]= c(output[[i]][213,1:5]$`2.5%`,output[[i]][213,1:5]$`25%`,output[[i]][213,1:5]$`50%`,output[[i]][213,1:5]$`75%`,output[[i]][213,1:5]$`97.5%`)
#put in tauS
M[i,14:18]= c(output[[i]][214,1:5]$`2.5%`,output[[i]][214,1:5]$`25%`,output[[i]][214,1:5]$`50%`,output[[i]][214,1:5]$`75%`,output[[i]][214,1:5]$`97.5%`)
#put in theta
M[i,20:24]= c(output[[i]][215,1:5]$`2.5%`,output[[i]][215,1:5]$`25%`,output[[i]][215,1:5]$`50%`,output[[i]][215,1:5]$`75%`,output[[i]][215,1:5]$`97.5%`)
}
  }

#lambda,I,Sp,C,CASES
#tauO(1), tauOJAGS(2-6), tauP(7), tauPJAGS(8-12)
#tauS(13), tauSJAGS(14-18), theta(19), thetaJAGS(20-24), rho(25), rho(26-30), gamma(31), gamma(32-36)
#Sp0(37), Sp0(38-42), CASES(43-94), CASES(95-354), I(355-406), IJAGS(407-666), Sp(667-718), Sp(719-978)
#lambda(979-1030), lambda(1031-1290), C(1291-1342), Births (1343-1394), Pop(1395-1446)


#need to make sure that 52 CASES are right
#52 I are right
#52 Sp are right
#Sp0 is right
#Gamma is right
#52 lambda are right
#rho is right
#tauO, tauP, tauS are right
#theta is right

#MAKE A TABLE
# Each entry should have value for that trial and that time 
#And then parantheses () with the Confidence Interval from JAGS.
library(data.table)
library(dplyr)
library(formattable)
library(tidyr)
library(glue)
Mround <- round(M, digits=4)
O <- vector(mode="list", length=ntrials)
#each element in list is a data frame with all the parameter values across time 
#columns will be the different parameters, rows will be the 52 weeks
for (i in 1:ntrials){
  
  #initialize data frame for each trial
  O[[i]] <- data.frame(matrix(ncol = 14, nrow = 52))
  # Rename column 
  names(O[[i]])[names(O[[i]]) == "X1"] <- "tauO"
  names(O[[i]])[names(O[[i]]) == "X2"] <- "tauP"
  names(O[[i]])[names(O[[i]]) == "X3"] <- "tauS"
  names(O[[i]])[names(O[[i]]) == "X4"] <- "theta"
  names(O[[i]])[names(O[[i]]) == "X5"] <- "rho"
  names(O[[i]])[names(O[[i]]) == "X6"] <- "gamma"
  names(O[[i]])[names(O[[i]]) == "X7"] <- "Sp0"
  names(O[[i]])[names(O[[i]]) == "X8"] <- "CASES"
  names(O[[i]])[names(O[[i]]) == "X9"] <- "I"
  names(O[[i]])[names(O[[i]]) == "X10"] <- "Sp"
  names(O[[i]])[names(O[[i]]) == "X11"] <- "lambda"
  names(O[[i]])[names(O[[i]]) == "X12"] <- "C"
  names(O[[i]])[names(O[[i]]) == "X13"] <- "Births"
  names(O[[i]])[names(O[[i]]) == "X14"] <- "Pop"
  O[[i]][,1]<-Mround[i,1] #tauO
  O[[i]][,1][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,1]), Mround[i,2], Mround[i,6]);O[[i]][,1]  #put in CI 
  
  O[[i]][,2]<-Mround[i,7] #tauP
  O[[i]][,2][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,2]), Mround[i,8], Mround[i,12]);O[[i]][,2]  #put in CI 
  
  O[[i]][,3]<-Mround[i,13] #tauS
  O[[i]][,3][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,3]), Mround[i,14], Mround[i,18]);O[[i]][,3]  #put in CI 
  
  O[[i]][,4]<-Mround[i,19] #theta
  O[[i]][,4][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,4]), Mround[i,20], Mround[i,24]);O[[i]][,4]  #put in CI 
  
  O[[i]][,5]<-Mround[i,25] #rho
  O[[i]][,5][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,5]), Mround[i,26], Mround[i,30]);O[[i]][,5]  #put in CI 
  
  O[[i]][,6]<-Mround[i,31] #gaMroundma
  O[[i]][,6][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,6]), Mround[i,32], Mround[i,36]);O[[i]][,6]  #put in CI 
  
  O[[i]][,7]<-Mround[i,37] #SpO
  O[[i]][,7][] <- sprintf('%s (%s, %s)', unlist(O[[i]][,7]), Mround[i,38], Mround[i,42]);O[[i]][,7]  #put in CI 

  O[[i]][,8] <-Mround[i,43:94]#CASES
  
  f8 <- function(trial,cell) {
    paste(cell, " (", Mround[trial,5*(which(O[[trial]][,8]==cell))+90], ",", Mround[trial,5*(which(O[[trial]][,8]==cell))+94], ")")
    }

  x8 <- vector(mode = "list", length = ntrials)
  for(j in 1:52){
     x8[[i]][j] <- f8(i, O[[i]][,8][j])
     }
  O[[i]][,8] <- x8[[i]]
   
 
  
   f9 <- function(trial, cell) {
    paste(cell, " (", Mround[trial,5*(which(O[[trial]][,9]==cell))+402], ",", Mround[trial,5*(which(O[[trial]][,9]==cell))+406], ")")
  }
  O[[i]][,9] <-Mround[i,355:406]#CASES
  x9 <- vector(mode = "list", length = ntrials)
  for(j in 1:52){
    x9[[i]][j] <- f9(i, O[[i]][,9][j])
  }
  O[[i]][,9] <- x9[[i]]
  
  
  f10 <- function(trial,cell) {
    paste(cell, " (", Mround[trial,5*(which(O[[trial]][,10]==cell))+714], ",", Mround[trial,5*(which(O[[trial]][,10]==cell))+718], ")")
  }
  O[[i]][,10] <-Mround[i,667:718]#CASES
  x10 <- vector(mode = "list", length = ntrials)
  for(j in 1:52){
    x10[[i]][j] <- f10(i, O[[i]][,10][j])
  }
  O[[i]][,10] <- x10[[i]]
  
  f11 <- function(trial,cell) {
    paste(cell, " (", Mround[trial,5*(which(O[[trial]][,11]==cell))+1026], ",", Mround[trial,5*(which(O[[trial]][,11]==cell))+1030], ")")
  }
  O[[i]][,11] <-Mround[i,979:1030]#CASES
  x11 <- vector(mode = "list", length = ntrials)
  for(j in 1:52){
    x11[[i]][j] <- f11(i, O[[i]][,11][j])
  }
  O[[i]][,11] <- x11[[i]]
  
  
  O[[i]][,12] <-Mround[i,1291:1342]#C

  O[[i]][,13] <-Mround[i,1343:1394]#Births
  
  O[[i]][,14] <-Mround[i,1395:1446]#Pop
  
  
   
}
#tauO(1), tauOJAGS(2-6), tauP(7), tauPJAGS(8-12)
#tauS(13), tauSJAGS(14-18), theta(19), thetaJAGS(20-24), rho(25), rho(26-30), gamma(31), gamma(32-36)
#Sp0(37), Sp0(38-42), CASES(43-94), CASES(95-354), I(355-406), IJAGS(407-666), Sp(667-718), Sp(719-978)
#lambda(979-1030), lambda(1031-1290), C(1291-1342), Births (1343-1394), Pop(1395-1446)

#eliminate random CASES column
for (i in 1:ntrials){
O[[i]]$CASES <- NULL
}

x <- vector()
for (i in 1:ntrials){
  x[i] <- paste("trial",i)
  
}
names(O) <- x
#figure out  better way to name O by trial
#figure out way to adapt writexl to have alltrials
library(writexl)
writexl::write_xlsx(O, "prelim5.xlsx")






