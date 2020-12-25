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

library(openxlsx)


library(rio)
data_list <- import_list("prelim40.xlsx")
df <- vector(mode="list", length=(length(data_list)-length(whichzero)))
# get only columns that have Confidence Interval
for (i in 1:(length(data_list)-length(whichzero))){
  df[[i]] <- data_list[[i]][,1:10]
}

  
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
newx <- vector(mode="list", length=length(data_list))
for (i in 1:length(data_list)){
newx[[i]] <- df[[i]] %>% mutate(id=1:n()) %>% 
  pivot_longer(-id) %>%
  separate_rows(value,sep = '\\(') %>%
  separate_rows(value,sep = ',') %>%
  mutate(value=trimws(gsub(')','',value))) %>%
  group_by(id,name) %>%
  mutate(name=ifelse(row_number()==1,paste0(name,'.est'),
                     ifelse(row_number()==2,paste0(name,'.low'),paste0(name,'.upper')))) %>%
  pivot_wider(names_from = name,values_from=value) %>% ungroup()
}
for (i in 1:length(data_list)){
  for(j in 1:30){
    newx[[i]][,j] <- as.numeric(unlist(newx[[i]][,j]))
  }
}



for (i in 1:length(data_list)){
newx[[i]] <- newx[[i]][,2:31]
}

for (i in 1:length(data_list)){
     zeros[i] <- sum(newx[[i]]==0)
     whichzero <- which(zeros!=0)
  }
newx <- newx[-whichzero]
#accuracy for each variable for each trial
x<- vector(mode="list", length=(length(data_list)-length(whichzero)))
for (i in 1:(length(data_list)-length(whichzero))){
for(k in 1:length(colnames(df[[1]]))){
  x[[i]][k] <- 100*(sum((newx[[i]][,(3*k)-2] > newx[[i]][,(3*k)-1])==TRUE & (newx[[i]][,(3*k)-2] < newx[[i]][,3*k])==TRUE ))/52
}
}

a<- vector(mode="list", length=(length(data_list)-length(whichzero)))

#get convergence for one time parameters across all trials
for (i in 1:(length(data_list)-length(whichzero))){
for(j in 1:7){
a[[j]][i] <- ((x[[i]][j] > 99.9) & (x[[i]][j] < 100.1) )
 tauOcon <- 100*sum(a[[1]]==TRUE)/(length(data_list)-length(whichzero))
 tauPcon <- 100*sum(a[[2]]==TRUE)/(length(data_list)-length(whichzero))
 tauScon <- 100*sum(a[[3]]==TRUE)/(length(data_list)-length(whichzero))
 thetacon <- 100*sum(a[[4]]==TRUE)/(length(data_list)-length(whichzero))
 rhocon <- 100*sum(a[[5]]==TRUE)/(length(data_list)-length(whichzero))
 gammacon <- 100*sum(a[[6]]==TRUE)/(length(data_list)-length(whichzero))
 Sp0con <- 100*sum(a[[7]]==TRUE)/(length(data_list)-length(whichzero))
 
 }
}

library(htmlTable)
library(magrittr)
tbl <- matrix(c(tauOcon, tauPcon, tauScon, thetacon, rhocon, gammacon,Sp0con),
       ncol = 1,
       dimnames = list(c("tauO", "tauP","tauS","theta","rho","gamma","Sp0"),
                       c("% Convergence from N= Trials"))) %>% 
  htmlTable
tbl

b<- vector(mode="list", length=(length(data_list)-length(whichzero)))
for (i in 1:3){
  for(j in 1: (length(data_list)-length(whichzero))){
  b[[1]][j] <- c(x[[j]][8]) #infected
  b[[2]][j] <- c(x[[j]][9]) #susceptibles
  b[[3]][j] <- c(x[[j]][10]) #lambda
  
infectedmean <-round(mean(b[[1]]), digits=1)
infectedsd <- round(sd(b[[1]]), digits=1)
infectedstring <- sprintf('%s +/- %s', infectedmean, infectedsd)  #put in CI 


susceptiblemean <- round(mean(b[[2]]), digits=1)
susceptiblesd <- round(sd(b[[2]]), digits=1)
susceptiblestring <- sprintf('%s +/- %s', susceptiblemean, susceptiblesd)  #put in CI 


lambdamean <- round(mean(b[[3]]), digits=1)
lambdasd <- round(sd(b[[3]]),digits=1)
lambdastring <- sprintf('%s +/- %s', lambdamean, lambdasd)  #put in CI 

}
}

tbl2 <- matrix(c(infectedstring, susceptiblestring, lambdastring),
              ncol = 1,
              dimnames = list(c("Infected", "Susceptible", "Lambda"),
                              c("% Mean +/- SD , N=35 "))) %>% 
  htmlTable

zeros <- vector()
for (i in 1:(length(data_list)-length(whichzero))){
  zeros[i] <- sum(newx[[i]]==0)
  whichzero <- which
}
