#This function calculates tissue:plasma partition coefficients according to Schmitt https://www.sciencedirect.com/science/article/pii/S0887233307002573?via%3Dihub

.libPaths("lib")
library(dplyr)

calcKp_Schmitt <- function(logP, pKa, fup, type = 1, dat){
  #logMA is the log of membrane affinity = phosphatidylcholine:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)
  
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Schmitt.csv")
  #dat <- read.csv("/data/internship-summer-2018/data/PKSim_tissue_comp_Schmitt.csv") #Physiology from PK-Sim
  
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark.csv") #Human physiology from Ruark et al.
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_rat_for_R&R.csv") #Rat physiology from Ruark et al.
  
  #dat <- read.csv("/data/internship-summer-2018/data/unified_tissue_comp.csv") # Unified physiology from Ruark, P&T, R&R, and PK-Sim
  ##dat <- read.csv("../data/unified_tissue_comp.csv")
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for adipose and RBCs
  
  
  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  #K_protein <- 0.163+0.0221*K_n_pl    #protein:water partition; Schmitt, Walter (2008)
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5
  pH <- dat_all$pH
  alpha <- 1e-3  #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)
  

  W <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH-pKa),
              #3-monoprotic base
              10^(pKa-pH),
              #4-diprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2]), 
              #5-diprotic base
              10^(pKa[2]-pH)+10^(pKa[1]+pKa[2]-2*pH), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH)+10^(pH-pKa[1]),
              #7-triprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2])+10^(3*pH-pKa[1]-pKa[2]-pKa[3]))
#             #8-triprotic base
#             10^(pKa[1]+pKa[2]+pKa[3]-3*pH),
  
  
  
  if(type==1 | type==2 | type==4 | type==7){ # neutral, monoprotic acid, diprotic acid, triprotic acid
    K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha)
    K_a_pl <- K_n_pl*((1/(1+W))+0.05*(1-(1/(1+W)))) 
  }
  else if(type==3){ # monoprotic base
    K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha) 
    K_a_pl <- K_n_pl*((1/(1+W))+20*(1-(1/(1+W)))) 
  }
  else if(type==5){ # diprotic base 
    F1 <- (1/(1+10^(pKa[1]-pH)))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
    K_a_pl <- K_n_pl*(F1*F2 + 20*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
  }
  else if(type==6){ # monoprotic acid monoprotic base (acid comes first)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
    K_a_pl <- K_n_pl*(F1*F2 + 0.05*(1-F1)*F2 + 20*(1-(F2))*F1 + (1-F1)*(1-F2))
  }
  else if(type==8){ # triprotic base 
    F1 <- (1/(1+10^(pKa[1]-pH)))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 20*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3)) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
  }
  else if(type==9){ # diprotic acid monoprotic base (first two are acid)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pH-pKa[2])))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*((1-F1)*F2*F3 + F1*(1-F2)*F3 + (1-F1)*(1-F2)*F3) + 20*F1*F2*(1-F3) + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
  }
  else if(type==10){ # diprotic base monoprotic acid (first one is acid)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*(1-F1)*F2*F3 + 20*(F1*(1-F2)*F3 + F1*F2*(1-F3) + F1*(1-F2)*(1-F3)) + (1-F1)*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(1-F2)*(1-F3))
  }
  
  kp <- (dat_all$f_water+(K_n_l*dat_all$f_n_l)+(K_n_pl*dat_all$f_n_pl)+(K_a_pl*dat_all$f_a_pl)+(K_protein*dat_all$f_proteins))*fup

  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name
  
  return(Kp)
}