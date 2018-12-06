#calculate tissue:plasma partition coefficients based on: Rodgers and Rowland http://jpharmsci.org/article/S0022-3549(16)31789-0/fulltext and http://jpharmsci.org/article/S0022-3549(16)32034-2/fulltext
#tissue data are derived from rats

.libPaths("lib")
library(dplyr)

calcKp_RR <- function(logP, pKa=0, fup, BP=1, type=1, dat){
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_P_R&R.csv") # For use when comparing with Peyret et al. (rat data)
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_R&R.csv") # Rat data from paper
  #dat <- read.csv("/data/internship-summer-2018/data/PKSim_tissue_comp_RR.csv") # Data taken from PK-Sim
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_for_R&R.csv") # human f_nl, f_n_pl, and f_n_ap from Ruark, rest from PK-Sim for R&R
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_rat_for_R&R.csv") #rat f_nl, f_n_pl, and f_n_ap from Ruark, rest from PK-Sim for R&R
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_rat_for_R&R.csv") #Rat physiology from Ruark et al.
  
  #dat <- read.csv("/data/internship-summer-2018/data/unified_tissue_comp.csv") # Unified physiology from Ruark, P&T, R&R, and PK-Sim
  #dat <- read.csv("../data/unified_tissue_comp.csv")
  
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma
  
  
  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit
  
  
  #Calculate Kpu values
  #Kpu_bc <- (HCT - 1 + BP)/(HCT*fup)
  
  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT)
  
  X <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_IW-pKa),
              #3-monoprotic base
              10^(pKa-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_IW)+10^(pH_IW-pKa[1]),  
              #7-triprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))       
  
  Y <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_P-pKa),
              #3-monoprotic base
              10^(pKa-pH_P), 
              #4-diprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_P)+10^(pH_P-pKa[1]),  
              #7-triprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_P)+10^(pKa[3]+pka[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))       
  
  Z <- switch(type,
              #1-neutral
              1,   
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa-pH_RBC), 
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),  
              #7-triprotic acid
              1,  
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC)) 
  
  
  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
  
  #PK-Sim requires the following:
  #Ka_PR <- max(0,Ka_PR) 
  #Ka_AP <- max(0,Ka_AP)
  
  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
  
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  
  
  # Calculate Kpu (reported in R&R paper)
  # if(type_calc==1){  #moderate to strong bases
  #   Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))  #non lipid
  #   Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))  #lipid
  # }else if(type_calc==2){   #acidic and zwitterions
  #   Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))  #non lipid
  #   Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y)) #lipid
  # }else{  #neutrals
  #   Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))  #non lipid
  #   Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))  #lipid
  # }
  # 
  # nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  # nms_all <- paste("Kp", nms_all, sep="")
  # nms <- c(nms_all,"Kpad")
  # Kp <- as.list(c(Kp_all,Kp_ad))
  # names(Kp) <- nms
  
  
  
  # Multiply by fup to get Kp rather than Kpu
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*fup  #lipid
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup #lipid
  }else{  #neutrals
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))*fup  #lipid
  }


  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms
  
  # nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  # nms_all <- paste("Kp", nms_all, sep="")
  # nms <- c(nms_all,"Kpad")
  # Kp <- as.list(c(Kp_all,Kp_ad))
  # names(Kp) <- nms
  
  return(Kp)
}




