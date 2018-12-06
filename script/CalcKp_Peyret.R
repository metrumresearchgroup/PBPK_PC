# Calculate Kp using the Peyret et al. model (aka the unified algorithm)

.libPaths("lib")
library(dplyr)


calcKp_Pey <- function(logP, pKa=0, fup, BP=1, type=1){
  
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Peyret.csv")
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_rat_for_R&R.csv") #Rat physiology from Ruark et al.
  #dat <- read.csv("/data/internship-summer-2018/data/unified_tissue_comp.csv") # Unified physiology from Ruark, P&T, R&R, and PK-Sim
  ##dat <- read.csv("../data/unified_tissue_comp.csv")
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma","Adipose"))  #df for all tissues except for adipose and RBCs
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs")
  dat_plas <- dat %>% filter(tissue == "Plasma")
  
  pH_IW <- 7       #pH of intracellular tissue water (cells)
  pH_IF <- 7.4     #pH of interstitial fluid
  pH_RBC <- 7.22    #pH of blood cells
  pH_P <- 7.4      #pH of plasma
  P <- 10^(logP)   # octonal:water partition coeff  
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  #Ka <- 10^(-pKa)
  
  HCT <- 0.45      
  Kpu_bc <- (HCT - 1 + BP)/(HCT*fup) #RBC:plasma water conc ratio (P_{ew} in the Peyret paper)
  #K_rbc <- Kpu_bc*fup 
  
  W <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_IF-pKa),
              #3-monoprotic base
              10^(pKa-pH_IF),
              #4-diprotic acid
              10^(pH_IF-pKa[1])+10^(2*pH_IF-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IF)+10^(pKa[1]+pKa[2]-2*pH_IF), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_IF)+10^(pH_IF-pKa[1]),  
              #7-triprotic acid
              10^(pH_IF-pKa[1])+10^(2*pH_IF-pKa[1]-pKa[2])+10^(3*pH_IF-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_IF)+10^(pKa[3]+pKa[2]-2*pH_IF)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IF),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IF)+10^(pH_IF-pKa[1])+10^(2*pH_IF-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IF-pKa[1])+10^(pKa[3]-pH_IF)+10^(pKa[2]+pKa[3]-2*pH_IF))     
  
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
  
  
  # Only for acids, zwitterions, weak bases, and neutrals, protein:water PC plasma
  plas_alb <- 0.029 #fraction of albumin in plasma
  plas_lip <- 0.0006 #fraction of lipoprotein in plasma (typo in Table 3 in paper, listed later as 0.0006)
  
  
  #Only for acids and weak bases
  Pprwp_alb <- (1/fup - 1 - P*(dat_plas$f_n_l+0.3*dat_plas$f_n_pl)/(1 + Y))*(1/(plas_alb))
  
  #Only for neutrals
  Pprwp_lip <- (1/fup - 1 - P*(dat_plas$f_n_l+0.3*dat_plas$f_n_pl)/(1 + Y))*(1/(plas_lip))
  
  
  # Only for pka>7, acidic phospholipid:water PC 
  Paplw <- (Kpu_bc - ((1 + Z)*dat_rbc$f_water + P*(dat_rbc$f_n_l+0.3*dat_rbc$f_n_pl))/(1+Y))*((1+Y)/(Z*dat_rbc$f_a_pl))

  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals

 
  # if(type_calc==1){  #moderate to strong bases
  #   Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + X*Paplw*dat_all$f_a_pl)/(1+W)
  #   Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + X*Paplw*dat_ad$f_a_pl)/(1+W)
  # }else if(type_calc==2){  #weak bases, acids and zwitterions (binds to albumin)
  #   Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + (1+X)*Pprwp_alb*(dat_all$fif_ap))/(1+W) #f_proteins:fraction volume of binding protein
  #   Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + (1+X)*Pprwp_alb*(dat_ad$fif_ap))/(1+W)
  # }else{  #neutrals (binds to lipoproteins)
  #   Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + (1+X)*Pprwp_lip*dat_all$fif_lp)/(1+W)
  #   Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + (1+X)*Pprwp_lip*dat_ad$fif_lp)/(1+W)
  # }
  # 
  
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + X*Paplw*dat_all$f_a_pl)/(1+W)
    Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + X*Paplw*dat_ad$f_a_pl)/(1+W)
  }else if(type_calc==2){  #weak bases, acids and zwitterions (binds to albumin)
    Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + (1+X)*Pprwp_alb*(dat_all$AR))/(1+W) #f_proteins:fraction volume of binding protein
    Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + (1+X)*Pprwp_alb*(dat_ad$AR))/(1+W)
  }else{  #neutrals (binds to lipoproteins)
    Kp_all <- ((1+X)*(dat_all$f_water + 0.7*dat_all$f_n_pl) + P*(dat_all$f_n_l + 0.3*dat_all$f_n_pl) + (1+X)*Pprwp_lip*dat_all$LR)/(1+W)
    Kp_ad <-  ((1+X)*(dat_ad$f_water + 0.7*dat_ad$f_n_pl) + P_OW*(dat_ad$f_n_l + 0.3*dat_ad$f_n_pl) + (1+X)*Pprwp_lip*dat_ad$LR)/(1+W)
  }
  
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms
  
  return(Kp)
  
}

