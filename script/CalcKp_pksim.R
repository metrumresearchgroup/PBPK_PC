#This function calculates the tissue:plasma partition coefficients according to PK-Sim https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159
.libPaths("lib")
library(dplyr)

calcKp_pksim <- function(logP, fup, dat){
  #logMA is the log of membrane affinity = phosphatidylcholin:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)
  #dat <- read.csv("../data/tissue_comp_Schmitt.csv")
  
  #dat <- read.csv("/data/internship-summer-2018/data/PKSim_tissue_comp_pksim.csv") #PK-Sim data
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark.csv")
  #dat <- read.csv("/data/internship-summer-2018/data/tissue_comp_Ruark_rat_for_R&R.csv") #Rat physiology from Ruark et al.
  
  #dat <- read.csv("/data/internship-summer-2018/data/unified_tissue_comp.csv") # Unified physiology from Ruark, P&T, R&R, and PK-Sim
  ##dat <- read.csv("../data/unified_tissue_comp.csv")
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for adipose and RBCs
  

  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  #K_protein <- 0.163 + 0.0221*K_n_pl    #protein:water partition; Schmitt, Walter (2008)
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5 # From PK-Sim (very similar value to the other method)
  
  
  kp <- (dat_all$f_water + (K_n_pl*dat_all$f_lipids) + (K_protein*dat_all$f_proteins))*fup

  #denom <- 0.945 + (10^logMA*0.00575) + (0.93*fup)  #plasma fractions
  #kp <- kp/denom  #according to Willmann et al. (2005)
  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name
  
  return(Kp)
}