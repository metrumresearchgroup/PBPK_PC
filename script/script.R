# This script reproduces manuscript figures
# Clear environment
rm(list=ls())

# Load libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(mrgsolve)
library(gridExtra)
library(PKNCA)
library(msm)
library(kableExtra)
source("sig.R")
source("CalcKp_P&T.R")
source("CalcKp_R&R.R")
source("CalcKpu_R&R.R")
source("CalcKp_Berez.R")
source("CalcKp_Schmitt.R")
source("CalcKp_pksim.R")
source("CalcVss_P&T.R")

source("getPartitionCoeff.R")

# Load physiological data
dat_uni <- read.csv("../data/unified_tissue_comp.csv")           # unified tissue composition
dat_PT <- read.csv("../data/tissue_comp_P&T.csv")                # data reported by Poulin and Theil
dat_Berez <- read.csv("../data/PKSim_tissue_comp_PT_Berez.csv")  # data used by PK-Sim for PT and Berez methods
dat_RR <- read.csv("../data/tissue_comp_R&R.csv")                # data reported by Rodgers and Rowland
dat_Schmitt_rep <- read.csv("../data/tissue_comp_Schmitt.csv")   # data reported by Schmitt
dat_Schmitt <- read.csv("../data/PKSim_tissue_comp_Schmitt.csv") # data used by PK-Sim for Schmitt method
dat_pksim <- read.csv("../data/PKSim_tissue_comp_pksim.csv")     # data used by PK-Sim for PK-Sim method

dat_PT_RR <- read.csv("../data/tissue_comp_R&R_for_P&T.csv")          # data reported by R&R formatted to use with PT method
dat_pksim_RR <- read.csv("../data/tissue_comp_R&R_for_pksim.csv")     # data reported by R&R formatted to use with PK-Sim method
dat_Schmitt_RR <- read.csv("../data/tissue_comp_R&R_for_Schmitt.csv") # data reported by R&R formatted to use with Schmitt method

# Set functions 
filter <- dplyr::filter
mutate <- dplyr::mutate
select <- dplyr::select

# Set figure default themes
# Theme 1 used for Fig. 2, 3 and 4
th1 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 10),
             legend.justification=c(0.025,0.95), legend.position=c(0.025,0.95), legend.key=element_blank(), 
             legend.title = element_blank(),
             plot.title = element_text(face="bold", size=15))

# Theme 5 used for Fig. 5
th5 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 10),
              legend.justification=c(1,1), legend.position=c(1,1), legend.key=element_blank(), 
              legend.title = element_blank(),
              plot.title = element_text(face="bold", size=15),
              axis.text.x  = element_text(angle=45, vjust=0.5))

# Theme 6 used for Fig. 6
th6 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10),
               legend.justification=c(1,1), legend.position=c(1,1.05), legend.key=element_blank(),
               plot.title = element_text(face="bold", size=15))

# Theme 7 used for Fig. 7
th7 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 10),
               legend.justification=c(1,1), legend.position=c(1,1), legend.key=element_blank(), legend.title = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(face="bold", size=15))

#########################################################################################################
######################################## CHUNK 1: Figure 2 ##############################################
#########################################################################################################
## This figure validates our code by comparing our Vss/Kp predictions (using literature tissue composition) 
## to the reported values in literature
## note: we will stick to the original reported values; eg: Vss for PT and Kpu for RR
## note: we will use as many of the reported drugs as possible for validation trying to include the drugs
## we investigated in simulation if their Kp/Vss are reported

# Define location of Pearson correlation coef. on every plot (fraction of maximum x and y value f plot)
pcc_x <- 0.75
pcc_y <- 0.35

## Poulin and Theil figure ##

### Acids ###
# Acetimenophen, Ascorbic acid, Isoniazid, Thiopental, Acetozolamid, Tetrahydrocam
pred_PT_a <- c(calcVss_PT(0.46, 9.38, 1, BP=1, type=2, dat_PT), 
               calcVss_PT(-1.85, 4.7, 0.76, BP=0.82, type=2, dat_PT), 
               calcVss_PT(-0.7, 1.82, 1, BP=1, type=2, dat_PT), 
               calcVss_PT(2.85, 7.45, 0.112, BP=1, type=2, dat_PT), 
               calcVss_PT(-0.26, 7.4, 0.044, BP=1, type=2, dat_PT), 
               calcVss_PT(6.97, 10.6, 0.115, BP=0.55, type=2, dat_PT))
rep_PT_a <- c(0.63,0.5,0.56,2.87, 0.32,5.53)  
acid <- data.frame(pred = unlist(pred_PT_a),rep_pred = rep_PT_a) %>%
  mutate(Drug_class = "Acid")

### Strong bases ###
# Caffeine, Acebutolol, Chlorpheniramine, Carvedilol, Mepivacaine, Trimethoprim
pred_PT_b <- c(calcVss_PT(-0.07, 10.4, 0.681, BP=0.98, type=3, dat_PT), 
               calcVss_PT(1.71, 9.67, 0.743, BP=1, type=3, dat_PT), 
               calcVss_PT(3.38, 9.13, 0.27, BP=1.34, type=3, dat_PT), 
               calcVss_PT(4.19, 8.4, 0.0054, BP=0.71, type=3, dat_PT), 
               calcVss_PT(1.95, 7.76, 0.225, BP=0.9, type=3, dat_PT), 
               calcVss_PT(0.91, 7.12, 0.304, BP=1, type=3, dat_PT))
rep_PT_b <- c(0.53,1.25,3.26,2.68,1.22,0.52)
sbase <- data.frame(pred = unlist(pred_PT_b),rep_pred = rep_PT_b) %>%
  mutate(Drug_class = "Strong base")

### Weak bases ###
# Alfentanil, Cimetidine, Etomidate
pred_PT_wb <- c(calcVss_PT(2.16, 6.5, 0.074, BP=1, type=3, dat_PT),
                calcVss_PT(0.4, 6.8, 0.8, BP=1, type=3, dat_PT),
                calcVss_PT(3.05, 4.5, 0.235, BP=1, type=3, dat_PT))
rep_PT_wb <- c(1.38,0.56,4.55)
wbase <- data.frame(pred = unlist(pred_PT_wb),rep_pred = rep_PT_wb) %>%
  mutate(Drug_class = "Weak base")

### Neutrals ###
# Coumarin, Griseofulvin, Cyclosporin
pred_PT_n <- c(calcVss_PT(1.39, 0, 0.853, BP=1, type=1, dat_PT),
               calcVss_PT(2.18, 0, 0.2, BP=1, type=1, dat_PT),
               calcVss_PT(2.92, 0, 0.0602, BP=1.36, type=1, dat_PT))
rep_PT_n <- c(1.08,1.72,2.77)
neut <- data.frame(pred = unlist(pred_PT_n),rep_pred = rep_PT_n) %>%
  mutate(Drug_class = "Neutral")

### Zwitterions ###
#Ciprofloxacin, Enoxacin 
pred_PT_z <-c(calcVss_PT(-0.2,c(6.3,8.6),0.485,1,type=6, dat_PT),
              calcVss_PT(0.28,c(6.3,8.3),0.52,1,type=6, dat_PT))
rep_PT_z <- c(0.47,0.44)
zwit <- data.frame(pred = unlist(pred_PT_z),rep_pred = rep_PT_z) %>%
  mutate(Drug_class = "Zwitterion")

drug_all_PT <- rbind(sbase,wbase,acid,neut,zwit)

# Calculate the Pearson correlation coefficient 
corr_coeff_PT <- round(cor(drug_all_PT$pred,drug_all_PT$rep_pred,method='pearson'), digits=2)
corr_coeff_PT <- paste("italic(r) %~~%", corr_coeff_PT)

fig2a <- ggplot() +
  geom_point(data=drug_all_PT, aes(x=pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_PT[1:2])) +
  ylim(0,max(drug_all_PT[1:2])) +
  xlab("Predicted Vss") +
  ylab("Reported Vss") +
  ggtitle("a    PT") + 
  th1 +
  geom_text(aes(x=max(drug_all_PT[1:2])*pcc_x, y=max(drug_all_PT[1:2])*pcc_y, label=corr_coeff_PT), parse=TRUE)


## Berezhkovskiy figure ##

### Strong base ###
# Acebutolol-R 
pred_ace <- calcKp_Berez(1.4,9.7,0.79,BP=1,type=3,dat_Berez)
pksim_pred_ace <- c(0.16,1.82,2.01,0.99,1.18,1.67,1.67,1.67,1.48,0.82,1.17,1.46,1.24,1.2)
ace <- data.frame(pred = unlist(pred_ace),rep_pred = pksim_pred_ace) %>%
  mutate(Drug_class = "Strong base")

### Weak base ###
# Alprazolam
pred_alp <- calcKp_Berez(2.5,2.4,0.35,BP=1,type=3,dat_Berez)
pksim_pred_alp <- c(7.75,6.06,5.72,1.66,2.39,4.56,4.56,4.56,3.7,0.84,2.41,3.76,2.84,2.43)
alp <- data.frame(pred = unlist(pred_alp),rep_pred = pksim_pred_alp) %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Phenobarbital 
pred_phen <- calcKp_Berez(0.8,7.35,0.78,BP=0.861,type=2,dat_Berez)
pksim_pred_phen <- c(0.26,0.78,1.09,0.79,0.86,0.94,0.94,0.94,0.92,0.78,0.84,0.88,0.83,0.87)
phen <- data.frame(pred = unlist(pred_phen),rep_pred = pksim_pred_phen) %>%
  mutate(Drug_class = "Acid")

### Neutral ###
pred_dig <- calcKp_Berez(1.48,0,0.87,BP=1,type=1,dat_Berez)
pksim_pred_dig <- c(1.61,2.24,2.4,1.09,1.34,1.98,1.98,1.98,1.73,0.87,1.32,1.72,1.43,1.36)
dig <- data.frame(pred = unlist(pred_dig),rep_pred = pksim_pred_dig) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion
# Ofloxacin
pred_oflo <- calcKp_Berez(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6,dat_Berez)
pksim_pred_oflo <- c(0.15,0.43,0.77,0.72,0.75,0.7,0.7,0.7,0.73,0.76,0.72,0.68,0.69,0.76)
oflo <- data.frame(pred = unlist(pred_oflo),rep_pred = pksim_pred_oflo) %>%
  mutate(Drug_class = "Zwitterion")

drug_all_Berez <- rbind(ace,alp,phen,dig,oflo)

# Calculate the Pearson correlation coefficient 
corr_coeff_Berez <- round(cor(drug_all_Berez$pred,drug_all_Berez$rep_pred,method='pearson'), digits=2)
corr_coeff_Berez <- paste("italic(r) %~~%", corr_coeff_Berez)

fig2b <- ggplot() +
  geom_point(data=drug_all_Berez, aes(x=pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_Berez[1:2])) +
  ylim(0,max(drug_all_Berez[1:2])) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("b    Berez") + 
  th1 +
  geom_text(aes(x=max(drug_all_Berez[1:2])*pcc_x, y=max(drug_all_Berez[1:2])*pcc_y, label=corr_coeff_Berez), parse=TRUE)


## Rodgers and Rowland figure ##

### Acids ###
# Cefazolin
pred_RR_a_all <- calcKpu_RR(0.3, 2.3, 0.16, BP=0.55, type=2, dat_RR)
pred_RR_a <- c(pred_RR_a_all[1], pred_RR_a_all[3:8], pred_RR_a_all[10])
rep_RR_a <- c(0.74,1.28,1.23,1.1,0.79,1.52,0.69,1.92)
cef <- data.frame(pred = unlist(pred_RR_a),rep = rep_RR_a) %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
pred_RR_a_all2 <- calcKpu_RR(1.6, 7.4, 0.64, BP=0.861, type=2, dat_RR)
pred_RR_a2 <- c(pred_RR_a_all2[3:8], pred_RR_a_all2[10], pred_RR_a_all2[13])
rep_RR_a2 <- c(1.43,0.88,0.97,0.92,1.07,0.79,1.84,1.23)
phen <- data.frame(pred = unlist(pred_RR_a2),rep = rep_RR_a2) %>%
  mutate(Drug_class = "Acid")


### Weak bases ###
# Alfentanil
pred_RR_b_all <- calcKpu_RR(2.2, 6.5, 0.11, BP=1, type=3, dat_RR)
pred_RR_b <- c(pred_RR_b_all[2:11], pred_RR_b_all[13])
rep_RR_b <- c(6.3,7.33,3.97,4.24,4.05,5.48,2.88,6.71,10.7,2.7,9.2)
alf <- data.frame(pred = unlist(pred_RR_b),rep = rep_RR_b) %>%
  mutate(Drug_class = "Weak base")

### Strong bases ###
# Bisoprolol-R
pred_RR_sb_all <- calcKpu_RR(1.87, 9.4, 0.85, BP=1.36, type=3, dat_RR)
pred_RR_sb <- c(pred_RR_sb_all[1:8], pred_RR_sb_all[10], pred_RR_sb_all[12:13])
rep_RR_sb <- c(3.61,3.36,11,10.4,21.4,19.6,17.2,7.73,6.37,10.8,1.79)
bis <- data.frame(pred = unlist(pred_RR_sb),rep = rep_RR_sb) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
pred_RR_sb_all2 <- calcKpu_RR(1.87, 9.7, 0.79, BP=1.09, type=3, dat_RR)
pred_RR_sb2 <- c(pred_RR_sb_all2[1:8], pred_RR_sb_all2[10], pred_RR_sb_all2[12:13])
rep_RR_sb2 <- c(2.34,2.57,6.45,6.1,11.8,11.0,9.67,4.84,3.85,6.45,0.82)
ace <- data.frame(pred = unlist(pred_RR_sb2),rep = rep_RR_sb2) %>%
  mutate(Drug_class = "Strong base")

### Neutral ###
# Ethoxybenzamide
pred_RR_n_all <- calcKpu_RR(0.77, 0, 0.5, BP=1, type=1, dat_RR)
pred_RR_n <- c(pred_RR_n_all[2:5], pred_RR_n_all[7:8], pred_RR_n_all[10:11], pred_RR_n_all[13])  
rep_RR_n <- c(1.03,1.1,0.79,0.91,0.86,0.85,1.11,0.78,0.48)
etho <- data.frame(pred = unlist(pred_RR_n),rep = rep_RR_n) %>%
  mutate(Drug_class = "Neutral")

### Zwitterions ###
# Enoxacin
pred_eno <- calcKpu_RR(log10(1.26), c(6.1,8.7), 0.66, BP=0.94, type=6, dat_RR)
pred_eno <- c(pred_eno[1], pred_eno[4:8], pred_eno[10:11])
rep_eno <- c(1.74,4.62,9.88,9.04,7.61,3.74,3.16,6.33)
eno <- data.frame(pred = unlist(pred_eno),rep = rep_eno) %>%
  mutate(Drug_class = "Zwitterion")


drug_all_RR <- rbind(alf,ace,bis,phen,cef,etho,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_RR <- round(cor(drug_all_RR$pred,drug_all_RR$rep,method='pearson'), digits=2)
corr_coeff_RR <- paste("italic(r) %~~%", corr_coeff_RR)

fig2c <- ggplot() +
  geom_point(data=drug_all_RR, aes(x=pred, y=rep, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_RR[1:2])) +
  ylim(0,max(drug_all_RR[1:2])) +
  xlab("Predicted Kpu") +
  ylab("Reported Kpu") +
  ggtitle("c    RR") + 
  th1 +
  geom_text(aes(x=max(drug_all_RR[1:2])*pcc_x, y=max(drug_all_RR[1:2])*pcc_y, label=corr_coeff_RR), parse=TRUE)


## Schmitt figure ##
# Using logMA below, could also use logP

### Strong base ###
# Quinidine 
pred_quin <- calcKp_Schmitt(2.3, 8.6, 0.33, type=3, dat_Schmitt)
pksim_pred_quin <- c(1.6,22.23,2.85,16.4,12.89,6.22,6.22,6.22,14.85,7.17,1.99,3.26,5.32,5.17)
quin <- data.frame(pred = unlist(pred_quin),pksim_pred = pksim_pred_quin) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
pred_ace <- calcKp_Schmitt(1.4, 9.7, 0.79, type=3, dat_Schmitt)
pksim_pred_ace <- c(0.69,7.38,0.38,5.52,4.51,2.45,2.45,2.45,5.08,2.69,1.15,1.48,1.94,2.13)
ace <- data.frame(pred = unlist(pred_ace),pksim_pred = pksim_pred_ace) %>%
  mutate(Drug_class = "Strong base")

### Weak base ###
# Alprazolam
pred_alp <- calcKp_Schmitt(2.5, 2.4, 0.35, type=3, dat_Schmitt)
pksim_pred_alp <- c(2.83,12.65,102.27,12.9,7.42,8.4,8.4,8.4,9.53,4.96,1.85,6.21,16.68,3.03)
alp <- data.frame(pred = unlist(pred_alp),pksim_pred = pksim_pred_alp) %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Phenobarbital
pred_phen <- calcKp_Schmitt(0.8, 7.35, 0.78, type=2, dat_Schmitt)
pksim_pred_phen <- c(0.32,1.08,2.95,1.02,0.87,0.85,0.85,0.85,0.9,0.78,0.68,0.73,0.96,0.72)
phen <- data.frame(pred = unlist(pred_phen),pksim_pred = pksim_pred_phen) %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
pred_dig <- calcKp_Schmitt(1.48, 0, 0.87, type=1, dat_Schmitt)
pksim_pred_dig <- c(0.9,3.64,24.31,3.64,2.36,2.63,2.63,2.63,2.82,1.77,1.06,2.05,4.38,1.34)
dig <- data.frame(pred = unlist(pred_dig),pksim_pred = pksim_pred_dig) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
pred_oflo <- calcKp_Schmitt(-0.4, c(5.97,9.28), 0.77, type=6, dat_Schmitt)
pksim_pred_oflo <- c(0.23,0.66,0.29,0.6,0.61,0.64,0.64,0.64,0.58,0.6,0.61,0.58,0.46,0.62)
oflo <- data.frame(pred = unlist(pred_oflo),pksim_pred = pksim_pred_oflo) %>%
  mutate(Drug_class = "Zwitterion")

drug_all_Schmitt <- rbind(quin,ace,alp,phen,dig,oflo)

# Calculate the Pearson correlation coefficient 
corr_coeff_Schmitt <- round(cor(drug_all_Schmitt$pred,drug_all_Schmitt$pksim_pred,method='pearson'), digits=2)
corr_coeff_Schmitt <- paste("italic(r) %~~%", corr_coeff_Schmitt)

fig2d <- ggplot() +
  geom_point(data=drug_all_Schmitt, aes(x=pred, y=pksim_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_Schmitt[1:2])) +
  ylim(0,max(drug_all_Schmitt[1:2])) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("d    Schmitt") + 
  th1 +
  geom_text(aes(x=max(drug_all_Schmitt[1:2])*pcc_x, y=max(drug_all_Schmitt[1:2])*pcc_y, label=corr_coeff_Schmitt), parse=TRUE)


## PK-Sim standard figure ##

### Strong Base ###
# Acebutolol-R
pred_ace <- calcKp_pksim(1.87, 0.79, dat_pksim)
pksim_pred_ace <- c(16.44,7.2,47.04,6.67,3.9,4.45,4.45,4.45,4.89,1.48,1.65,5.43,6.75,1.83)
ace <- data.frame(pred = unlist(pred_ace),pksim_pred = pksim_pred_ace) %>%
  mutate(Drug_class = "Strong base")

### Weak base ###
# Triazolam
pred_tri <- calcKp_pksim(2.4, 0.28, dat_pksim)
pksim_pred_tri <- c(19.41,8.09,56.39,7.51,4.15,4.79,4.79,4.79,5.36,1.22,1.42,6.03,7.66,1.65)
tri <- data.frame(pred = unlist(pred_tri),pksim_pred = pksim_pred_tri) %>%
  mutate(Drug_class = "Weak base")

### Acids ###
# Phenobarbital
pred_phen <- calcKp_pksim(1.6, 0.64, dat_pksim)
pksim_pred_phen <- c(7.3,3.37,20.51,3.13,1.93,2.18,2.18,2.18,2.36,0.89,0.97,2.58,3.13,1.03)
phen <- data.frame(pred = unlist(pred_phen),pksim_pred = pksim_pred_phen) %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Ethoxybenzamide
pred_etho <- calcKp_pksim(0.77, 0.5, dat_pksim)
pksim_pred_etho <- c(1.06,0.74,2.44,0.68,0.57,0.6,0.6,0.6,0.6,0.46,0.47,0.61,0.64,0.46)
etho <- data.frame(pred = unlist(pred_etho),pksim_pred = pksim_pred_etho) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Enoxacin
pred_eno <- calcKp_pksim(0.100, 0.66, dat_pksim)
pksim_pred_eno <- c(0.56,0.63,0.77,0.59,0.58,0.59,0.59,0.59,0.57,0.56,0.57,0.55,0.52,0.55)
eno <- data.frame(pred = unlist(pred_eno),pksim_pred = pksim_pred_eno) %>%
  mutate(Drug_class = "Zwitterion")

drug_all_pksim <- rbind(ace,alp,phen,dig,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_pksim <- round(cor(drug_all_pksim$pred,drug_all_pksim$pksim_pred,method='pearson'), digits=2)
corr_coeff_pksim <- paste("italic(r) %~~%", corr_coeff_pksim)

fig2e <- ggplot() +
  geom_point(data=drug_all_pksim, aes(x=pred, y=pksim_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_pksim[1:2])) +
  ylim(0,max(drug_all_pksim[1:2])) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("e    PK-Sim") + 
  th1 +
  geom_text(aes(x=max(drug_all_pksim[1:2])*pcc_x, y=max(drug_all_pksim[1:2])*pcc_y, label=corr_coeff_pksim), parse=TRUE)


fig2 <- grid.arrange(fig2a, fig2b, fig2c, fig2d, fig2e, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig2.png", fig2, width=8, height=12)

#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 2: Figure 3 ##############################################
#########################################################################################################
## This figure swaps physiologies between calculation methods to see impact on Kp/Vss predictions

# Define location of Pearson correlation coef. on every plot (fraction of maximum x and y value f plot)
pcc_x <- 0.75
pcc_y <- 0.35


# Use Rodgers and Rowland for Poulin and Theil method 
# Sum extracellular and intracellular water; sum neutral and acidic phospholipids

### Strong Base ###
# Metoprolol
met_PT <- calcKp_PT(2.15, 9.7, 0.8, BP=1.52, type=3, dat_PT) 
met_RR <- calcKp_PT(2.15, 9.7, 0.8, BP=1.52, type=3, dat_PT_RR)
met <- data.frame(rep_pred = unlist(met_PT), swap_pred = unlist(met_RR)) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace_PT <- calcKp_PT(1.4, 9.7, 0.79, BP=1, type=3, dat_PT) 
ace_RR <- calcKp_PT(1.4, 9.7, 0.79, BP=1, type=3, dat_PT_RR)
ace <- data.frame(rep_pred = unlist(ace_PT), swap_pred = unlist(ace_RR)) %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori_PT <- calcKp_PT(2.56, 1.76, 0.42, BP=1, type=3, dat_PT) 
vori_RR <- calcKp_PT(2.56, 1.76, 0.42, BP=1, type=3, dat_PT_RR)
vori <- data.frame(rep_pred = unlist(vori_PT), swap_pred = unlist(vori_RR)) %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp_PT <- calcKp_PT(2.5, 2.4, 0.35, BP=1, type=3, dat_PT) 
alp_RR <- calcKp_PT(2.5, 2.4, 0.35, BP=1, type=3, dat_PT_RR)
alp <- data.frame(rep_pred = unlist(alp_PT), swap_pred = unlist(alp_RR)) %>%
  mutate(Drug_class = "Weak base")


### Acid ###
# Thiopental
thio_PT <- calcKp_PT(2.9, 7.5, 0.13, BP=1, type=2, dat_PT) 
thio_RR <- calcKp_PT(2.9, 7.5, 0.13, BP=1, type=2, dat_PT_RR)
thio <- data.frame(rep_pred = unlist(thio_PT), swap_pred = unlist(thio_RR)) %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen_PT <- calcKp_PT(0.8, 7.35, 0.78, BP=1, type=2, dat_PT) 
phen_RR <- calcKp_PT(0.8, 7.35, 0.78, BP=1, type=2, dat_PT_RR)
phen <- data.frame(rep_pred = unlist(phen_PT), swap_pred = unlist(phen_RR)) %>%
  mutate(Drug_class = "Acid")


### Neutral ###
# Digoxin
dig_PT <- calcKp_PT(1.48, 0, 0.87, BP=1, type=1, dat_PT) 
dig_RR <- calcKp_PT(1.48, 0, 0.87, BP=1, type=1, dat_PT_RR)
dig <- data.frame(rep_pred = unlist(dig_PT), swap_pred = unlist(dig_RR)) %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho_PT <- calcKp_PT(0.77, 0, 0.5, BP=1, type=1, dat_PT) 
etho_RR <- calcKp_PT(0.77, 0, 0.5, BP=1, type=1, dat_PT_RR)
etho <- data.frame(rep_pred = unlist(etho_PT), swap_pred = unlist(etho_RR)) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo_PT <- calcKp_PT(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_PT) 
oflo_RR <- calcKp_PT(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_PT_RR)
oflo <- data.frame(rep_pred = unlist(oflo_PT), swap_pred = unlist(oflo_RR)) %>%
  mutate(Drug_class = "Zwitterion")
# Enoxacin
eno_PT <- calcKp_PT(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_PT) 
eno_RR <- calcKp_PT(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_PT_RR)
eno <- data.frame(rep_pred = unlist(eno_PT), swap_pred = unlist(eno_RR)) %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_PT <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient
corr_coeff_PT <- round(cor(drug_all_PT$rep_pred,drug_all_PT$swap_pred,method='pearson'), digits=2)
corr_coeff_PT <- paste("italic(r) %~~%", corr_coeff_PT)

fig3a <- ggplot() +
  geom_point(data=drug_all_PT, aes(x=swap_pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_PT[1:2])) +
  ylim(0,max(drug_all_PT[1:2])) +
  xlab("Predicted Kp (RR tissue comp.)") +
  ylab("Predicted Kp (reported tissue comp.)") +
  ggtitle("a    PT with RR tissue comp.") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_PT[1:2])*pcc_x, y=max(drug_all_PT[1:2])*pcc_y, label=corr_coeff_PT), parse=TRUE)


# Use R&R for Berez. (sum extracellular and intracellular water; sum neutral and acidic phospholipids)
### Strong Base ###
# Metoprolol
met_PT <- calcKp_Berez(2.15, 9.7, 0.8, BP=1.52, type=3, dat_PT) 
met_RR <- calcKp_Berez(2.15, 9.7, 0.8, BP=1.52, type=3, dat_PT_RR)
met <- data.frame(rep_pred = unlist(met_PT), swap_pred = unlist(met_RR)) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace_PT <- calcKp_Berez(1.4, 9.7, 0.79, BP=1, type=3, dat_PT) 
ace_RR <- calcKp_Berez(1.4, 9.7, 0.79, BP=1, type=3, dat_PT_RR)
ace <- data.frame(rep_pred = unlist(ace_PT), swap_pred = unlist(ace_RR)) %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori_PT <- calcKp_Berez(2.56, 1.76, 0.42, BP=1, type=3, dat_PT) 
vori_RR <- calcKp_Berez(2.56, 1.76, 0.42, BP=1, type=3, dat_PT_RR)
vori <- data.frame(rep_pred = unlist(vori_PT), swap_pred = unlist(vori_RR)) %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp_PT <- calcKp_Berez(2.5, 2.4, 0.35, BP=1, type=3, dat_PT) 
alp_RR <- calcKp_Berez(2.5, 2.4, 0.35, BP=1, type=3, dat_PT_RR)
alp <- data.frame(rep_pred = unlist(alp_PT), swap_pred = unlist(alp_RR)) %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio_PT <- calcKp_Berez(2.9, 7.5, 0.13, BP=1, type=2, dat_PT) 
thio_RR <- calcKp_Berez(2.9, 7.5, 0.13, BP=1, type=2, dat_PT_RR)
thio <- data.frame(rep_pred = unlist(thio_PT), swap_pred = unlist(thio_RR)) %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen_PT <- calcKp_Berez(0.8, 7.35, 0.78, BP=1, type=2, dat_PT) 
phen_RR <- calcKp_Berez(0.8, 7.35, 0.78, BP=1, type=2, dat_PT_RR)
phen <- data.frame(rep_pred = unlist(phen_PT), swap_pred = unlist(phen_RR)) %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig_PT <- calcKp_Berez(1.48, 0, 0.87, BP=1, type=1, dat_PT) 
dig_RR <- calcKp_Berez(1.48, 0, 0.87, BP=1, type=1, dat_PT_RR)
dig <- data.frame(rep_pred = unlist(dig_PT), swap_pred = unlist(dig_RR)) %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho_PT <- calcKp_Berez(0.77, 0, 0.5, BP=1, type=1, dat_PT) 
etho_RR <- calcKp_Berez(0.77, 0, 0.5, BP=1, type=1, dat_PT_RR)
etho <- data.frame(rep_pred = unlist(etho_PT), swap_pred = unlist(etho_RR)) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo_PT <- calcKp_Berez(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_PT) 
oflo_RR <- calcKp_Berez(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_PT_RR)
oflo <- data.frame(rep_pred = unlist(oflo_PT), swap_pred = unlist(oflo_RR)) %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno_PT <- calcKp_Berez(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_PT) 
eno_RR <- calcKp_Berez(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_PT_RR)
eno <- data.frame(rep_pred = unlist(eno_PT), swap_pred = unlist(eno_RR)) %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_Berez <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient
corr_coeff_Berez <- round(cor(drug_all_Berez$rep_pred,drug_all_Berez$swap_pred,method='pearson'), digits=2)
corr_coeff_Berez <- paste("italic(r) %~~%", corr_coeff_Berez)

fig3b <- ggplot() +
  geom_point(data=drug_all_Berez, aes(x=swap_pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_Berez[1:2])) +
  ylim(0,max(drug_all_Berez[1:2])) +
  xlab("Predicted Kp (RR tissue comp.)") +
  ylab("Predicted Kp (reported tissue comp.)") +
  ggtitle("b    Berez with RR tissue comp.") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_Berez[1:2])*pcc_x, y=max(drug_all_Berez[1:2])*pcc_y, label=corr_coeff_Berez), parse=TRUE)


# Use R&R for Schmitt (sum extracellular and intracellular water) 
### Strong Base ###
# Metoprolol
met_Schmitt <- calcKp_Schmitt(2.15, 9.7, 0.8, type=3, dat_Schmitt_rep) 
met_RR <- calcKp_Schmitt(2.15, 9.7, 0.8, type=3, dat_Schmitt_RR)
met <- data.frame(rep_pred = unlist(met_Schmitt), swap_pred = unlist(met_RR)) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace_Schmitt <- calcKp_Schmitt(1.4, 9.7, 0.79, type=3, dat_Schmitt_rep) 
ace_RR <- calcKp_Schmitt(1.4, 9.7, 0.79, type=3, dat_Schmitt_RR)
ace <- data.frame(rep_pred = unlist(ace_Schmitt), swap_pred = unlist(ace_RR)) %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori_Schmitt <- calcKp_Schmitt(2.56, 1.76, 0.42, type=3, dat_Schmitt_rep) 
vori_RR <- calcKp_Schmitt(2.56, 1.76, 0.42, type=3, dat_Schmitt_RR)
vori <- data.frame(rep_pred = unlist(vori_Schmitt), swap_pred = unlist(vori_RR)) %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp_Schmitt <- calcKp_Schmitt(2.5, 2.4, 0.35, type=3, dat_Schmitt_rep) 
alp_RR <- calcKp_Schmitt(2.5, 2.4, 0.35, type=3, dat_Schmitt_RR)
alp <- data.frame(rep_pred = unlist(alp_Schmitt), swap_pred = unlist(alp_RR)) %>%
  mutate(Drug_class = "Weak base")


### Acid ###
# Thiopental
thio_Schmitt <- calcKp_Schmitt(2.9, 7.5, 0.13, type=2, dat_Schmitt_rep) 
thio_RR <- calcKp_Schmitt(2.9, 7.5, 0.13, type=2, dat_Schmitt_RR)
thio <- data.frame(rep_pred = unlist(thio_Schmitt), swap_pred = unlist(thio_RR)) %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen_Schmitt <- calcKp_Schmitt(0.8, 7.35, 0.78, type=2, dat_Schmitt_rep) 
phen_RR <- calcKp_Schmitt(0.8, 7.35, 0.78, type=2, dat_Schmitt_RR)
phen <- data.frame(rep_pred = unlist(phen_Schmitt), swap_pred = unlist(phen_RR)) %>%
  mutate(Drug_class = "Acid")


### Neutral ###
# Digoxin
dig_Schmitt <- calcKp_Schmitt(1.48, 0, 0.87, type=1, dat_Schmitt_rep) 
dig_RR <- calcKp_Schmitt(1.48, 0, 0.87, type=1, dat_Schmitt_RR)
dig <- data.frame(rep_pred = unlist(dig_Schmitt), swap_pred = unlist(dig_RR)) %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho_Schmitt <- calcKp_Schmitt(0.77, 0, 0.5, type=1, dat_Schmitt_rep) 
etho_RR <- calcKp_Schmitt(0.77, 0, 0.5, type=1, dat_Schmitt_RR)
etho <- data.frame(rep_pred = unlist(etho_Schmitt), swap_pred = unlist(etho_RR)) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo_Schmitt <- calcKp_Schmitt(-0.4, c(5.97,9.28), 0.77, type=6, dat_Schmitt_rep) 
oflo_RR <- calcKp_Schmitt(-0.4, c(5.97,9.28), 0.77, type=6, dat_Schmitt_RR)
oflo <- data.frame(rep_pred = unlist(oflo_Schmitt), swap_pred = unlist(oflo_RR)) %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno_Schmitt <- calcKp_Schmitt(0.1, c(6.1,8.7), 0.66, type=6, dat_Schmitt_rep) 
eno_RR <- calcKp_Schmitt(0.1, c(6.1,8.7), 0.66, type=6, dat_Schmitt_RR)
eno <- data.frame(rep_pred = unlist(eno_Schmitt), swap_pred = unlist(eno_RR)) %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_Schmitt <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient
corr_coeff_Schmitt <- round(cor(drug_all_Schmitt$rep_pred,drug_all_Schmitt$swap_pred,method='pearson'), digits=2)
corr_coeff_Schmitt <- paste("italic(r) %~~%", corr_coeff_Schmitt)

fig3c <- ggplot() +
  geom_point(data=drug_all_Schmitt, aes(x=swap_pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_Schmitt[1:2])) +
  ylim(0,max(drug_all_Schmitt[1:2])) +
  xlab("Predicted Kp (RR tissue comp.)") +
  ylab("Predicted Kp (reported tissue comp.)") +
  ggtitle("c    Schmitt with RR tissue comp.") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_Schmitt[1:2])*pcc_x, y=max(drug_all_Schmitt[1:2])*pcc_y, label=corr_coeff_Schmitt), parse=TRUE)


# Use R&R for PK-Sim standard (sum extracellular and intracellular water for water; replace lipids with neutral 
# lipids; sum albumin and lipoprotein for proteins)
### Strong Base ###
# Metoprolol
met_pksim <- calcKp_pksim(2.15, 0.8, dat_pksim) 
met_RR <- calcKp_pksim(2.15, 0.8, dat_pksim_RR)
met <- data.frame(rep_pred = unlist(met_pksim), swap_pred = unlist(met_RR)) %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace_pksim <- calcKp_pksim(1.4, 0.79, dat_pksim)
ace_RR <- calcKp_pksim(1.4, 0.79, dat_pksim_RR)
ace <- data.frame(rep_pred = unlist(ace_pksim), swap_pred = unlist(ace_RR)) %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori_pksim <- calcKp_pksim(2.56, 0.42, dat_pksim)
vori_RR <- calcKp_pksim(2.56, 0.42, dat_pksim_RR)
vori <- data.frame(rep_pred = unlist(vori_pksim), swap_pred = unlist(vori_RR)) %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp_pksim <- calcKp_pksim(2.5, 0.35, dat_pksim)
alp_RR <- calcKp_pksim(2.5, 0.35, dat_pksim_RR)
alp <- data.frame(rep_pred = unlist(alp_pksim), swap_pred = unlist(alp_RR)) %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio_pksim <- calcKp_pksim(2.9, 0.13, dat_pksim)
thio_RR <- calcKp_pksim(2.9, 0.13, dat_pksim_RR)
thio <- data.frame(rep_pred = unlist(thio_pksim), swap_pred = unlist(thio_RR)) %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen_pksim <- calcKp_pksim(0.8, 0.78, dat_pksim)
phen_RR <- calcKp_pksim(0.8, 0.78, dat_pksim_RR)
phen <- data.frame(rep_pred = unlist(phen_pksim), swap_pred = unlist(phen_RR)) %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig_pksim <- calcKp_pksim(1.48, 0.87, dat_pksim)
dig_RR <- calcKp_pksim(1.48, 0.87, dat_pksim_RR)
dig <- data.frame(rep_pred = unlist(dig_pksim), swap_pred = unlist(dig_RR)) %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho_pksim <- calcKp_pksim(0.77, 0.5, dat_pksim)
etho_RR <- calcKp_pksim(0.77, 0.5, dat_pksim_RR)
etho <- data.frame(rep_pred = unlist(etho_pksim), swap_pred = unlist(etho_RR)) %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo_pksim <- calcKp_pksim(-0.4, 0.77, dat_pksim)
oflo_RR <- calcKp_pksim(-0.4, 0.77, dat_pksim_RR)
oflo <- data.frame(rep_pred = unlist(oflo_pksim), swap_pred = unlist(oflo_RR)) %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno_pksim <- calcKp_pksim(0.1, 0.66, dat_pksim)
eno_RR <- calcKp_pksim(0.1, 0.66, dat_pksim_RR)
eno <- data.frame(rep_pred = unlist(eno_pksim), swap_pred = unlist(eno_RR)) %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_pksim <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient
corr_coeff_pksim <- round(cor(drug_all_pksim$rep_pred,drug_all_pksim$swap_pred,method='pearson'), digits=2)
corr_coeff_pksim <- paste("italic(r) %~~%", corr_coeff_pksim)

fig3d <- ggplot() +
  geom_point(data=drug_all_pksim, aes(x=swap_pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_pksim[1:2])) +
  ylim(0,max(drug_all_pksim[1:2])) +
  xlab("Predicted Kp (RR tissue comp.)") +
  ylab("Predicted Kp (reported tissue comp.") +
  ggtitle("d    PK-Sim with RR tissue comp.") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_pksim[1:2])*pcc_x, y=max(drug_all_pksim[1:2])*pcc_y, label=corr_coeff_pksim), parse=TRUE)


fig3 <- grid.arrange(fig3a, fig3b, fig3c, fig3d, ncol=2, nrow=2)
#ggsave(file="../deliv/figure/fig3.png", fig3, width=8, height=8)

#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 3: Figure 4 ##############################################
#########################################################################################################
## This figure compares Kp predictions using the unified tissue composition vs the reported tissue composition

# General function to calculate Kp values for tissues in both unified and reported physiologies
drug_pred <- function(logP, pKa, fup, BP, type, dat_uni, dat_method, method){
  if(method=="P&T"){
    uni <- calcKp_PT(logP, pKa, fup, BP, type, dat_uni)        #prediction using unified tissue composition
    pred <- calcKp_PT(logP, pKa, fup, BP, type, dat_method)    #prediction using reported tissue composition
    pred <- c(pred[1:3], pred[5:6], pred[4], pred[7:9], pred[11], pred[10])
  }
  else if(method=="Berezhkovskiy"){
    uni <- calcKp_Berez(logP, pKa, fup, BP, type, dat_uni)     #prediction using unified tissue composition
    pred <- calcKp_Berez(logP, pKa, fup, BP, type, dat_method) #prediction using reported tissue composition
    pred <- c(pred[1:3], pred[5:6], pred[4], pred[7:9], pred[11], pred[10])
  }
  else if(method=="R&R"){
    uni <- calcKp_RR(logP, pKa, fup, BP, type, dat_uni)        #prediction using unified tissue composition
    pred <- calcKp_RR(logP, pKa, fup, BP, type, dat_method)    #prediction using reported tissue composition
    pred <- c(pred[1:3],pred[5:6],pred[4],pred[7:9],pred[11:12])
  }
  else if(method=="Schmitt"){
    uni <- calcKp_Schmitt(logP, pKa, fup, type, dat_uni)       #prediction using unified tissue composition
    pred <- calcKp_Schmitt(logP, pKa, fup, type, dat_method)   #prediction using PK-Sim reported tissue composition
    pred <- c(pred[3], pred[1:2], pred[4:6], pred[9:11], pred[13:14])
  }
  else if(method=="Schmitt2"){
    uni <- calcKp_Schmitt(logP, pKa, fup, type, dat_uni)       #prediction using unified tissue composition
    pred <- calcKp_Schmitt(logP, pKa, fup, type, dat_method)   #prediction using reported tissue composition
    pred <- c(pred[1:3], pred[5:6], pred[4],pred[7:11])
  }
  else if(method=="PK-Sim"){
    uni <- calcKp_pksim(logP, fup, dat_uni)                    #prediction using unified tissue composition
    pred <- calcKp_pksim(logP, fup, dat_method)                #prediction using reported tissue composition
    pred <- c(pred[3], pred[1:2], pred[4:6], pred[9:11], pred[13:14])
  }
  drug_vals <- data.frame(uni_pred = unlist(uni),pred = unlist(pred))
  return(drug_vals)
}

# Define location of Pearson correlation coef. on every plot (fraction of maximum x and y value f plot)
pcc_x <- 0.75
pcc_y <- 0.35

# Poulin and Theil method
# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: adipose, bone, brain, gut, heart, kidney, liver, lung, muscle, spleen, skin

### Strong Base ###
#Metoprolol
met <- drug_pred(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace <- drug_pred(1.4, 9.7, 0.79, BP=1, type=3, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori <- drug_pred(2.56, 1.76, 0.42, BP=1, type=3, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp <- drug_pred(2.5, 2.4, 0.35, BP=1, type=3, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio <- drug_pred(2.9, 7.5, 0.13, BP=1, type=2, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen <- drug_pred(0.8, 7.35, 0.78, BP=1, type=2, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig <- drug_pred(1.48, 0, 0.87, BP=1, type=1, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho <- drug_pred(0.77, 0, 0.5, BP=1, type=1, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo <- drug_pred(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno <- drug_pred(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_uni, dat_PT, "P&T") %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_PT <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_PT <- round(cor(drug_all_PT$uni_pred,drug_all_PT$pred,method='pearson'), digits=2)
corr_coeff_PT <- paste("italic(r) %~~%", corr_coeff_PT)

fig4a <- ggplot() +
  geom_point(data=drug_all_PT, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_PT[1:2])) +
  ylim(0,max(drug_all_PT[1:2])) +
  xlab("Predicted Kp (standardized tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  ggtitle("a    PT") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_PT[1:2])*pcc_x, y=max(drug_all_PT[1:2])*pcc_y, label=corr_coeff_PT), parse=TRUE)


# Berezkovskiy method
# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: adipose, bone, brain, gut, heart, kidney, liver, lung, muscle, spleen, skin

# We use dat_PT as the reported tissue composition because Berezkovskiy uses the same tissue composition as Poulin and Theil

### Strong Base ###
#Metoprolol
met <- drug_pred(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace <- drug_pred(1.4, 9.7, 0.79, BP=1, type=3, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori <- drug_pred(2.56, 1.76, 0.42, BP=1, type=3, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp <- drug_pred(2.5, 2.4, 0.35, BP=1, type=3, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio <- drug_pred(2.9, 7.5, 0.13, BP=1, type=2, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen <- drug_pred(0.8, 7.35, 0.78, BP=1, type=2, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig <- drug_pred(1.48, 0, 0.87, BP=1, type=1, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho <- drug_pred(0.77, 0, 0.5, BP=1, type=1, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo <- drug_pred(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno <- drug_pred(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_uni, dat_PT, "Berezhkovskiy") %>%
  mutate(Drug_class = "Zwitterion")


# Bind all the outputs together
drug_all_Berez <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_Berez <- round(cor(drug_all_Berez$uni_pred,drug_all_Berez$pred,method='pearson'),digits=2)
corr_coeff_Berez <- paste("italic(r) %~~%", corr_coeff_Berez)

fig4b <- ggplot() +
  geom_point(data=drug_all_Berez, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_Berez[1:2])) +
  ylim(0,max(drug_all_Berez[1:2])) +
  xlab("Predicted Kp (standardized tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  ggtitle("b    Berez") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_Berez[1:2])*pcc_x, y=max(drug_all_Berez[1:2])*pcc_y, label=corr_coeff_Berez), parse=TRUE)


# Rodgers and Rowland method
# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: adipose, bone, brain, gut, heart, kidney liver, lung, muscle, pancreas, skin, spleen, thymus

### Strong Base ###
#Metoprolol
met <- drug_pred(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace <- drug_pred(1.4, 9.7, 0.79, BP=1, type=3, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori <- drug_pred(2.56, 1.76, 0.42, BP=1, type=3, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp <- drug_pred(2.5, 2.4, 0.35, BP=1, type=3, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio <- drug_pred(2.9, 7.5, 0.13, BP=1, type=2, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen <- drug_pred(0.8, 7.35, 0.78, BP=1, type=2, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig <- drug_pred(1.48, 0, 0.87, BP=1, type=1, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho <- drug_pred(0.77, 0, 0.5, BP=1, type=1, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo <- drug_pred(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Zwitterion")
# Enoxacin
eno <- drug_pred(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_uni, dat_RR, "R&R") %>%
  mutate(Drug_class = "Zwitterion")


# Bind all the outputs together
drug_all_RR <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_RR <- round(cor(drug_all_RR$uni_pred,drug_all_RR$pred,method='pearson'), digits=2)
corr_coeff_RR <- paste("italic(r) %~~%", corr_coeff_RR)

fig4c <- ggplot() +
  geom_point(data=drug_all_RR, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_RR[1:2])) +
  ylim(0,max(drug_all_RR[1:2])) +
  xlab("Predicted Kp (standardized tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  ggtitle("c   RR") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_RR[1:2])*pcc_x, y=max(drug_all_RR[1:2])*pcc_y, label=corr_coeff_RR), parse=TRUE)


# Schmitt method
# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: adipose, bone, brain, gut, heart, kidney, liver, lung, muscle, skin, spleen, testes, RBCs

### Strong Base ###
# Metoprolol
met <- drug_pred(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace <- drug_pred(1.4, 9.7, 0.79, BP=1, type=3, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori <- drug_pred(2.56, 1.76, 0.42, BP=1, type=3, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp <- drug_pred(2.5, 2.4, 0.35, BP=1, type=3, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio <- drug_pred(2.9, 7.5, 0.13, BP=1, type=2, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen <- drug_pred(0.8, 7.35, 0.78, BP=1, type=2, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig <- drug_pred(1.48, 0, 0.87, BP=1, type=1, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho <- drug_pred(0.77, 0, 0.5, BP=1, type=1, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo <- drug_pred(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Zwitterion")
#Enoxacin
eno <- drug_pred(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_uni, dat_Schmitt, "Schmitt") %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_schmitt <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_Schmitt <- round(cor(drug_all_schmitt$uni_pred,drug_all_schmitt$pred,method='pearson'), digits=2)
corr_coeff_Schmitt <- paste("italic(r) %~~%", corr_coeff_Schmitt)

fig4d <- ggplot() +
  geom_point(data=drug_all_schmitt, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_schmitt[1:2])) +
  ylim(0,max(drug_all_schmitt[1:2])) +
  xlab("Predicted Kp (standardized tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  ggtitle("d   Schmitt") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_schmitt[1:2])*pcc_x, y=max(drug_all_schmitt[1:2])*pcc_y, label=corr_coeff_Schmitt), parse=TRUE)


#PK-Sim standard method
# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: bone, brain, adipose, heart, kidney, stomach, small intestine, large intestine, liver, lung, muscle, pancreas, skin, spleen

### Strong Base ###
# Metoprolol
met <- drug_pred(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Strong base")
# Acebutolol-R
ace <- drug_pred(1.4, 9.7, 0.79, BP=1, type=3, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Strong base")

### Weak Base ###
# Voriconazole
vori <- drug_pred(2.56, 1.76, 0.42, BP=1, type=3, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Weak base")
#Alprazolam
alp <- drug_pred(2.5, 2.4, 0.35, BP=1, type=3, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Weak base")

### Acid ###
# Thiopental
thio <- drug_pred(2.9, 7.5, 0.13, BP=1, type=2, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Acid")
# Phenobarbital
phen <- drug_pred(0.8, 7.35, 0.78, BP=1, type=2, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Acid")

### Neutral ###
# Digoxin
dig <- drug_pred(1.48, 0, 0.87, BP=1, type=1, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Neutral")
# Ethoxybenzamide
etho <- drug_pred(0.77, 0, 0.5, BP=1, type=1, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Neutral")

### Zwitterion ###
# Ofloxacin
oflo <- drug_pred(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Zwitterion")
# Enoxacin
eno <- drug_pred(0.1, c(6.1,8.7), 0.66, BP=0.94, type=6, dat_uni, dat_pksim, "PK-Sim") %>%
  mutate(Drug_class = "Zwitterion")

# Bind all the outputs together
drug_all_pksim <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff_pksim <- round(cor(drug_all_pksim$uni_pred,drug_all_pksim$pred,method='pearson'), digits=2)
corr_coeff_pksim <- paste("italic(r) %~~%", corr_coeff_pksim)

fig4e <- ggplot() +
  geom_point(data=drug_all_pksim, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all_pksim[1:2])) +
  ylim(0,max(drug_all_pksim[1:2])) +
  xlab("Predicted Kp (standardized tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  ggtitle("e   PK-Sim") +
  theme(legend.title =  element_blank()) +
  th1 +
  geom_text(aes(x=max(drug_all_pksim[1:2])*pcc_x, y=max(drug_all_pksim[1:2])*pcc_y, label=corr_coeff_pksim), parse=TRUE)


fig4 <- grid.arrange(fig4a, fig4b, fig4c, fig4d, fig4e, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig4.png", fig4, width=8, height=12)

#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 4: Figure 5 ##############################################
#########################################################################################################
## This figure shows the range of Kp predictions for each drug class for each tissue using different 
## calculation methods

# Calculate Kp values using each method for any drug
kp_pred_comp <- function(logP, pKa, fup, BP, type, dat_uni){
  drug_PT <- as.numeric(calcKp_PT(logP, pKa, fup, BP, type, dat_uni))
  drug_Berez <- as.numeric(calcKp_Berez(logP, pKa, fup, BP, type, dat_uni)) 
  drug_RR <- as.numeric(calcKp_RR(logP, pKa, fup, BP, type, dat_uni)) 
  drug_Schmitt <- as.numeric(calcKp_Schmitt(logP, pKa, fup, type, dat_uni)) 
  drug_pksim <- as.numeric(calcKp_pksim(logP, fup, dat_uni)) 
  
  drug_kp_all <- cbind(drug_PT, drug_Berez, drug_RR, drug_Schmitt, drug_pksim)
  drug_var <- apply(drug_kp_all, 1, var)
  
  drug_kp_all <- round(cbind(drug_kp_all, drug_var), digits=3)
  rownames(drug_kp_all) <- c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")
  
  drug_PT_new <- drug_PT %>%
    as.data.frame() %>%
    mutate(Method="PT") %>%
    mutate(Tissue=c(1:11))
  
  drug_Berez_new <- drug_Berez %>%
    as.data.frame() %>%
    mutate(Method="Berez") %>%
    mutate(Tissue=c(1:11))
  
  drug_RR_new <- drug_RR %>%
    as.data.frame() %>%
    mutate(Method="RR") %>%
    mutate(Tissue=c(1:11))
  
  drug_Schmitt_new <- drug_Schmitt %>%
    as.data.frame() %>%
    mutate(Method="Schmitt") %>%
    mutate(Tissue=c(1:11))
  
  drug_pksim_new <- drug_pksim %>%
    as.data.frame() %>%
    mutate(Method="PK-Sim") %>%
    mutate(Tissue=c(1:11))
  
  drug_kp <- rbind(drug_PT_new,drug_Berez_new,drug_RR_new,drug_Schmitt_new,drug_pksim_new)
  names(drug_kp)[1] <- "Kp"
  return(drug_kp)
}

### Strong Bases ###
# Metoprolol
met_kp <- kp_pred_comp(2.15, 9.7, 0.8, BP=1.52, type=3, dat_uni)
fig5_met <- ggplot(data=met_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=met_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(met_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("a    Metoprolol") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5

# Caffeine
caf_kp <- kp_pred_comp(-0.07, 10.4, 0.681, BP=0.98, type=3, dat_uni)
fig5_caf <- ggplot(data=caf_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=caf_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(caf_kp$Kp)) +
  xlab("") +
  ylab("Predicted Kp") +
  ggtitle("a    Caffeine") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5


### Weak Bases ###
# Voriconazole
vori_kp <- kp_pred_comp(2.56, 1.76, 0.42, BP=1, type=3, dat_uni)
fig5_vori <- ggplot(data=vori_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=vori_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(vori_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("b    Voriconazole") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5

# Alfentanil
alf_kp <- kp_pred_comp(2.2, 6.5, 0.11, BP=0.63, type=3, dat_uni)
fig5_alf <- ggplot(data=alf_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=alf_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(alf_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("b    Alfentanil") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5
fig5_alf

# Midazolam
mid_kp <- kp_pred_comp(3.1, 6, 0.059, BP=1, type=3, dat_uni)
fig5_mid <- ggplot(data=mid_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=mid_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(mid_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("c    Midazolam") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5

# Nevirapine
nev_kp <- kp_pred_comp(1.93, 2.8, 0.4, BP=1.04, type=3, dat_uni)
fig5_nev <- ggplot(data=nev_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=nev_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(nev_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("d    Nevirapine") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5


### Acid ###
# Thiopental
thio_kp <- kp_pred_comp(2.9, 7.5, 0.13, BP=1, type=2, dat_uni)
fig5_thio <- ggplot(data=thio_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=thio_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(thio_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("e    S-Thiopental") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5

# Nifedipine
nif_kp <- kp_pred_comp(2.2, 3.93, 0.04, BP=0.73, type=2, dat_uni)
fig5_nif <- ggplot(data=nif_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=nif_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(nif_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("c    Nifedipine") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5


### Neutral ###
# Digoxin
dig_kp <- kp_pred_comp(1.48, 0, 0.87, BP=1, type=1, dat_uni)
fig5_dig <- ggplot(data=dig_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=dig_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(dig_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("d    Digoxin") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5

# Artemether
art_kp <- kp_pred_comp(3.28, 0, 0.046, BP=0.8, type=1, dat_uni)
fig5_art <- ggplot(data=art_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=art_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(art_kp$Kp)) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("f    Artemether") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5


### Zwitterion ###
# Ofloxacin
oflo_kp <- kp_pred_comp(-0.4, c(5.97,9.28), 0.77, BP=0.92, type=6, dat_uni) 
fig5_oflo <- ggplot(data=oflo_kp, aes(x=Tissue, y=Kp)) +
  geom_point(data=oflo_kp, aes(x=Tissue, y=Kp, shape=Method),size=2.5, stroke=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels = c("Adipose","Bone","Brain","Heart","Kidney","Gut","Liver","Lung","Muscle","Skin","Spleen")) +
  ylim(0,max(oflo_kp$Kp)+1.4) +
  xlab("Tissue") +
  ylab("Predicted Kp") +
  ggtitle("e    Ofloxacin") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, fatten=0.1) +
  th5 

fig5 <- grid.arrange(fig5_met, fig5_vori, fig5_nif, fig5_dig, fig5_oflo, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig5.png", fig5, width=8, height=12)

figS1 <- grid.arrange(fig5_caf, fig5_alf, fig5_mid, fig5_nev, fig5_thio, fig5_art, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/figS1.png", figS1, width=8, height=12)

#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 5: Figure 6 ##############################################
#########################################################################################################
## This figure shows the impact of calculation methods on the model predictions from each drug class 
## representative

### Strong bases ###
# Metoprolol
source("PBPK_sim_metoprolol.R")
plot_met <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=mean), size=2) +
  geom_errorbar(data=df, aes(x=time, ymin=mean-sd, ymax=mean+sd), width=.1) +
  xlim(0,13) +
  scale_y_log10() +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("a  Metoprolol") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6 


# Caffeine
source("PBPK_sim_caffeine.R")
plot_caf <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=conc), size=2) +
  geom_errorbar(data=df, aes(x=time, ymin=pmax(conc-(sd_max-conc),0), ymax=sd_max), width=.1) +
  xlim(-0.1,24.1) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("a  Caffeine") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


### Weak bases ###
# Voriconazole 
source("PBPK_sim_voriconazole.R")
plot_vori <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("b  Voriconazole") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6

# Alfentanil
source("PBPK_sim_alfentanil.R")
plot_alf <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method), size=2, stroke=0.6) +
  geom_point(data=df_mean, aes(x=time, y=mean), size=2) +
  geom_errorbar(data=df_mean, aes(x=time, ymin=pmax(mean-SD,0), ymax=mean+SD), width=.1) +
  xlim(-0.1,12.1) +
  scale_y_log10() +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("b  Alfentanil") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6

# Midazolam 
source("PBPK_sim_midazolam.R")
plot_mid <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.5) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin = conc-(max_sd-conc), ymax = max_sd), width=.1) +
  xlim(0,6.1) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("c  Midazolam") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6

# Nevirapine
source( "PBPK_sim_nevirapine.R")
plot_nev <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data = df, aes(time, conc), size=2.5) +
  xlim(0,100) +
  #ylim(0,350) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("d  Nevirapine") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


### Acids ###
# Thiopental
source("PBPK_sim_thiopental.R")
plotS <-  ggplot() +
  geom_line(data=out_all_S,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced_S,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  xlim(0,15) +
  scale_y_log10() +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("e  S-Thiopental") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6
plotR <-  ggplot() +
  geom_line(data=out_all_R,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced_R,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df_R, aes(x=time, y=conc), size=2.5) +
  xlim(0,15) +
  scale_y_log10() +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("f  R-Thiopental") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6
# Plot both R-thiopental and S-thiopental together 
plot_thio <- grid.arrange(plotR, plotS, ncol=2)

# Nifedipine 
source("PBPK_sim_nifedipine.R")
plot_nif <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data = df, aes(time, conc), size=2.5) +
  geom_errorbar(data = df, aes(time, ymin = conc-(max_sd-conc), ymax = max_sd), width=.1) +
  xlim(0,5.7) +
  ylim(0,1.65) +
  xlab("Time (h)") +
  ylab("Plasma concentraction (mg/L)") +
  ggtitle("c  Nifedipine") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


### Neutrals ###
# Digoxin
source("PBPK_sim_digoxin.R")
plot_dig <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  xlim(0,6.12) +
  scale_y_log10() +
  xlab("Time (h)") +
  ylab("Percentage dose/L") +
  ggtitle("d  Digoxin") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


# Artemether
source("PBPK_sim_artemether.R")
plot_art <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data = df, aes(time, conc), size=2.5) +
  geom_errorbar(data = df, aes(time, ymin = pmax(conc-(max_sd-conc),0), ymax = max_sd), width=.1) +
  xlim(0,16.2) +
  ylim(0,0.35) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("g  Artemether") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


### Zwitterion ###
# Ofloxacin
source("PBPK_sim_ofloxacin.R")
plot_oflo <- ggplot() +
  geom_line(data=out_all,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
  geom_point(data=out_reduced,aes(x=time, y=Cplasma, shape=Method),size=2, stroke=0.6) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data=df, aes(x=time, ymin=conc-(sd_max-conc), ymax=sd_max), width=.1) +
  xlim(0,25) + 
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("e  Ofloxacin") +
  scale_linetype_manual("", values=c(1,1,1,1,1,1,1)) +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  scale_color_manual("", values=c('grey60','grey60','grey60','grey60','grey60','grey60')) +
  th6


fig6 <- grid.arrange(plot_met, plot_vori, plot_nif, plot_dig, plot_oflo, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig6.png", fig6, width=8, height=12)

figS2 <- grid.arrange(plot_caf, plot_alf, plot_mid, plot_nev, plotS, plotR, plot_art, ncol=2, nrow=4) 
#ggsave(file="../deliv/figure/figS2.png", figS2, width=8, height=16)

#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 6: Figure 8 ##############################################
#########################################################################################################
## This figure shows prediction ranges for representative drugs after incorporating inter-individual 
## variability and uncertainty in physiological and drug-related parameters, respectively

# Voriconazole
source("Variability_sim_voriconazole.R")
fig8 <- grid.arrange(pred_PT, pred_Berez, pred_RR, pred_Schmitt, pred_pksim, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig8.png", fig8, width=8, height=12)


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 7: Table 1 and Figure 7 ##################################
#########################################################################################################
## This table shows the error estimates for PK parameters (AUC and half-life) and RMSE from the model 
## predictions from each drug class representative

source("PBPK_sim_voriconazole.R") # pk_vori
source("PBPK_sim_alfentanil.R")   # pk_alf
source("PBPK_sim_nevirapine.R")   # pk_nev
source("PBPK_sim_midazolam.R")    # pk_mid
source("PBPK_sim_metoprolol.R")   # pk_met
source("PBPK_sim_caffeine.R")     # pk_caf
source("PBPK_sim_thiopental.R")   # pk_thio_S and pk_thio_R
source("PBPK_sim_nifedipine.R")   # pk_nif
source("PBPK_sim_digoxin.R")      # pk_dig
source("PBPK_sim_artemether.R")   # pk_art
source("PBPK_sim_ofloxacin.R")    # pk_oflo


# Make a figure, but could also convert to table
size <- 2.5
stroke <- 0.5


# RelRMSE = scale the root mean square error by the range of the observations
fig7a <- ggplot() +
  geom_point(data=pk_met, aes(x=1, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_caf, aes(x=2, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_vori, aes(x=4, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_alf, aes(x=5, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nev, aes(x=6, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_mid, aes(x=7, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_S, aes(x=9, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_R, aes(x=10, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nif, aes(x=11, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_dig, aes(x=13, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_art, aes(x=14, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_oflo, aes(x=16, y=RelRMSE, shape=Method),size=size, stroke=stroke) +
  scale_x_continuous(breaks=c(1,2,4,5,6,7,9,10,11,13,14,16),
                     labels = c("Metoprolol","Caffeine","Voriconazole","Alfentanil",
                                "Nevirapine","Midazolam","S-Thiopental","R-Thiopental",
                                "Nifedipine","Digoxin","Artemether","Ofloxacin")) +
  #ylim(0,18) + #Using error
  #ylim(0, 1500) +
  scale_y_log10(limits = c(1e-1,1e6)) +
  xlab("") +
  ylab("Percent error") +
  ggtitle("a  RMSE") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  th7
fig7a

### AUC error
fig7b <- ggplot() +
  geom_point(data=pk_met, aes(x=1, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_caf, aes(x=2, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_vori, aes(x=4, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_alf, aes(x=5, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nev, aes(x=6, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_mid, aes(x=7, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_S, aes(x=9, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_R, aes(x=10, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nif, aes(x=11, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_dig, aes(x=13, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_art, aes(x=14, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_oflo, aes(x=16, y=AUCerror, shape=Method),size=size, stroke=stroke) +
  #geom_abline(intercept = 0, slope = 0) + # Using error
  scale_x_continuous(breaks=c(1,2,4,5,6,7,9,10,11,13,14,16),
                     labels = c("Metoprolol","Caffeine","Voriconazole","Alfentanil",
                                "Nevirapine","Midazolam","S-Thiopental","R-Thiopental",
                                "Nifedipine","Digoxin","Artemether","Ofloxacin")) +
  #ylim(-27,26) + #Using error
  #ylim(0, 4200) +
  scale_y_log10(limits = c(1e-1,1e6)) +
  xlab("") +
  ylab("Percent error") +
  ggtitle("b  AUC error") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  th7
fig7b

### Half-life (semi-log plot)
fig7c <- ggplot() +
  geom_point(data=pk_met, aes(x=1, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_caf, aes(x=2, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_vori, aes(x=4, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_alf, aes(x=5, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nev, aes(x=6, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_mid, aes(x=7, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_S, aes(x=9, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_thio_R, aes(x=10, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_nif, aes(x=11, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_dig, aes(x=13, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_art, aes(x=14, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  geom_point(data=pk_oflo, aes(x=16, y=abs(hlerror), shape=Method),size=size, stroke=stroke) +
  scale_x_continuous(breaks=c(1,2,4,5,6,7,9,10,11,13,14,16),
                     labels = c("Metoprolol","Caffeine","Voriconazole","Alfentanil",
                                "Nevirapine","Midazolam","S-Thiopental","R-Thiopental",
                                "Nifedipine","Digoxin","Artemether","Ofloxacin")) +
  scale_y_log10(limits = c(1e-1,1e6)) +
  xlab("") +
  ylab("Percent error") +
  ggtitle("c  Half-life error") +
  scale_shape_manual("", values=c(0,2,3,4,5,8)) +
  th7
fig7c

fig7 <- grid.arrange(fig7a, fig7b, fig7c, ncol=3, nrow=1)
#ggsave(file="../deliv/figure/fig7.png", fig7, width=8, height=6)

### Half-life
# fig7c_no_log <- ggplot() +
#   geom_point(data=pk_met, aes(x=1, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_caf, aes(x=2, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_vori, aes(x=4, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_alf, aes(x=5, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_nev, aes(x=6, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_mid, aes(x=7, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_thio_S, aes(x=9, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_thio_R, aes(x=10, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_nif, aes(x=11, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_dig, aes(x=13, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_art, aes(x=14, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_point(data=pk_oflo, aes(x=16, y=hlerror, shape=Method),size=size, stroke=stroke) +
#   geom_abline(intercept = 0, slope = 0) +
#   scale_x_continuous(breaks=c(1,2,4,5,6,7,9,10,11,13,14,16),
#                      labels = c("Metoprolol","Caffeine","Voriconazole","Alfentanil",
#                                 "Nevirapine","Midazolam","S-Thiopental","R-Thiopental",
#                                 "Nifedipine","Digoxin","Artemether","Ofloxacin")) +
#   ylim(-80,80) +
#   xlab("") +
#   ylab("Residual error") +
#   ggtitle("c  Half-life error") +
#   scale_shape_manual("", values=c(0,2,3,4,5,8)) +
#   th7
# 
# fig7_no_log <- grid.arrange(fig7a, fig7b, fig7c_no_log, ncol=3, nrow=1)
# ggsave(file="../deliv/figure/fig7_no_log.png", fig7_no_log, width=8, height=6)




### Generate table (drugs in the rows, columns of RMSE, AUC, and half-life error)
pk_met_mod <- c(pk_met$RelRMSE,pk_met$AUCerror,pk_met$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Metoprolol")
pk_caf_mod <- c(pk_caf$RelRMSE,pk_caf$AUCerror,pk_caf$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Caffeine")

pk_vori_mod <- c(pk_vori$RelRMSE,pk_vori$AUCerror,pk_vori$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Voriconazole")
pk_alf_mod <- c(pk_alf$RelRMSE,pk_alf$AUCerror,pk_alf$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Alfentanil")
pk_nev_mod <- c(pk_nev$RelRMSE,pk_nev$AUCerror,pk_nev$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Nevirapine")
pk_mid_mod <- c(pk_mid$RelRMSE,pk_mid$AUCerror,pk_mid$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Midazolam")

pk_thio_S_mod <- c(pk_thio_S$RelRMSE,pk_thio_S$AUCerror,pk_thio_S$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("S-Thiopental")
pk_thio_R_mod <- c(pk_thio_R$RelRMSE,pk_thio_R$AUCerror,pk_thio_R$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("R-Thiopental")
pk_nif_mod <- c(pk_nif$RelRMSE,pk_nif$AUCerror,pk_nif$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Nifedipine")

pk_dig_mod <- c(pk_dig$RelRMSE,pk_dig$AUCerror,pk_dig$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Digoxin")
pk_art_mod <- c(pk_art$RelRMSE,pk_art$AUCerror,pk_art$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Artemether")

pk_oflo_mod <- c(pk_oflo$RelRMSE,pk_oflo$AUCerror,pk_oflo$hlerror) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Ofloxacin")

pk_mod <- rbind(pk_met_mod, pk_caf_mod, pk_vori_mod, pk_alf_mod, pk_nev_mod, pk_mid_mod, 
                pk_thio_S_mod, pk_thio_R_mod, pk_nif_mod, pk_dig_mod, pk_art_mod, pk_oflo_mod)

# pk_mod <- rbind(pk_met_mod, pk_caf_mod, pk_vori_mod, pk_alf_mod, pk_nev_mod, pk_mid_mod, 
#                 pk_thio_S_mod, pk_thio_R_mod, pk_nif_mod, pk_dig_mod, pk_art_mod, pk_oflo_mod) %>%
#   'colnames<-'(c("PT","Berez","RR","Schmitt","PK-Sim",
#                  "PT","Berez","RR","Schmitt","PK-Sim",
#                  "PT","Berez","RR","Schmitt","PK-Sim"))

# Include row of average values for each method (over all drugs)
pk_mod_avg <- rbind(c(pk_met$RelRMSE,pk_met$AUCerror,pk_met$hlerror), c(pk_caf$RelRMSE,pk_caf$AUCerror,pk_caf$hlerror),
                    c(pk_vori$RelRMSE,pk_vori$AUCerror,pk_vori$hlerror), c(pk_alf$RelRMSE,pk_alf$AUCerror,pk_alf$hlerror),
                    c(pk_nev$RelRMSE,pk_nev$AUCerror,pk_nev$hlerror), c(pk_mid$RelRMSE,pk_mid$AUCerror,pk_mid$hlerror),
                    c(pk_thio_S$RelRMSE,pk_thio_S$AUCerror,pk_thio_S$hlerror), c(pk_thio_R$RelRMSE,pk_thio_R$AUCerror,pk_thio_R$hlerror),
                    c(pk_nif$RelRMSE,pk_nif$AUCerror,pk_nif$hlerror), c(pk_dig$RelRMSE,pk_dig$AUCerror,pk_dig$hlerror),
                    c(pk_art$RelRMSE,pk_art$AUCerror,pk_art$hlerror), c(pk_oflo$RelRMSE,pk_oflo$AUCerror,pk_oflo$hlerror))

pk_avg <-  colMeans(pk_mod_avg) %>%
  sig() %>%
  t() %>%
  'rownames<-'("Mean")

pk_mod_all <- rbind(pk_mod, pk_avg) %>%
   'colnames<-'(c("PT","Berez","RR","Schmitt","PK-Sim",
                  "PT","Berez","RR","Schmitt","PK-Sim",
                  "PT","Berez","RR","Schmitt","PK-Sim"))

table1 <- kable(pk_mod_all,"html") %>%
  kable_styling(full_width=F) %>%
  add_header_above(c("","Relative percent RMSE"=5,"AUC percent error"=5, "Half-life percent error"=5)) %>%
  group_rows("Strong bases", 1, 2) %>%
  group_rows("Weak bases", 3, 6) %>%
  group_rows("Acids", 7, 9) %>%
  group_rows("Neutrals", 10, 11) %>%
  group_rows("Zwitterion", 12, 12) 


#########################################################################################################
#########################################################################################################


