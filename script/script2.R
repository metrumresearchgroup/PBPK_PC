# This script reproduces manuscript figures (Ahmed)
#clear environment
rm(list=ls())
# load libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(mrgsolve)
library(gridExtra)
source("CalcKp_P&T.R")
source("CalcKp_R&R.R")
source("CalcKpu_R&R.R")
source("CalcKp_Berez.R")
source("CalcKp_Schmitt.R")
source("CalcKp_pksim.R")
source("CalcVss_P&T.R")

## set functions ##
filter <- dplyr::filter
mutate <- dplyr::mutate
select <- dplyr::select

#########################################################################################################
######################################## CHUNK 1: Figure 1 ##############################################
#########################################################################################################
## This figure validates our code by comparing our Vss/Kp predictions (using literature physiology) 
## to the reported values in literature
## note: we will stick to the original reported values; eg: Vss for PT and Kpu for RR
## note: we will use as many of the reported drugs as possible for validation trying to include the drugs
## we investigated in simulation if their Kp/Vss are reported

# set figure default theme
th1 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 10),
             legend.justification=c(0.025,0.95), legend.position=c(0.025,0.95), legend.key=element_blank(), 
             legend.title = element_blank(),
             plot.title = element_text(face="bold", size=15))


## Poulin and Theil figure ##
# load data from paper
dat <- read.csv("../data/tissue_comp_P&T.csv")

# Acids:
# Acetimenophen, Ascorbic acid, Isoniazid, Thiopental, Acetozolamid, Tetrahydrocam
pred_PT_a <- c(calcVss_PT(0.46, 9.38, 1, BP=1, type=2), 
               calcVss_PT(-1.85, 4.7, 0.76, BP=0.82, type=2), 
               calcVss_PT(-0.7, 1.82, 1, BP=1, type=2), 
               calcVss_PT(2.85, 7.45, 0.112, BP=1, type=2), 
               calcVss_PT(-0.26, 7.4, 0.044, BP=1, type=2), 
               calcVss_PT(6.97, 10.6, 0.115, BP=0.55, type=2))

exp_PT_a  <- c(0.85, 0.72, 0.56, 2, 0.38,8.95)  # experiemntal
rep_PT_a <- c(0.63,0.5,0.56,2.87, 0.32,5.53)  # reported predicted

# Strong Bases: 
# Caffeine, Acebutolol, Chlorpheniramine, Carvedilol, Mepivacaine, Trimethoprim
pred_PT_b <- c(calcVss_PT(-0.07, 10.4, 0.681, BP=0.98, type=3), calcVss_PT(1.71, 9.67, 0.743, BP=1, type=3), calcVss_PT(3.38, 9.13, 0.27, BP=1.34, type=3), calcVss_PT(4.19, 8.4, 0.0054, BP=0.71, type=3), calcVss_PT(1.95, 7.76, 0.225, BP=0.9, type=3), calcVss_PT(0.91, 7.12, 0.304, BP=1, type=3))

#exp_PT_b <- c(0.73,1.17, 3.17, 1.54, 1.1, 1.96)    
rep_PT_b <- c(0.53,1.25,3.26,2.68,1.22,0.52)


#Weak Bases
# Alfentanil, Cimetidine, Etomidate
pred_PT_wb <- c(calcVss_PT(2.16, 6.5, 0.074, BP=1, type=3),calcVss_PT(0.4, 6.8, 0.8, BP=1, type=3),calcVss_PT(3.05, 4.5, 0.235, BP=1, type=3))
exp_PT_wb <- c(0.226,1,2.5)
rep_PT_wb <- c(1.38,0.56,4.55)

#Neutrals:
# Coumarin, Griseofulvin, Cyclosporin
pred_PT_n <- c(calcVss_PT(1.39, 0, 0.853, BP=1, type=1),calcVss_PT(2.18, 0, 0.2, BP=1, type=1),calcVss_PT(2.92, 0, 0.0602, BP=1.36, type=1))

exp_PT_n <- c(1.29, 1.4, 1.89)  
rep_PT_n <- c(1.08,1.72,2.77)


#Zwitterions:
#Ciprofloxacin, Enoxacin 
pred_PT_z <-c(calcVss_PT(-0.2,c(6.3,8.6),0.485,1,type=6),calcVss_PT(0.28,c(6.3,8.3),0.52,1,type=6))
rep_PT_z <- c(0.47,0.44)


# Data processing
sbase <- data.frame(pred = unlist(pred_PT_b),rep_pred = rep_PT_b)
sbase <- mutate(sbase, Drug_class = "Strong base")

wbase <- data.frame(pred = unlist(pred_PT_wb),rep_pred = rep_PT_wb)
wbase <- mutate(wbase, Drug_class = "Weak base")

acid <- data.frame(pred = unlist(pred_PT_a),rep_pred = rep_PT_a)
acid <- mutate(acid, Drug_class = "Acid")

neut <- data.frame(pred = unlist(pred_PT_n),rep_pred = rep_PT_n)
neut <- mutate(neut, Drug_class = "Neutral")

zwit <- data.frame(pred = unlist(pred_PT_z),rep_pred = rep_PT_z)
zwit <- mutate(zwit, Drug_class = "Zwitterion")

drug_all <- rbind(sbase,wbase,acid,neut,zwit)


# Calculate the Pearson correlation coefficient 
corr_coeff <- round(cor(drug_all$pred,drug_all$rep_pred,method='pearson'), digits=2)
corr_coeff <- paste("italic(r) %~~%", corr_coeff)


gp1 <- ggplot() +
  geom_point(data=drug_all, aes(x=pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$rep_pred)) +
  xlab("Predicted Vss") +
  ylab("Reported Vss") +
  ggtitle("a    Poulin and Theil") + 
  th1 +
  geom_text(aes(x=3, y=2, label=corr_coeff), parse=TRUE)
gp1



## Berez figure ##
# load data from paper
dat <- read.csv("../data/PKSim_tissue_comp_PT_Berez.csv")  # the tissue composition as reported by PK-Sim

### Strong base
#Quinidine (logMA=2.3)
#pred_quin <- calcKp_Berez(2.3, 8.6, 0.33,type=3)
#pksim_pred_quin <- c()

#Acebutolol-R
pred_ace <- calcKp_Berez(1.4,9.7,0.79,type=3)
pksim_pred_ace <- c(0.16,1.82,2.01,0.99,1.18,1.67,1.67,1.67,1.48,0.82,1.17,1.46,1.24,1.2)

#pred_ace_p <- calcKp_Berez(1.87,9.7,0.79,type=3)
#pksim_pred_ace_p <- c(5.06,1.57,20.55,1.06,15.17,12.13,5.99,5.99,5.99,13.88,6.78,2.17,3.23,4.91) #using logP=1.87


### Weak base
#Alprazolam
pred_alp <- calcKp_Berez(2.5,2.4,0.35,type=3)
pksim_pred_alp <- c(7.75,6.06,5.72,1.66,2.39,4.56,4.56,4.56,3.7,0.84,2.41,3.76,2.84,2.43)


### Acid
#Phenobarbital (logMA)
pred_phen <- calcKp_Berez(0.8,7.35,0.78, BP=0.861, type=2)
pksim_pred_phen <- c(0.26,0.78,1.09,0.79,0.86,0.94,0.94,0.94,0.92,0.78,0.84,0.88,0.83,0.87)


### Neutral
pred_dig <- calcKp_Berez(1.48,0,0.87,type=1)
pksim_pred_dig <- c(1.61,2.24,2.4,1.09,1.34,1.98,1.98,1.98,1.73,0.87,1.32,1.72,1.43,1.36)


### Zwitterion
#Ofloxacin
pred_oflo <- calcKp_Berez(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6)
pksim_pred_oflo <- c(0.15,0.43,0.77,0.72,0.75,0.7,0.7,0.7,0.73,0.76,0.72,0.68,0.69,0.76)



# Data processing
ace <- data.frame(pred = unlist(pred_ace),rep_pred = pksim_pred_ace)
ace <- mutate(ace, Drug_class = "Strong base")

alp <- data.frame(pred = unlist(pred_alp),rep_pred = pksim_pred_alp)
alp <- mutate(alp, Drug_class = "Weak base")

phen <- data.frame(pred = unlist(pred_phen),rep_pred = pksim_pred_phen)
phen <- mutate(phen, Drug_class = "Acid")

dig <- data.frame(pred = unlist(pred_dig),rep_pred = pksim_pred_dig)
dig <- mutate(dig, Drug_class = "Neutral")

oflo <- data.frame(pred = unlist(pred_oflo),rep_pred = pksim_pred_oflo)
oflo <- mutate(oflo, Drug_class = "Zwitterion")

drug_all <- rbind(ace,alp,phen,dig,oflo)

# Calculate the Pearson correlation coefficient 
corr_coeff <- round(cor(drug_all$pred,drug_all$rep_pred,method='pearson'), digits=2)
corr_coeff <- paste("italic(r) %~~%", corr_coeff)


gp2 <- ggplot() +
  geom_point(data=drug_all, aes(x=pred, y=rep_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$rep_pred)) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("b    Berezhkovskiy") + 
  th1 +
  geom_text(aes(x=4, y=2, label=corr_coeff), parse=TRUE)
gp2



## Rodgers and Rowland figure ##
# load data
dat <- read.csv("../data/tissue_comp_R&R.csv")

# Acids:
# Cefazolin
pred_RR_a_all <- calcKpu_RR(0.3, 2.3, 0.16, BP=0.55, type=2)
pred_RR_a <- c(pred_RR_a_all[1],pred_RR_a_all[3:8],pred_RR_a_all[10])

#exp_RR_a  <- c(0.57,1.11,0.43,17.1,4.84,1.03,0.4,1.89)
rep_RR_a <- c(0.74,1.28,1.23,1.1,0.79,1.52,0.69,1.92)

# Phenobarbital
pred_RR_a_all2 <- calcKpu_RR(1.6,7.4,0.64,BP=0.861,type=2)
pred_RR_a2 <- c(pred_RR_a_all2[3:8],pred_RR_a_all2[10],pred_RR_a_all2[13])

#exp_RR_a2  <- c(2.82,1.42,1.14,2.82,1.21,1.55,1.88,0.47)
rep_RR_a2 <- c(1.43,0.88,0.97,0.92,1.07,0.79,1.84,1.23)

#Phenytoin
#pred_RR_a_all3 <- calcKpu_RR(2.5,8.3,0.19,BP=0.19,type=2)
#pred_RR_a3 <- c(pred_RR_a_all3[2:8],pred_RR_a_all3[10],pred_RR_a_all3[13])

#exp_RR_a3  <- c(4.7,12.9,6.2,8.27,11.9,4.91,5.68,8.79,9.3)
#rep_RR_a3 <- c(12.0,12.9,5.81,6.55,6.74,8.39,4.33,19.0,21.2)



# Bases: 
#Alfentanil (weak)
pred_RR_b_all <- calcKpu_RR(2.2,6.5,0.11,BP=1, type=3)
pred_RR_b <- c(pred_RR_b_all[2:11],pred_RR_b_all[13])

exp_RR_b <- c(1.14,19.2,5.0,7.47,9.05,7.03,2.78,8.65,1.65,6.65,19.1)   
rep_RR_b <- c(6.3,7.33,3.97,4.24,4.05,5.48,2.88,6.71,10.7,2.7,9.2)

#Strong bases:
# Bisoprolol-R
pred_RR_sb_all <- calcKpu_RR(1.87,9.4,0.85,BP=1.36,type=3)
pred_RR_sb <- c(pred_RR_sb_all[1:8],pred_RR_sb_all[10],pred_RR_sb_all[12:13])

exp_RR_sb <- c(5.47,1.93,31.2,7.63,29.3,26.8,49.2,6.35,2.56,10.1,1.21)   
rep_RR_sb <- c(3.61,3.36,11,10.4,21.4,19.6,17.2,7.73,6.37,10.8,1.79)

#Acebutolol-R
pred_RR_sb_all2 <- calcKpu_RR(1.87,9.7,0.79,BP=1.09,type=3)
pred_RR_sb2 <- c(pred_RR_sb_all2[1:8],pred_RR_sb_all2[10],pred_RR_sb_all2[12:13])

exp_RR_sb2 <- c(0.07,0.61,112,7.23,29.7,39.7,13.1,6.28,3.83,6.85,1.38)   
rep_RR_sb2 <- c(2.34,2.57,6.45,6.1,11.8,11.0,9.67,4.84,3.85,6.45,0.82)


#Neutral:
# Ethoxybenzamide
pred_RR_n_all <- calcKpu_RR(0.77, 0, 0.5, BP=1, type=1)

pred_RR_n <- c(pred_RR_n_all[2:5],pred_RR_n_all[7:8],pred_RR_n_all[10:11],pred_RR_n_all[13])  
exp_RR_n <- c(1.58,0.97,1.74,2.19,1.57,1.32,1.73,1.47,1.26)  
rep_RR_n <- c(1.03,1.1,0.79,0.91,0.86,0.85,1.11,0.78,0.48)


#Zwitterions:
#Enoxacin
pred_eno <- calcKpu_RR(log10(1.26),c(6.1,8.7),0.66,BP=0.94,type=6)
pred_eno <- c(pred_eno[1],pred_eno[4:8],pred_eno[10:11])
rep_eno <- c(1.74,4.62,9.88,9.04,7.61,3.74,3.16,6.33)
exp_eno <- c(2.18,1.62,6.99,4.87,1.73,2.2,2.06,2.47)

# #Nalidixic acid
# pred_nal <- calcKpu_RR(log10(12.6),c(5.1,3.3),0.29,type=6)
# pred_nal <- c(pred_nal[1:8],pred_nal[10:11],pred_nal[13])
# rep_nal <- c(0.47,0.52,0.86,0.8,0.75,0.57,0.95,0.52,1.17,0.59,0.27)
# 
# #Ceftazidime
# pred_ceft <- calcKpu_RR(log10(0.316),c(2.5,1.9,3.8),0.1,type=9)
# pred_ceft <- c(pred_ceft[3:8],pred_ceft[10],pred_ceft[13])
# rep_ceft <- c(0.48,0.44,0.45,0.36,0.45,0.37,0.52,0.14)


# Data processing
alf <- data.frame(pred = unlist(pred_RR_b),rep = rep_RR_b)
alf <- mutate(alf, Drug_class = "Weak base")

ace <- data.frame(pred = unlist(pred_RR_sb2),rep = rep_RR_sb2)
ace <- mutate(ace, Drug_class = "Strong base")

bis <- data.frame(pred = unlist(pred_RR_sb),rep = rep_RR_sb)
bis <- mutate(bis, Drug_class = "Strong base")

phen <- data.frame(pred = unlist(pred_RR_a2),rep = rep_RR_a2)
phen <- mutate(phen, Drug_class = "Acid")

cef <- data.frame(pred = unlist(pred_RR_a),rep = rep_RR_a)
cef <- mutate(cef, Drug_class = "Acid")

etho <- data.frame(pred = unlist(pred_RR_n),rep = rep_RR_n)
etho <- mutate(etho, Drug_class = "Neutral")

eno <- data.frame(pred = unlist(pred_eno),rep = rep_eno)
eno <- mutate(eno, Drug_class = "Zwitterion")

drug_all <- rbind(alf,ace,bis,phen,cef,etho,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff <- round(cor(drug_all$pred,drug_all$rep,method='pearson'), digits=2)
corr_coeff <- paste("italic(r) %~~%", corr_coeff)


gp3 <- ggplot() +
  geom_point(data=drug_all, aes(x=pred, y=rep, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$rep)) +
  xlab("Predicted Kpu") +
  ylab("Reported Kpu") +
  ggtitle("c    Rodgers and Rowland") + 
  th1 +
  geom_text(aes(x=12, y=5, label=corr_coeff), parse=TRUE)
gp3


## Schmitt fiure ##
# load data
# dat <- read.csv("../data/tissue_comp_Schmitt.csv")  # this doesn't work but pksim physiology does; ask Kiersten
dat <- read.csv("../data/PKSim_tissue_comp_Schmitt.csv")

#Using logMA below, could also use logP
### Strong base
#Quinidine (logMA=2.3)
pred_quin <- calcKp_Schmitt(2.3, 8.6, 0.33, type=3)
pksim_pred_quin <- c(1.6,22.23,2.85,16.4,12.89,6.22,6.22,6.22,14.85,7.17,1.99,3.26,5.32,5.17)

#Acebutolol-R
pred_ace <- calcKp_Schmitt(1.4,9.7,0.79,type=3)
pksim_pred_ace <- c(0.69,7.38,0.38,5.52,4.51,2.45,2.45,2.45,5.08,2.69,1.15,1.48,1.94,2.13)

#pred_ace_p <- calcKp_Schmitt(1.87,9.7,0.79,type=3)
#pksim_pred_ace_p <- c(1.57,20.55,1.06,15.17,12.13,5.99,5.99,5.99,13.88,6.78,2.17,3.23,4.91,5.06) #using logP=1.87

#pred_ace_RR <- calcKp_RR(1.4,9.7,0.79,BP=1.19,type=3)
#pksim_pred_ace_RR <- c(1.68,1.89,0.72,4,8.27,4.62,4.62,4.62,7.67,6.39,4.88,3.55,2.77,5.46)


### Weak base
#Alprazolam
pred_alp <- calcKp_Schmitt(2.5,2.4,0.35,type=3)
pksim_pred_alp <- c(2.83,12.65,102.27,12.9,7.42,8.4,8.4,8.4,9.53,4.96,1.85,6.21,16.68,3.03)


### Acid
#Phenobarbital (logMA)
pred_phen <- calcKp_Schmitt(0.8,7.35,0.78, type=2)
pksim_pred_phen <- c(0.32,1.08,2.95,1.02,0.87,0.85,0.85,0.85,0.9,0.78,0.68,0.73,0.96,0.72)

### Neutral
pred_dig <- calcKp_Schmitt(1.48,0,0.87,type=1)
pksim_pred_dig <- c(0.9,3.64,24.31,3.64,2.36,2.63,2.63,2.63,2.82,1.77,1.06,2.05,4.38,1.34)


### Zwitterion
# Ofloxacin
# logP = -0.4, pKa = c(5.97,9.28) , fup = 0.77, BP = 0.92, type = 6 (monoprotic acid, monoprotic base)
#source("../CalcKp_Schmitt.R")
pred_oflo <- calcKp_Schmitt(-0.4,c(5.97,9.28),0.77,type=6)
pksim_pred_oflo <- c(0.23,0.66,0.29,0.6,0.61,0.64,0.64,0.64,0.58,0.6,0.61,0.58,0.46,0.62)


# Data processing
quin <- data.frame(pred = unlist(pred_quin),pksim_pred = pksim_pred_quin)
quin <- mutate(quin, Drug_class = "Strong base")

ace <- data.frame(pred = unlist(pred_ace),pksim_pred = pksim_pred_ace)
ace <- mutate(ace, Drug_class = "Strong base")

alp <- data.frame(pred = unlist(pred_alp),pksim_pred = pksim_pred_alp)
alp <- mutate(alp, Drug_class = "Weak base")

phen <- data.frame(pred = unlist(pred_phen),pksim_pred = pksim_pred_phen)
phen <- mutate(phen, Drug_class = "Acid")

dig <- data.frame(pred = unlist(pred_dig),pksim_pred = pksim_pred_dig)
dig <- mutate(dig, Drug_class = "Neutral")

oflo <- data.frame(pred = unlist(pred_oflo),pksim_pred = pksim_pred_oflo)
oflo <- mutate(oflo, Drug_class = "Zwitterion")

drug_all <- rbind(quin,ace,alp,phen,dig,oflo)

# Calculate the Pearson correlation coefficient 
corr_coeff <- round(cor(drug_all$pred,drug_all$pksim_pred,method='pearson'), digits=2)
corr_coeff <- paste("italic(r) %~~%", corr_coeff)


gp4 <- ggplot() +
  geom_point(data=drug_all, aes(x=pred, y=pksim_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$pksim_pred)) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("d    Schmitt") + 
  th1 +
  geom_text(aes(x=50, y=25, label=corr_coeff), parse=TRUE)
gp4



## PK-Sim standard figure ##
# load data
dat <- read.csv("../data/PKSim_tissue_comp_pksim.csv")

### Strong Base
# Acebutolol-R
pred_ace <- calcKp_pksim(1.87,0.79)
pksim_pred_ace <- c(16.44,7.2,47.04,6.67,3.9,4.45,4.45,4.45,4.89,1.48,1.65,5.43,6.75,1.83)

### Weak Base
# Triazolam
pred_tri <- calcKp_pksim(2.4,0.28)
pksim_pred_tri <- c(19.41,8.09,56.39,7.51,4.15,4.79,4.79,4.79,5.36,1.22,1.42,6.03,7.66,1.65)

### Acids
# Phenobarbital
pred_phen <- calcKp_pksim(1.6,0.64)
pksim_pred_phen <- c(7.3,3.37,20.51,3.13,1.93,2.18,2.18,2.18,2.36,0.89,0.97,2.58,3.13,1.03)

###Neutral
# Ethoxybenzamide
pred_etho <- calcKp_pksim(0.77,0.5)
pksim_pred_etho <- c(1.06,0.74,2.44,0.68,0.57,0.6,0.6,0.6,0.6,0.46,0.47,0.61,0.64,0.46)


###Zwitterion
#Enoxacin
pred_eno <- calcKp_pksim(0.100,0.66)
pksim_pred_eno <- c(0.56,0.63,0.77,0.59,0.58,0.59,0.59,0.59,0.57,0.56,0.57,0.55,0.52,0.55)

# m <- 21
# plot(pred_phen,pksim_pred_phen,xlab = "Predicted Kp", ylab = "PK-Sim Kp", xlim = c(0,m), ylim = c(0,m), pch=0)
# par(new=T)
# plot(pred_tri,pksim_pred_tri,xlab = "Predicted Kp", ylab = "PK-Sim Kp", xlim = c(0,m), ylim = c(0,m), pch=2)
# par(new=T)
# plot(pred_ace,pksim_pred_ace,xlab = "Predicted Kp", ylab = "PK-Sim Kp", xlim = c(0,m), ylim = c(0,m), pch=1)
# par(new=T)
# plot(pred_etho,pksim_pred_etho,xlab = "Predicted Kp", ylab = "PK-Sim Kp", xlim = c(0,m), ylim = c(0,m), pch=3)
# par(new=T)
# plot(pred_eno,pksim_pred_eno,xlab = "Predicted Kp", ylab = "PK-Sim Kp", xlim = c(0,m), ylim = c(0,m), pch=8)
# title(main="Comparison of Kp")
# legend("topleft", legend = c("Acids", "Weak Bases", "Strong Bases", "Neutrals", "Zwitterions"), pch = c(0, 2, 1, 3, 8))
# abline(0,1)


# Data processing
ace <- data.frame(pred = unlist(pred_ace),pksim_pred = pksim_pred_ace)
ace <- mutate(ace, Drug_class = "Strong base")

tri <- data.frame(pred = unlist(pred_tri),pksim_pred = pksim_pred_tri)
tri <- mutate(tri, Drug_class = "Weak base")

phen <- data.frame(pred = unlist(pred_phen),pksim_pred = pksim_pred_phen)
phen <- mutate(phen, Drug_class = "Acid")

etho <- data.frame(pred = unlist(pred_etho),pksim_pred = pksim_pred_etho)
etho <- mutate(etho, Drug_class = "Neutral")

eno <- data.frame(pred = unlist(pred_eno),pksim_pred = pksim_pred_eno)
eno <- mutate(eno, Drug_class = "Zwitterion")

drug_all <- rbind(ace,alp,phen,dig,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff <- round(cor(drug_all$pred,drug_all$pksim_pred,method='pearson'), digits=2)
corr_coeff <- paste("italic(r) %~~%", corr_coeff)


gp5 <- ggplot() +
  geom_point(data=drug_all, aes(x=pred, y=pksim_pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$pksim_pred)) +
  xlab("Predicted Kp") +
  ylab("Reported Kp (PK-Sim)") +
  ggtitle("e    PK-Sim Standard") + 
  th1 +
  geom_text(aes(x=50, y=25, label=corr_coeff), parse=TRUE)
gp5

gp <- grid.arrange(gp1, gp2, gp3, gp4, gp5, ncol=2, nrow=3)
ggsave(file="../deliv/figure/fig1.pdf", gp, width=8, height=12)
#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 2: Figure 2 ##############################################
#########################################################################################################
## This figure swaps physiologies between calculation methods to see impact on Kp/Vss predictions


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 3: Figure 3 ##############################################
#########################################################################################################
## This figure compares Kp/Vss predictions using the unified physiology vs those with reported physiologies

### Poulin and Theil ###
# ## Predictions with unified physiology ##
# ### Strong Base ###
# #Metoprolol
# dat <- read.csv("../data/unified_tissue_comp.csv")
# met_uni <- calcKp_PT(2.15,9.7,0.8,BP=1.52,type=3) #prediction using unified physiology
# 
# # Acebutolol-R
# ace_uni <- calcKp_PT(1.4,9.7,0.79,type=3)
# 
# 
# ### Weak Base ###
# # Voriconazole
# vori_uni <- calcKp_PT(2.56,1.76,0.42,BP=1,type=3) 
# 
# #Alprazolam
# alp_uni <- calcKp_PT(2.5,2.4,0.35,type=3) 
# 
# 
# ### Acid ###
# # Thiopental
# thio_uni <- calcKp_PT(2.9,7.5,0.13, BP=1,type=2)
# 
# # Phenobarbital
# phen_uni <- calcKp_PT(0.8,7.35,0.78, type=2)
# 
# 
# ### Neutral ###
# # Digoxin
# dig_uni <- calcKp_PT(1.48,0,0.87,type=1)
# 
# # Ethoxybenzamide
# etho_uni <- calcKp_PT(0.77,0,0.5,BP=1,type=1)
# 
# 
# ### Zwitterion ###
# # Ofloxacin
# oflo_uni <- calcKp_PT(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6)
# 
# #Enoxacin
# eno_uni <- calcKp_PT(0.1,c(6.1,8.7),0.66,BP=0.94,type=6)
# 
# 
# 
# ## Predictions with reported physiology ##
# ### Strong Base ###
# #Metoprolol
# dat <- read.csv("../data/tissue_comp_P&T.csv")
# met_pred <- calcKp_PT(2.15,9.7,0.8,BP=1.52,type=3) #prediction using reported physiology
# met_pred <- c(met_pred[1:3],met_pred[5:6],met_pred[4],met_pred[7:9],met_pred[11],met_pred[10])
# 
# # Acebutolol-R
# ace_pred <- calcKp_PT(1.4,9.7,0.79,type=3) 
# ace_pred <- c(ace_pred[1:3],ace_pred[5:6],ace_pred[4],ace_pred[7:9],ace_pred[11],ace_pred[10])
# 
# 
# ### Weak Base ###
# # Voriconazole
# vori_pred <- calcKp_PT(2.56,1.76,0.42,BP=1,type=3) 
# vori_pred <- c(vori_pred[1:3],vori_pred[5:6],vori_pred[4],vori_pred[7:9],vori_pred[11],vori_pred[10])
# 
# #Alprazolam
# alp_pred <- calcKp_PT(2.5,2.4,0.35,type=3) 
# alp_pred <- c(alp_pred[1:3],alp_pred[5:6],alp_pred[4],alp_pred[7:9],alp_pred[11],alp_pred[10])
# 
# 
# ### Acid ###
# # Thiopental
# thio_pred <- calcKp_PT(2.9,7.5,0.13, BP=1,type=2) 
# thio_pred <- c(thio_pred[1:3],thio_pred[5:6],thio_pred[4],thio_pred[7:9],thio_pred[11],thio_pred[10])
# 
# # Phenobarbital
# phen_pred <- calcKp_PT(0.8,7.35,0.78, type=2) 
# phen_pred <- c(phen_pred[1:3],phen_pred[5:6],phen_pred[4],phen_pred[7:9],phen_pred[11],phen_pred[10])
# 
# 
# ### Neutral ###
# # Digoxin
# dig_pred <- calcKp_PT(1.48,0,0.87,type=1) 
# dig_pred <- c(dig_pred[1:3],dig_pred[5:6],dig_pred[4],dig_pred[7:9],dig_pred[11],dig_pred[10])
# 
# # Ethoxybenzamide
# etho_pred <- calcKp_PT(0.77,0,0.5,BP=1,type=1) 
# etho_pred <- c(etho_pred[1:3],etho_pred[5:6],etho_pred[4],etho_pred[7:9],etho_pred[11],etho_pred[10])
# 
# 
# ### Zwitterion ###
# # Ofloxacin
# oflo_pred <- calcKp_PT(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6) 
# oflo_pred <- c(oflo_pred[1:3],oflo_pred[5:6],oflo_pred[4],oflo_pred[7:9],oflo_pred[11],oflo_pred[10])
# 
# #Enoxacin
# eno_pred <- calcKp_PT(0.1,c(6.1,8.7),0.66,BP=0.94,type=6) 
# eno_pred <- c(eno_pred[1:3],eno_pred[5:6],eno_pred[4],eno_pred[7:9],eno_pred[11],eno_pred[10])
# 
# 
# ##prepare dataframes
# met <- data.frame(uni_pred = unlist(met_uni),pred = unlist(met_pred))
# met <- mutate(met, Drug_class = "Strong base")
# 
# ace <- data.frame(uni_pred = unlist(ace_uni),pred = unlist(ace_pred))
# ace <- mutate(ace, Drug_class = "Strong base")
# 
# vori <- data.frame(uni_pred = unlist(vori_uni),pred = unlist(vori_pred))
# vori <- mutate(vori, Drug_class = "Weak base")
# 
# alp <- data.frame(uni_pred = unlist(alp_uni),pred = unlist(alp_pred))
# alp <- mutate(alp, Drug_class = "Weak base")
# 
# thio <- data.frame(uni_pred = unlist(thio_uni),pred = unlist(thio_pred))
# thio <- mutate(thio, Drug_class = "Acid")
# 
# phen <- data.frame(uni_pred = unlist(phen_uni),pred = unlist(phen_pred))
# phen <- mutate(phen, Drug_class = "Acid")
# 
# dig <- data.frame(uni_pred = unlist(dig_uni),pred = unlist(dig_pred))
# dig <- mutate(dig, Drug_class = "Neutral")
# 
# etho <- data.frame(uni_pred = unlist(etho_uni),pred = unlist(etho_pred))
# etho <- mutate(etho, Drug_class = "Neutral")
# 
# oflo <- data.frame(uni_pred = unlist(oflo_uni),pred = unlist(oflo_pred))
# oflo <- mutate(oflo, Drug_class = "Zwitterion")
# 
# eno <- data.frame(uni_pred = unlist(eno_uni),pred = unlist(eno_pred))
# eno <- mutate(eno, Drug_class = "Zwitterion")
# 
# 
# # Bind all the outputs together
# drug_all <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)
# 
# # Calculate the Pearson correlation coefficient 
# corr_coeff <- cor(drug_all$uni_pred,drug_all$pred,method='pearson') 
# 
# plot_PT <- ggplot() +
#   geom_point(data=drug_all, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
#   scale_shape_discrete(solid=F) +
#   geom_abline(intercept = 0, slope = 1) +
#   xlim(0,15) +
#   ylim(0,15) +
#   xlab("Predicted Kp (unified physiology)") +
#   ylab("Predicted Kp (reported physiology)") +
#   # ggtitle("Poulin and Theil") +
#   ggtitle("(a)") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 10))  +
#   theme(legend.justification=c(0.025,0.95), legend.position=c(0.025,0.95), legend.key=element_blank()) + theme(legend.title =  element_blank()) #+
# #theme(plot.title = element_text(hjust = 0.5))


source("CalcKp_P&T.R")
source("Copies/CalcKp_P&T_copy.R")

# uni output order: adipose, bone, brain, heart, kidney, gut, liver, lung, muscle, skin, spleen
# pred output order: adipose, bone, brain, gut, heart, kidney, liver, lung, muscle, spleen, skin

### Strong Base ###
#Metoprolol
met_uni <- calcKp_PT(2.15,9.7,0.8,BP=1.52,type=3) #prediction using unified physiology
met_pred <- calcKp_PT_dat(2.15,9.7,0.8,BP=1.52,type=3) #prediction using reported physiology
met_pred <- c(met_pred[1:3],met_pred[5:6],met_pred[4],met_pred[7:9],met_pred[11],met_pred[10])
met <- data.frame(uni_pred = unlist(met_uni),pred = unlist(met_pred))
met <- mutate(met, Drug_class = "Strong base")

# Acebutolol-R
ace_uni <- calcKp_PT(1.4,9.7,0.79,type=3)
ace_pred <- calcKp_PT_dat(1.4,9.7,0.79,type=3) 
ace_pred <- c(ace_pred[1:3],ace_pred[5:6],ace_pred[4],ace_pred[7:9],ace_pred[11],ace_pred[10])
ace <- data.frame(uni_pred = unlist(ace_uni),pred = unlist(ace_pred))
ace <- mutate(ace, Drug_class = "Strong base")


### Weak Base ###
# Voriconazole
vori_uni <- calcKp_PT(2.56,1.76,0.42,BP=1,type=3) 
vori_pred <- calcKp_PT_dat(2.56,1.76,0.42,BP=1,type=3) 
vori_pred <- c(vori_pred[1:3],vori_pred[5:6],vori_pred[4],vori_pred[7:9],vori_pred[11],vori_pred[10])
vori <- data.frame(uni_pred = unlist(vori_uni),pred = unlist(vori_pred))
vori <- mutate(vori, Drug_class = "Weak base")

#Alprazolam
alp_uni <- calcKp_PT(2.5,2.4,0.35,type=3) 
alp_pred <- calcKp_PT_dat(2.5,2.4,0.35,type=3) 
alp_pred <- c(alp_pred[1:3],alp_pred[5:6],alp_pred[4],alp_pred[7:9],alp_pred[11],alp_pred[10])
alp <- data.frame(uni_pred = unlist(alp_uni),pred = unlist(alp_pred))
alp <- mutate(alp, Drug_class = "Weak base")


### Acid ###
# Thiopental
thio_uni <- calcKp_PT(2.9,7.5,0.13, BP=1,type=2)
thio_pred <- calcKp_PT_dat(2.9,7.5,0.13, BP=1,type=2) 
thio_pred <- c(thio_pred[1:3],thio_pred[5:6],thio_pred[4],thio_pred[7:9],thio_pred[11],thio_pred[10])
thio <- data.frame(uni_pred = unlist(thio_uni),pred = unlist(thio_pred))
thio <- mutate(thio, Drug_class = "Acid")

# Phenobarbital
phen_uni <- calcKp_PT(0.8,7.35,0.78, type=2)
phen_pred <- calcKp_PT_dat(0.8,7.35,0.78, type=2) 
phen_pred <- c(phen_pred[1:3],phen_pred[5:6],phen_pred[4],phen_pred[7:9],phen_pred[11],phen_pred[10])
phen <- data.frame(uni_pred = unlist(phen_uni),pred = unlist(phen_pred))
phen <- mutate(phen, Drug_class = "Acid")


### Neutral ###
# Digoxin
dig_uni <- calcKp_PT(1.48,0,0.87,type=1)
dig_pred <- calcKp_PT_dat(1.48,0,0.87,type=1) 
dig_pred <- c(dig_pred[1:3],dig_pred[5:6],dig_pred[4],dig_pred[7:9],dig_pred[11],dig_pred[10])
dig <- data.frame(uni_pred = unlist(dig_uni),pred = unlist(dig_pred))
dig <- mutate(dig, Drug_class = "Neutral")

# Ethoxybenzamide
etho_uni <- calcKp_PT(0.77,0,0.5,BP=1,type=1)
etho_pred <- calcKp_PT_dat(0.77,0,0.5,BP=1,type=1) 
etho_pred <- c(etho_pred[1:3],etho_pred[5:6],etho_pred[4],etho_pred[7:9],etho_pred[11],etho_pred[10])
etho <- data.frame(uni_pred = unlist(etho_uni),pred = unlist(etho_pred))
etho <- mutate(etho, Drug_class = "Neutral")


### Zwitterion ###
# Ofloxacin
oflo_uni <- calcKp_PT(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6)
oflo_pred <- calcKp_PT_dat(-0.4,c(5.97,9.28),0.77,BP=0.92,type=6) 
oflo_pred <- c(oflo_pred[1:3],oflo_pred[5:6],oflo_pred[4],oflo_pred[7:9],oflo_pred[11],oflo_pred[10])
oflo <- data.frame(uni_pred = unlist(oflo_uni),pred = unlist(oflo_pred))
oflo <- mutate(oflo, Drug_class = "Zwitterion")

#Enoxacin
eno_uni <- calcKp_PT(0.1,c(6.1,8.7),0.66,BP=0.94,type=6)
eno_pred <- calcKp_PT_dat(0.1,c(6.1,8.7),0.66,BP=0.94,type=6) 
eno_pred <- c(eno_pred[1:3],eno_pred[5:6],eno_pred[4],eno_pred[7:9],eno_pred[11],eno_pred[10])
eno <- data.frame(uni_pred = unlist(eno_uni),pred = unlist(eno_pred))
eno <- mutate(eno, Drug_class = "Zwitterion")


# Bind all the outputs together
drug_all <- rbind(met,ace,vori,alp,thio,phen,dig,etho,oflo,eno)

# Calculate the Pearson correlation coefficient 
corr_coeff <- cor(drug_all$uni_pred,drug_all$pred,method='pearson')
corr_coeff <- paste("italic(r) %~~%", corr_coeff)

gp1 <- ggplot() +
  geom_point(data=drug_all, aes(x=uni_pred, y=pred, shape=Drug_class), size=3) +
  scale_shape_discrete(solid=F) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,max(drug_all$pred)) +
  ylim(0,max(drug_all$pred)) +
  ggtitle("a    Poulin and Theil") + 
  xlab("Predicted Kp (unified tissue composition)") +
  ylab("Predicted Kp (reported tissue composition)") +
  th1 +
  geom_text(aes(x=5, y=7.5, label=corr_coeff), parse=TRUE)


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 4: Figure 4 ##############################################
#########################################################################################################
## This figure shows the range of Kp/Vss predcitions for each drug class for each tissue using different 
## calculation methods


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 5: Figure 5 ##############################################
#########################################################################################################
## This figure sows the impact of calculation methods on the model predictions from each drug class 
## representatives


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 6: Figure 6 ##############################################
#########################################################################################################
## This figure shows prediction ranges for representative drugs after incorporating inter-individual 
## variability and uncertainty in physiological and drug-related parameters, respectively


#########################################################################################################
#########################################################################################################

#########################################################################################################
######################################## CHUNK 7: Table 1 ###############################################
#########################################################################################################
## This table shows the error estimates for PK parameters and RMSE from the model predictions from each
## drug class representatives


#########################################################################################################
#########################################################################################################



