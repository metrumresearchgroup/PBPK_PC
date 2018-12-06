######################################################################################################
######################################### Artemether - neutral #######################################
######################################################################################################
# Reference: https://jpharmsci.org/article/S0022-3549(16)41525-X/fulltext
# Compile model
mod <- mread("../model/artemetherPBPK_adult")

# Load data
df <- read.csv("../data/obs_artemether.csv")

# Convert output and data to common units: ng/mL -> mg/L
df <- df %>%
  mutate(conc = conc*(10^-6)*10^3) %>%
  mutate(max_sd = max_sd*(10^-6)*10^3)

# Calculate Kps
type <- 1  #neutral
logP <- 3.28
pKa <- 0  
fup <- 0.046
BP <- 0.8  

Kps1 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="P&T", dat_uni)
Kps2 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Berez", dat_uni)
Kps3 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="R&R", dat_uni)
Kps4 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Schmitt", dat_uni)
Kps5 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="pksim", dat_uni)

# ICRP parameters for men
men <- list(BW = 73, Vad = 18.2, Vbo = 10.5, Vbl = 5.6, Vlu = 0.5, Vbr = 1.45, 
            Vhe = 0.33, Vki = 0.31, Vmu = 29, Vsk = 3.3, Vli = 1.8, Vsp = 0.15, Vpa = 0.14, 
            Vgu = 1.3, Qc = 390, Qad = 19.5, Qbo = 19.5, Qbr = 46.8, Qhe = 15.6, Qki = 74.1, 
            Qmu = 66.3, Qsk = 19.5, Qli = 99.45, Qsp = 11.7, Qpa = 3.9, Qgu = 58.5)

# New parameter list
pars1 <- c(Kps1, men)
pars2 <- c(Kps2, men)
pars3 <- c(Kps3, men)
pars4 <- c(Kps4, men)
pars5 <- c(Kps5, men)

# 2 doses of 80 mg each, the second at t = 8
e <- ev(amt = 80, ii = 8, addl = 1, cmt= "D")

outFun <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         ev(e) %>%
                         param(pars) %>%
                         mrgsim(end = 16, delta = 0.01) %>%
                         #mutate(Cplasma = Cplasma * 1000) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 2))
  return(out)
}

out1 <- outFun(pars1) %>%
  mutate(Method = "Poulin and Theil")
out2 <- outFun(pars2) %>%
  mutate(Method = "Berezhkovskiy")
out3 <- outFun(pars3) %>%
  mutate(Method = "Rodgers and Rowland")
out4 <- outFun(pars4) %>%
  mutate(Method = "Schmitt")
out5 <- outFun(pars5) %>%
  mutate(Method = "PK-Sim Standard")

# Bind all outputs into one matrix
out_all <- rbind(out1,out2,out3,out4,out5,out5)

# Select outputs every num_pts rows to add the geom_point shapes layer over the geom_line
num_pts <- 60  # num_pts indicates the number of rows to skip, this influences the density of the points in the plot
out1_reduced <- out1[seq(1, nrow(out1), num_pts), ]
out2_reduced <- out2[seq(1, nrow(out2), num_pts), ]
out3_reduced <- out3[seq(1, nrow(out3), num_pts), ]
out4_reduced <- out4[seq(1, nrow(out4), num_pts), ]
out5_reduced <- out5[seq(1, nrow(out5), num_pts), ]
out_reduced <- rbind(out1_reduced,out2_reduced,out3_reduced,out4_reduced,out5_reduced)


# Time points from observed data
sample <- df$time

outFun_new <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end=-1, add=sample) %>%
                         # mutate(Cplasma = Cplasma * 1000) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 1))
  return(out)
}
out1_new <- outFun_new(pars1)
out2_new <- outFun_new(pars2)
out3_new <- outFun_new(pars3)
out4_new <- outFun_new(pars4)
out5_new <- outFun_new(pars5)

#---------------
# Calculate the relative root mean square error for each method
rmse <- function(pred,obs){
  rmsd <- sqrt(mean((pred - obs)^2))/(max(obs)-min(obs))
  return(rmsd)
}

rmse_PT <- rmse(out1_new$Cplasma,df$conc)
rmse_Berez <- rmse(out2_new$Cplasma,df$conc)
rmse_RR <- rmse(out3_new$Cplasma,df$conc)
rmse_Schmitt <- rmse(out4_new$Cplasma,df$conc)
rmse_pksim <- rmse(out5_new$Cplasma,df$conc)

rel_rmse <- c(rmse_PT,rmse_Berez,rmse_RR,rmse_Schmitt,rmse_pksim)
rel_rmse

#---------------------
# Calculate AUC for each curve

auc_obs <- round(pk.calc.auc(df$conc,sample, interval=c(sample[1],last(sample))), digits=3)
auc_PT <- pk.calc.auc(out1$Cplasma,out1$time,interval=c(sample[1],last(sample)))
auc_Berez <-pk.calc.auc(out2$Cplasma,out2$time,interval=c(sample[1],last(sample)))
auc_RR <- pk.calc.auc(out3$Cplasma,out3$time,interval=c(sample[1],last(sample)))
auc_Schmitt <- pk.calc.auc(out4$Cplasma,out4$time,interval=c(sample[1],last(sample)))
auc_pksim <- pk.calc.auc(out5$Cplasma,out5$time,interval=c(sample[1],last(sample)))

auc_all <- round(c(auc_PT,auc_Berez,auc_RR,auc_Schmitt,auc_pksim), digits=3)

auc_error <- round(auc_obs - auc_all, digits=3)

# # Combine in data frame
auc_art <- round(cbind(auc_obs,auc_all,auc_error), digits=3)


#---------------------
# Calculate half life for each curve
hl_obs <- pk.calc.half.life(df$conc,sample)
hl_PT <- pk.calc.half.life(out1$Cplasma,out1$time)
hl_Berez <- pk.calc.half.life(out2$Cplasma,out2$time)
hl_RR <- pk.calc.half.life(out3$Cplasma,out3$time)
hl_Schmitt <- pk.calc.half.life(out4$Cplasma,out4$time)
hl_pksim <- pk.calc.half.life(out5$Cplasma,out5$time)


hl_all <- round(c(hl_PT$half.life,hl_Berez$half.life,hl_RR$half.life,hl_Schmitt$half.life,hl_pksim$half.life), digits=3)

hl_error <- round(hl_obs$half.life - hl_all, digits=3)

hl_art <- cbind(hl_obs$half.life,hl_all,hl_error)


#---------------------
# Combine PK info in data frame and store as csv file
pk_art <- cbind(rel_rmse,auc_art,hl_art)
colnames(pk_art) <- c("RelRMSE","AUCobs","AUCpred", "AUCerror", "hlobs", "hlpred", "hlerror")
pk_art <- mutate(as.data.frame(pk_art), Method=c("Poulin and Theil", "Berezhkovskiy", "Rodgers and Rowland", "Schmitt", "PK-Sim Standard"))

pk_art <- pk_art %>% mutate(Type="Neutral")
