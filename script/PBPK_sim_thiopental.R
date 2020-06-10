####################################################################################################
######################################### Thiopental - acid ########################################
####################################################################################################
# Reference: https://insights.ovid.com/pubmed?pmid=8780280
# Drug given was a racemic mixture: R-thiopental and S-thiopental measured and modeled separately 
# Compile model
mod <- mread("../model/thiopentalPBPK_Adult")

# Load observed data
df <- read.csv("../data/obs_thiopental") #S-thiopental

getThiPK <- function(unified = T, config){

#df <- df[order(df$time),] # re-order the data by time
CLhepatic = 13.8 #S-thiopental

# Calculate Kps
type <- 2  #acid
logP <- 2.9  
pKa <- 7.5  
fup <- 0.13
BP <- 1  


if(unified){
  Kps1 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="P&T", dat_uni)
  Kps2 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Berez", dat_uni)
  Kps3 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="R&R", dat_uni)
  Kps4 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Schmitt", dat_uni)
  Kps5 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="pksim", dat_uni) 
}else{
  Kps1 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="P&T", dat_PT)
  Kps2 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Berez", dat_Berez)
  Kps3 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="R&R", dat_RR)
  Kps4 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Schmitt", dat_Schmitt_rep)
  Kps5 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="pksim", dat_pksim)
}

# New parameter list
pars1 <- c(Kps1, list(CLhepatic=CLhepatic))
pars2 <- c(Kps2, list(CLhepatic=CLhepatic))
pars3 <- c(Kps3, list(CLhepatic=CLhepatic))
pars4 <- c(Kps4, list(CLhepatic=CLhepatic))
pars5 <- c(Kps5, list(CLhepatic=CLhepatic))

# Simulate 250 mg IV dose
cmt <- "VEN"
dose = 250
e <- ev(amt=dose, cmt=cmt) #set event

outFun <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end=15, delta=0.01) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 2))
  return(out)
}

out1 <- outFun(pars1) %>%
  mutate(Method = "PT")
out2 <- outFun(pars2) %>%
  mutate(Method = "Berez")
out3 <- outFun(pars3) %>%
  mutate(Method = "RR")
out4 <- outFun(pars4) %>%
  mutate(Method = "Schmitt")
out5 <- outFun(pars5) %>%
  mutate(Method = "PK-Sim")

# Bind all outputs into one matrix
out_all_S <- rbind(out1,out2,out3,out4,out5)

# Select outputs every num_pts rows to add the geom_point shapes layer over the geom_line
num_pts <- 110
out1_reduced <- out1[seq(1, nrow(out1), num_pts), ]
out2_reduced <- out2[seq(1, nrow(out2), num_pts), ]
out3_reduced <- out3[seq(1, nrow(out3), num_pts), ]
out4_reduced <- out4[seq(1, nrow(out4), num_pts), ]
out5_reduced <- out5[seq(1, nrow(out5), num_pts), ]
out_reduced_S <- rbind(out1_reduced,out2_reduced,out3_reduced,out4_reduced,out5_reduced)




######################################### R-thiopental - acid ########################################
######################################################################################################

# Load observed data (R-thiopental)
df_R <- read.csv("../data/obs_thiopental_R.csv") 

CLhepatic = 17.7 #R-thiopental

# New parameter list
pars1_R <- c(Kps1, list(CLhepatic=CLhepatic))
pars2_R <- c(Kps2, list(CLhepatic=CLhepatic))
pars3_R <- c(Kps3, list(CLhepatic=CLhepatic))
pars4_R <- c(Kps4, list(CLhepatic=CLhepatic))
pars5_R <- c(Kps5, list(CLhepatic=CLhepatic))

out1_R <- outFun(pars1_R) %>%
  mutate(Method = "PT")
out2_R <- outFun(pars2_R) %>%
  mutate(Method = "Berez")
out3_R <- outFun(pars3_R) %>%
  mutate(Method = "RR")
out4_R <- outFun(pars4_R) %>%
  mutate(Method = "Schmitt")
out5_R <- outFun(pars5_R) %>%
  mutate(Method = "PK-Sim")

# Bind all outputs into one matrix
out_all_R <- rbind(out1_R,out2_R,out3_R,out4_R,out5_R)

# Select outputs every num_pts rows to add the geom_point shapes layer over the geom_line
num_pts <- 110
out1_reduced_R <- out1_R[seq(1, nrow(out1_R), num_pts), ]
out2_reduced_R <- out2_R[seq(1, nrow(out2_R), num_pts), ]
out3_reduced_R <- out3_R[seq(1, nrow(out3_R), num_pts), ]
out4_reduced_R <- out4_R[seq(1, nrow(out4_R), num_pts), ]
out5_reduced_R <- out5_R[seq(1, nrow(out5_R), num_pts), ]
out_reduced_R <- rbind(out1_reduced_R,out2_reduced_R,out3_reduced_R,out4_reduced_R,out5_reduced_R)


# Calculate errors for S-thiopental
# Time points from observed data
sample <- df$time 

#e <- ev(amt=dose, cmt=cmt, rate=dose) #set event

outFun_new <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end=-1, add=sample) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 1))
  return(out)
}
out1_new <- outFun_new(pars1)
out2_new <- outFun_new(pars2)
out3_new <- outFun_new(pars3)
out4_new <- outFun_new(pars4)
out5_new <- outFun_new(pars5)


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

rel_rmse <- c(rmse_PT,rmse_Berez,rmse_RR,rmse_Schmitt,rmse_pksim)*100

#---------------------
# Calculate AUC for each curve
auc_obs <- pk.calc.auc(df$conc,df$time, interval=c(df$time[1],last(df$time)))
auc_PT <- pk.calc.auc(out1$Cplasma,out1$time,interval=c(df$time[1],last(df$time)))
auc_Berez <-pk.calc.auc(out2$Cplasma,out2$time,interval=c(df$time[1],last(df$time)))
auc_RR <- pk.calc.auc(out3$Cplasma,out3$time,interval=c(df$time[1],last(df$time)))
auc_Schmitt <- pk.calc.auc(out4$Cplasma,out4$time,interval=c(df$time[1],last(df$time)))
auc_pksim <- pk.calc.auc(out5$Cplasma,out5$time,interval=c(df$time[1],last(df$time)))

auc_all <- c(auc_PT,auc_Berez,auc_RR,auc_Schmitt,auc_pksim)

# Error
#auc_error <- auc_obs - auc_all

# Relative error
auc_error <- (abs(auc_obs - auc_all)/auc_obs)*100

# # Combine in data frame
auc_thio <- cbind(auc_obs,auc_all,auc_error)


#---------------------
# Calculate half life for each curve
hl_obs <- pk.calc.half.life(df$conc,df$time)
hl_PT <- pk.calc.half.life(out1$Cplasma,out1$time)
hl_Berez <- pk.calc.half.life(out2$Cplasma,out2$time)
hl_RR <- pk.calc.half.life(out3$Cplasma,out3$time)
hl_Schmitt <- pk.calc.half.life(out4$Cplasma,out4$time)
hl_pksim <- pk.calc.half.life(out5$Cplasma,out5$time)


hl_all <- c(hl_PT$half.life,hl_Berez$half.life,hl_RR$half.life,hl_Schmitt$half.life,hl_pksim$half.life)

# Error
#hl_error <- hl_obs$half.life - hl_all

# Relative error
hl_error <- (abs(hl_obs$half.life - hl_all)/hl_obs$half.life)*100

pk_thio <- cbind(rel_rmse,auc_thio,hl_obs$half.life,hl_all,hl_error)
colnames(pk_thio) <- c("RelRMSE","AUCobs","AUCpred", "AUCerror", "hlobs", "hlpred", "hlerror")
pk_thio <- mutate(as.data.frame(pk_thio), Method=c("PT", "Berez", "RR", "Schmitt", "PK-Sim"))

# Store the PK info
pk_thio_S <- pk_thio %>% mutate(Type="Acid")


# Calculate erros for R-thiopental
out1_new_R <- outFun_new(pars1_R)
out2_new_R <- outFun_new(pars2_R)
out3_new_R <- outFun_new(pars3_R)
out4_new_R <- outFun_new(pars4_R)
out5_new_R <- outFun_new(pars5_R)

rmse_PT <- rmse(out1_new_R$Cplasma,df_R$conc)
rmse_Berez <- rmse(out2_new_R$Cplasma,df_R$conc)
rmse_RR <- rmse(out3_new_R$Cplasma,df_R$conc)
rmse_Schmitt <- rmse(out4_new_R$Cplasma,df_R$conc)
rmse_pksim <- rmse(out5_new_R$Cplasma,df_R$conc)

rel_rmse <- c(rmse_PT,rmse_Berez,rmse_RR,rmse_Schmitt,rmse_pksim)*100

#---------------------
# Calculate AUC for each curve
auc_obs <- pk.calc.auc(df_R$conc,df_R$time, interval=c(df_R$time[1],last(df_R$time)))
auc_PT <- pk.calc.auc(out1_R$Cplasma,out1$time,interval=c(df_R$time[1],last(df_R$time)))
auc_Berez <-pk.calc.auc(out2_R$Cplasma,out2$time,interval=c(df_R$time[1],last(df_R$time)))
auc_RR <- pk.calc.auc(out3_R$Cplasma,out3$time,interval=c(df_R$time[1],last(df_R$time)))
auc_Schmitt <- pk.calc.auc(out4_R$Cplasma,out4$time,interval=c(df_R$time[1],last(df_R$time)))
auc_pksim <- pk.calc.auc(out5_R$Cplasma,out5$time,interval=c(df_R$time[1],last(df_R$time)))

auc_all <- c(auc_PT,auc_Berez,auc_RR,auc_Schmitt,auc_pksim)

# Error
#auc_error <- auc_obs - auc_all

# Relative error
auc_error <- (abs(auc_obs - auc_all)/auc_obs)*100

# # Combine in data frame
auc_thio <- cbind(auc_obs,auc_all,auc_error)


#---------------------
# Calculate half life for each curve
hl_obs <- pk.calc.half.life(df_R$conc,df_R$time)
hl_PT <- pk.calc.half.life(out1_R$Cplasma,out1$time)
hl_Berez <- pk.calc.half.life(out2_R$Cplasma,out2$time)
hl_RR <- pk.calc.half.life(out3_R$Cplasma,out3$time)
hl_Schmitt <- pk.calc.half.life(out4_R$Cplasma,out4$time)
hl_pksim <- pk.calc.half.life(out5_R$Cplasma,out5$time)


hl_all <- c(hl_PT$half.life,hl_Berez$half.life,hl_RR$half.life,hl_Schmitt$half.life,hl_pksim$half.life)

# Error
#hl_error <- hl_obs$half.life - hl_all

# Relative error
hl_error <- (abs(hl_obs$half.life - hl_all)/hl_obs$half.life)*100

pk_thio <- cbind(rel_rmse,auc_thio,hl_obs$half.life,hl_all,hl_error)
colnames(pk_thio) <- c("RelRMSE","AUCobs","AUCpred", "AUCerror", "hlobs", "hlpred", "hlerror")
pk_thio <- mutate(as.data.frame(pk_thio), Method=c("PT", "Berez", "RR", "Schmitt", "PK-Sim"))

# Store the PK info
pk_thio_R <- pk_thio %>% mutate(Type="Acid")

if(config == "R"){
  return(pk_thio_R) 
}else{
  return(pk_thio_S)
}

}

