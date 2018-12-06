######################################################################################################
########################################## Caffeine - strong base ####################################
######################################################################################################
# Reference: https://bpspubs.onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2125.2012.04363.x
# Compile model
mod <- mread("../model/caffeinePBPK_adult")

# Load data
df <- read.csv("../data/obs_caffeine.csv")

# Calculate Kps
type <- 3 #base  
logP <- -0.07
pKa <- 10.4
fup <- 0.681
BP <- 0.98

Kps1 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="P&T", dat_uni)
Kps2 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Berez", dat_uni)
Kps3 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="R&R", dat_uni)
Kps4 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="Schmitt", dat_uni)
Kps5 <- pcoeffs(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, pred="pksim", dat_uni)

# Simulate P.O dose
cmt <- "D"
dose = 150
e <- ev(amt=dose, cmt=cmt) #set event

outFun <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end = 24, delta = 0.01) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 2))
  return(out)
}

out1 <- outFun(Kps1) %>%
  mutate(Method = "Poulin and Theil")
out2 <- outFun(Kps2) %>%
  mutate(Method = "Berezhkovskiy")
out3 <- outFun(Kps3) %>%
  mutate(Method = "Rodgers and Rowland")
out4 <- outFun(Kps4) %>%
  mutate(Method = "Schmitt")
out5 <- outFun(Kps5) %>%
  mutate(Method = "PK-Sim Standard")

# Bind all outputs into one matrix
out_all <- rbind(out1,out2,out3,out4,out5)

# Select outputs every num_pts rows to add the geom_point shapes layer over the geom_line
num_pts <- 110 # num_pts indicates the number of rows to skip, this influences the density of the points in the plot
out1_reduced <- out1[seq(1, nrow(out1), num_pts), ]
out2_reduced <- out2[seq(1, nrow(out2), num_pts), ]
out3_reduced <- out3[seq(1, nrow(out3), num_pts), ]
out4_reduced <- out4[seq(1, nrow(out4), num_pts), ]
out5_reduced <- out5[seq(1, nrow(out5), num_pts), ]
out_reduced <- rbind(out1_reduced,out2_reduced,out3_reduced,out4_reduced,out5_reduced)



# Simulate only at the time points in the data set

# Remove t=0 from data set
df <- df %>%
  dplyr::filter(dplyr::row_number() > 1)

# Time points from observed data
sample <- df$time


outFun_new <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end=-1, add=sample) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 1)) # use 2 to filter out first data point (very close to t=0)
  return(out)
}

out1_new <- outFun_new(Kps1)
out2_new <- outFun_new(Kps2)
out3_new <- outFun_new(Kps3)
out4_new <- outFun_new(Kps4)
out5_new <- outFun_new(Kps5)


#---------------
#Calculate the relative root mean square error for each method

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

auc_obs <- pk.calc.auc(df$conc,df$time, interval=c(sample[2],last(sample)))
auc_PT <- pk.calc.auc(out1$Cplasma,out1$time,interval=c(sample[2],last(sample)))
auc_Berez <- pk.calc.auc(out2$Cplasma,out2$time,interval=c(sample[2],last(sample)))
auc_RR <- pk.calc.auc(out3$Cplasma,out3$time,interval=c(sample[2],last(sample)))
auc_Schmitt <- pk.calc.auc(out4$Cplasma,out4$time,interval=c(sample[2],last(sample)))
auc_pksim <- pk.calc.auc(out5$Cplasma,out5$time,interval=c(sample[2],last(sample)))

auc_all <- c(auc_PT,auc_Berez,auc_RR,auc_Schmitt,auc_pksim)
auc_all

auc_error <- auc_obs - auc_all

# # Combine in data frame
auc_caf <- cbind(auc_obs,auc_all,auc_error)


#---------------------
# Calculate half life for each curvedf_2$time
hl_obs <- pk.calc.half.life(df$conc,df$time)
hl_PT <- pk.calc.half.life(out1$Cplasma,out1$time)
hl_Berez <- pk.calc.half.life(out2$Cplasma,out2$time)
hl_RR <- pk.calc.half.life(out3$Cplasma,out3$time)
hl_Schmitt <- pk.calc.half.life(out4$Cplasma,out4$time)
hl_pksim <- pk.calc.half.life(out5$Cplasma,out5$time)

hl_all <- c(hl_PT$half.life,hl_Berez$half.life,hl_RR$half.life,hl_Schmitt$half.life,hl_pksim$half.life)

hl_error <- hl_obs$half.life - hl_all
hl_caf <- cbind(hl_obs$half.life,hl_all,hl_error)


#---------------------
# Combine PK info in data frame and store as csv file
pk_caf <- cbind(rel_rmse,auc_caf,hl_caf)
colnames(pk_caf) <- c("RelRMSE","AUCobs","AUCpred", "AUCerror", "hlobs", "hlpred", "hlerror")
pk_caf <- mutate(as.data.frame(pk_caf), Method=c("Poulin and Theil", "Berezhkovskiy", "Rodgers and Rowland", "Schmitt", "PK-Sim Standard"))

pk_caf <- pk_caf %>% mutate(Type="Strong base")

