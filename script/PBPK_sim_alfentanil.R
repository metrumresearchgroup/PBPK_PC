######################################### Alfentanil - weak base ########################################

# Compile model
mod <- mread("../model/alfentanilPBPK_Adult")

# Load data
df <- read.csv("../data/obs_alfentanil.csv")

getAlfPK <- function(unified = T){
  # Reference: http://anesthesiology.pubs.asahq.org/article.aspx?articleid=1948182
  
  # Calculate mean and SD for each time point
  df <- df %>%
    group_by(time) %>%
    mutate(avg=mean(conc)) %>%
    mutate(stand_dev=sd(conc))
  df_mean <- cbind(df$time[!duplicated(df$time)],df$avg[!duplicated(df$avg)],df$stand_dev[!duplicated(df$stand_dev)])
  colnames(df_mean) <- c("time","mean","SD")
  df_mean <- as.data.frame(df_mean)
  
  # Convert output and data to common units: ng/mL -> mg/L
  df_mean <- df_mean %>%
    mutate(mean = mean*(10^-6)*10^3) %>%
    mutate(SD = SD*(10^-6)*10^3)
  
  # Calculate Kps
  type <- 3  #base
  logP <- 2.2
  pKa <- 6.5  
  fup <- 0.11
  BP <- 0.63  #http://anesthesiology.pubs.asahq.org/article.aspx?articleid=1948182
  
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
  
  cmt <- "VEN"
  dose = 50*70/1000  #http://anesthesiology.pubs.asahq.org/article.aspx?articleid=1948182
  e <- ev(amt=dose, cmt=cmt) #set event
  outFun <- function(pars){
    out <- as.data.frame(out <- mod %>%
                           param(pars) %>%
                           ev(e) %>%
                           mrgsim(end=12, delta=0.01) %>%
                           mutate(Cplasma = Cplasma * 10^-3) %>% # Convert units: ng/mL -> mg/L
                           as.data.frame() %>%
                           dplyr::filter(dplyr::row_number() > 2))
    return(out)
  }
  
  out1 <- outFun(Kps1) %>%
    mutate(Method = "PT")
  out2 <- outFun(Kps2) %>%
    mutate(Method = "Berez")
  out3 <- outFun(Kps3) %>%
    mutate(Method = "RR")
  out4 <- outFun(Kps4) %>%
    mutate(Method = "Schmitt")
  out5 <- outFun(Kps5) %>%
    mutate(Method = "PK-Sim")
  
  # Bind all outputs into one matrix
  out_all <- rbind(out1,out2,out3,out4,out5)
  
  # Select outputs every num_pts rows to add the geom_point shapes layer over the geom_line
  num_pts <- 70
  out1_reduced <- out1[seq(1, nrow(out1), num_pts), ]
  out2_reduced <- out2[seq(1, nrow(out2), num_pts), ]
  out3_reduced <- out3[seq(1, nrow(out3), num_pts), ]
  out4_reduced <- out4[seq(1, nrow(out4), num_pts), ]
  out5_reduced <- out5[seq(1, nrow(out5), num_pts), ]
  out_reduced <- rbind(out1_reduced,out2_reduced,out3_reduced,out4_reduced,out5_reduced)
  
  
  # Simulate only at the time points in the data set
  
  # Time points from observed data
  sample <- df_mean$time
  
  
  outFun_new <- function(pars){
    out <- as.data.frame(out <- mod %>%
                           param(pars) %>%
                           ev(e) %>%
                           mrgsim(end=-1, add=sample) %>%
                           mutate(Cplasma = Cplasma * 10^-3) %>% # Convert units: ng/mL -> mg/L
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
  # Calculate the relative root mean square error for each method
  
  rmse <- function(pred,obs){
    rmsd <- sqrt(mean((pred - obs)^2))/(max(obs)-min(obs))
    return(rmsd)
  }
  
  rmse_PT <- rmse(out1_new$Cplasma,df_mean$mean)
  rmse_Berez <- rmse(out2_new$Cplasma,df_mean$mean)
  rmse_RR <- rmse(out3_new$Cplasma,df_mean$mean)
  rmse_Schmitt <- rmse(out4_new$Cplasma,df_mean$mean)
  rmse_pksim <- rmse(out5_new$Cplasma,df_mean$mean)
  
  rel_rmse <- c(rmse_PT,rmse_Berez,rmse_RR,rmse_Schmitt,rmse_pksim)*100
  
  #---------------------
  # Calculate AUC for each curve
  
  auc_obs <- pk.calc.auc(df_mean$mean,df_mean$time, interval=c(sample[1],last(sample)))
  auc_PT <- pk.calc.auc(out1$Cplasma,out1$time,interval=c(sample[1],last(sample)))
  auc_Berez <- pk.calc.auc(out2$Cplasma,out2$time,interval=c(sample[1],last(sample)))
  auc_RR <- pk.calc.auc(out3$Cplasma,out3$time,interval=c(sample[1],last(sample)))
  auc_Schmitt <- pk.calc.auc(out4$Cplasma,out4$time,interval=c(sample[1],last(sample)))
  auc_pksim <- pk.calc.auc(out5$Cplasma,out5$time,interval=c(sample[1],last(sample)))
  
  auc_all <- c(auc_PT,auc_Berez,auc_RR,auc_Schmitt,auc_pksim)
  
  # Error
  #auc_error <- auc_obs - auc_all
  
  # Relative error
  auc_error <- (abs(auc_obs - auc_all)/auc_obs)*100
  
  # # Combine in data frame
  auc_alf <- cbind(auc_obs,auc_all,auc_error)
  
  
  #---------------------
  # Calculate half life for each curvedf_2$time
  hl_obs <- pk.calc.half.life(df_mean$mean,df_mean$time)
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
  
  pk_alf <- cbind(rel_rmse,auc_alf,hl_obs$half.life,hl_all,hl_error)
  colnames(pk_alf) <- c("RelRMSE","AUCobs","AUCpred", "AUCerror", "hlobs", "hlpred", "hlerror")
  pk_alf <- mutate(as.data.frame(pk_alf), Method=c("PT", "Berez", "RR", "Schmitt", "PK-Sim"))
  
  
  # Store the PK info
  pk_alf <- pk_alf %>% mutate(Type="Weak base")
  
  return(pk_alf)
  
}

