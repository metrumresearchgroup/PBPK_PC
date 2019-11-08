######################################### Voriconazole - weak base ###################################
# Reference: 
# Compile model
mod <- mread("../model/voriPBPK_Adult")

# Load data (digitized observations and then calculated mean and sd)
df <- read.csv("../data/obs_voriconazole.csv")
sample <- df$time 

# Unified tissue composition with coefficient of variation
dat_uni_cv <- read.csv("../data/unified_tissue_comp_cv.csv")  
dat_uni <- read.csv("../data/unified_tissue_comp.csv")

# Number of simulated patients
num_patients <- 1000 

# Drug-related parameters (will sample for logP, pKa, fup, and BP)
type <- 3  # base
logP <- 2.56
pKa <- 1.76  
fup <- 0.42
BP <- 1

bw <- 73
e <- ev(amt = 4*bw, cmt = "VEN", ii = 12, addl = 13, rate = 4*bw, ss = 1)
outFun <- function(pars){
  out <- as.data.frame(out <- mod %>%
                         param(pars) %>%
                         ev(e) %>%
                         mrgsim(end=10, delta=0.01) %>%
                         as.data.frame() %>%
                         dplyr::filter(dplyr::row_number() > 1))
  return(out)
}


dat_uni_samp <- dat_uni # replace values in the unified tissue composition data frame with sampled values

# Create empty data frames to store simulation results
out_PT <- data.frame(ID=double(),time=double(),LUNG=double(),ADIPOSE=double(),BONE=double(),BRAIN=double(),
                     HEART=double(),KIDNEY=double(),MUSCLE=double(),SKIN=double(),LIVER=double(),SPLEEN=double(),
                     GUT=double(),ART=double(),VEN=double(),D=double(),PANCREAS=double(),REST=double(),
                     Cplasma=double(),Cvenous=double(),Method=character())
out_Berez <- out_PT
out_RR <- out_PT  
out_Schmitt <- out_PT
out_pksim <- out_PT

# Store PK parameters (AUC and Cmax)
PK_PT <- matrix(0,num_patients, 2)
PK_Berez <- PK_PT
PK_RR <- PK_PT
PK_Schmitt <- PK_PT
PK_pksim <- PK_PT

for (j in 1:num_patients){
  
  # Sample tissue composition values from a truncated normal distribution with mean given in the unified_tissue_comp.csv file
  # and standard deviation calculated from the coefficient of variation reported by Ruark et al. (2014)
  # Truncated normal distribution: min = mean - 2*stdev, max = mean + 2*stdev
  # Tissue composition not in Ruark et al. set sd = 5% of reported reported value
  
  suppressWarnings( # suppress warnings for the columns of tissue composition that include NA
   for (i in 1:13){
     dat_uni_samp[i,2] <- rtnorm(1, mean=dat_uni_cv$f_water[i], sd=dat_uni_cv$f_water_cv[i]*dat_uni_cv$f_water[i]/100, 
                                 lower=max(0,dat_uni_cv$f_water[i]-2*dat_uni_cv$f_water_cv[i]*dat_uni_cv$f_water[i]/100), 
                                 upper=dat_uni_cv$f_water[i]+2*dat_uni_cv$f_water_cv[i]*dat_uni_cv$f_water[i]/100)
     dat_uni_samp[i,3] <- rtnorm(1, mean=dat_uni_cv$f_lipids[i], sd=dat_uni_cv$f_lipids_cv[i]*dat_uni_cv$f_lipids[i]/100,
                                 lower=max(0,dat_uni_cv$f_lipids[i]-2*dat_uni_cv$f_lipids_cv[i]*dat_uni_cv$f_lipids[i]/100),
                                 upper=dat_uni_cv$f_lipids[i]+2*dat_uni_cv$f_lipids_cv[i]*dat_uni_cv$f_lipids[i]/100)
     dat_uni_samp[i,4] <- rtnorm(1, mean=dat_uni_cv$f_proteins[i], sd=dat_uni_cv$f_proteins_cv[i]*dat_uni_cv$f_proteins[i]/100,
                                 lower=max(0,dat_uni_cv$f_proteins[i]-2*dat_uni_cv$f_proteins_cv[i]*dat_uni_cv$f_proteins[i]/100),
                                 upper=dat_uni_cv$f_proteins[i]+2*dat_uni_cv$f_proteins_cv[i]*dat_uni_cv$f_proteins[i]/100)
     dat_uni_samp[i,6] <- rtnorm(1, mean=dat_uni_cv$f_n_l[i], sd=dat_uni_cv$f_n_l_cv[i]*dat_uni_cv$f_n_l[i]/100,
                                 lower=max(0,dat_uni_cv$f_n_l[i]-2*dat_uni_cv$f_n_l_cv[i]*dat_uni_cv$f_n_l[i]/100),
                                 upper=dat_uni_cv$f_n_l[i]+2*dat_uni_cv$f_n_l_cv[i]*dat_uni_cv$f_n_l[i]/100)
     dat_uni_samp[i,7] <- rtnorm(1, mean=dat_uni_cv$f_n_pl[i], sd=dat_uni_cv$f_n_pl_cv[i]*dat_uni_cv$f_n_pl[i]/100,
                                 lower=max(0,dat_uni_cv$f_n_pl[i]-2*dat_uni_cv$f_n_pl_cv[i]*dat_uni_cv$f_n_pl[i]/100),
                                 upper=dat_uni_cv$f_n_pl[i]+2*dat_uni_cv$f_n_pl_cv[i]*dat_uni_cv$f_n_pl[i]/100)
     dat_uni_samp[i,8] <- rtnorm(1, mean=dat_uni_cv$f_a_pl[i], sd=dat_uni_cv$f_a_pl_cv[i]*dat_uni_cv$f_a_pl[i]/100,
                                 lower=max(0,dat_uni_cv$f_a_pl[i]-2*dat_uni_cv$f_a_pl_cv[i]*dat_uni_cv$f_a_pl[i]/100),
                                 upper=dat_uni_cv$f_a_pl[i]+2*dat_uni_cv$f_a_pl_cv[i]*dat_uni_cv$f_a_pl[i]/100)
     dat_uni_samp[i,9] <- rtnorm(1, mean=dat_uni_cv$pH[i], sd=dat_uni_cv$pH_cv[i]*dat_uni_cv$pH[i]/100,
                                 lower=max(0,dat_uni_cv$pH[i]-2*dat_uni_cv$pH_cv[i]*dat_uni_cv$pH[i]/100),
                                 upper=dat_uni_cv$pH[i]+2*dat_uni_cv$pH_cv[i]*dat_uni_cv$pH[i]/100)
     
     dat_uni_samp[i,5] <- rtnorm(1, mean=dat_uni_cv$f_pl[i], sd=0.05*dat_uni_cv$f_pl[i],
                                 lower=max(0,dat_uni_cv$f_pl[i]-2*0.05*dat_uni_cv$f_pl[i]),
                                 upper=dat_uni_cv$f_pl[i]+2*0.05*dat_uni_cv$f_pl[i])
     dat_uni_samp[i,10] <- rtnorm(1, mean=dat_uni_cv$f_ew[i], sd=0.05*dat_uni_cv$f_ew[i],
                                 lower=max(0,dat_uni_cv$f_ew[i]-2*0.05*dat_uni_cv$f_ew[i]),
                                 upper=dat_uni_cv$f_ew[i]+2*0.05*dat_uni_cv$f_ew[i])
     dat_uni_samp[i,11] <- rtnorm(1, mean=dat_uni_cv$f_iw[i], sd=0.05*dat_uni_cv$f_iw[i],
                                  lower=max(0,dat_uni_cv$f_iw[i]-2*0.05*dat_uni_cv$f_iw[i]),
                                  upper=dat_uni_cv$f_iw[i]+2*0.05*dat_uni_cv$f_iw[i])
     dat_uni_samp[i,12] <- rtnorm(1, mean=dat_uni_cv$AR[i], sd=0.05*dat_uni_cv$AR[i],
                                  lower=max(0,dat_uni_cv$AR[i]-2*0.05*dat_uni_cv$AR[i]),
                                  upper=dat_uni_cv$AR[i]+2*0.05*dat_uni_cv$AR[i])
     dat_uni_samp[i,13] <- rtnorm(1, mean=dat_uni_cv$LR[i], sd=0.05*dat_uni_cv$LR[i],
                                  lower=max(0,dat_uni_cv$LR[i]-2*0.05*dat_uni_cv$LR[i]),
                                  upper=dat_uni_cv$LR[i]+2*0.05*dat_uni_cv$LR[i])
   }
  )
   
  # Sample drug parameters from a uniform distribution with 100% increase from reported value as max and 
  # 50% decrease from reported value as min.
  # Modify to handle negative logP values
  if (logP<0){
    logP_samp <- runif(1, min=logP+logP/2, max=logP/2)
  }else{
    logP_samp <- runif(1, min=logP/2, max=logP+logP/2)
  }
  pKa_samp <- runif(1, min=pKa/2, max=pKa+pKa/2) # restrict to be non-negative, note: this will need to be adjusted if simulating zwitterions
  fup_samp <- runif(1, min=fup/2, max=fup+fup/2)
  BP_samp <- runif(1, min=BP/2, max=BP+BP/2)
  
  # Calculate Kp values for sampled parameters
  Kps1 <- pcoeffs(logP=logP_samp, pKa=pKa_samp, fup=fup_samp, BP=BP_samp, type=type, pred="P&T", dat_uni_samp)
  Kps2 <- pcoeffs(logP=logP_samp, pKa=pKa_samp, fup=fup_samp, BP=BP_samp, type=type, pred="Berez", dat_uni_samp)
  Kps3 <- pcoeffs(logP=logP_samp, pKa=pKa_samp, fup=fup_samp, BP=BP_samp, type=type, pred="R&R", dat_uni_samp)
  Kps4 <- pcoeffs(logP=logP_samp, pKa=pKa_samp, fup=fup_samp, BP=BP_samp, type=type, pred="Schmitt", dat_uni_samp)
  Kps5 <- pcoeffs(logP=logP_samp, pKa=pKa_samp, fup=fup_samp, BP=BP_samp, type=type, pred="pksim", dat_uni_samp)

  # Run PBPK model predictions to generate range of predictions for each method
  out1 <- outFun(Kps1) %>%
    mutate(Method = "PT") %>%
    mutate(ID = j)
  out2 <- outFun(Kps2) %>%
    mutate(Method = "Berez") %>%
    mutate(ID = j)
  out3 <- outFun(Kps3) %>%
    mutate(Method = "RR") %>%
    mutate(ID = j)
  out4 <- outFun(Kps4) %>%
    mutate(Method = "Schmitt") %>%
    mutate(ID = j)
  out5 <- outFun(Kps5) %>%
    mutate(Method = "PK-Sim") %>%
    mutate(ID = j)
  
  # Row-bind the predictions for each patient using each method
  out_PT <- rbind(out_PT, out1)
  out_Berez <- rbind(out_Berez, out2)
  out_RR <- rbind(out_RR, out3)
  out_Schmitt <- rbind(out_Schmitt, out4)
  out_pksim <- rbind(out_pksim, out5)
  
  # Calcuate and record the AUC and Cmax values for each method
  PK_PT[j,] <- c(pk.calc.auc(out1$Cplasma,out1$time,interval=c(sample[1],last(sample))), max(out1$Cplasma)) 
  PK_Berez[j,] <- c(pk.calc.auc(out2$Cplasma,out2$time,interval=c(sample[1],last(sample))), max(out2$Cplasma)) 
  PK_RR[j,] <- c(pk.calc.auc(out3$Cplasma,out3$time,interval=c(sample[1],last(sample))), max(out3$Cplasma)) 
  PK_Schmitt[j,] <- c(pk.calc.auc(out4$Cplasma,out4$time,interval=c(sample[1],last(sample))), max(out4$Cplasma)) 
  PK_pksim[j,] <- c(pk.calc.auc(out5$Cplasma,out5$time,interval=c(sample[1],last(sample))), max(out5$Cplasma)) 
}

# Store output of simulation that is used for the paper (figure 8 and CV calculations)
# Cplasma_all <- cbind(out_PT$time[1:1001], out_PT$Cplasma, out_Berez$Cplasma, out_RR$Cplasma, out_Schmitt$Cplasma, out_pksim$Cplasma)
# PK_all <- cbind(PK_PT, PK_Berez, PK_RR, PK_Schmitt, PK_pksim)
# write.csv(Cplasma_all, file = "../data/Cplasma_all.csv")
# write.csv(PK_all, file = "../data/PK_all.csv")

# Call output of simulation that is used for the paper
# PK_all <- read.csv("../data/PK_all.csv")
# PK_PT <- PK_all[,2:3]
# PK_Berez <- PK_all[,4:5]
# PK_RR <- PK_all[,6:7]
# PK_Schmitt <- PK_all[,8:9]
# PK_pksim <- PK_all[,10:11]

# Calculate the percent coefficient of variation for AUC and Cmax values for each method
CV_calc <- function(x){
  CV <- sig((sd(x) / mean(x))*100)
  return(CV)
}
CV_PT <- c(CV_calc(PK_PT[,1]), CV_calc(PK_PT[,2]))
CV_Berez <- c(CV_calc(PK_Berez[,1]), CV_calc(PK_Berez[,2]))
CV_RR <- c(CV_calc(PK_RR[,1]), CV_calc(PK_RR[,2]))
CV_Schmitt <- c(CV_calc(PK_Schmitt[,1]), CV_calc(PK_Schmitt[,2]))
CV_pksim <- c(CV_calc(PK_pksim[,1]), CV_calc(PK_pksim[,2]))

CV_all <- rbind(CV_PT, CV_Berez, CV_RR, CV_Schmitt, CV_pksim) %>%
  'colnames<-' (c("AUC CV", "Cmax CV")) %>%
  'rownames<-' (c("PT", "Berez", "RR", "Schmitt", "pksim"))
  
CV_table <- kable(CV_all,"html") %>%
  kable_styling(full_width=F) 
  

# Calculate the 95% prediction interval for each method
time_all <- out_PT$time[1:1001] # time points from the numerical simulation in mrgsolve

quant_calc <- function(Cplasma){
  Cplasma_all <- matrix(Cplasma, nrow = 10, ncol = 1001, byrow = TRUE)
  quant <- matrix(0, 1001, 2) 
  for (i in 1:1001){
    quant[i,] <- quantile(Cplasma_all[,i], probs=c(0.025,0.975))
  } 
  quant <- cbind(time_all, quant) %>%
    'colnames<-' (c("time_all", "lower_b", "upper_b")) %>%
    as.data.frame()
  return(quant)
}

quant_PT <- quant_calc(out_PT$Cplasma)
quant_Berez <- quant_calc(out_Berez$Cplasma)
quant_RR <- quant_calc(out_RR$Cplasma)
quant_Schmitt <- quant_calc(out_Schmitt$Cplasma)
quant_pksim <- quant_calc(out_pksim$Cplasma)


# Generate the figure
th6 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15),
             legend.position="none",
             plot.title = element_text(face="bold", size=15))

pred_PT <- ggplot() +
  geom_ribbon(data=quant_PT, aes(x=time_all, ymax=upper_b, ymin=lower_b), fill="grey60") +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("a    PT") +
  #ggtitle("PT") +
  th6

pred_Berez <- ggplot() +
  geom_ribbon(data=quant_Berez, aes(x=time_all, ymax=upper_b, ymin=lower_b), fill="grey60") +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("b    Berez") +
  #ggtitle("Berez") +
  th6

pred_RR <- ggplot() +
  geom_ribbon(data=quant_RR, aes(x=time_all, ymax=upper_b, ymin=lower_b), fill="grey60") +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("c    RR") +
  #ggtitle("RR") +
  scale_color_manual(values=c('grey60')) +
  th6


pred_Schmitt <- ggplot() +
  geom_ribbon(data=quant_Schmitt, aes(x=time_all, ymax=upper_b, ymin=lower_b), fill="grey60") +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("d    Schmitt") +
  #ggtitle("Schmitt") +
  scale_color_manual(values=c('grey60')) +
  th6


pred_pksim <- ggplot() +
  geom_ribbon(data=quant_pksim, aes(x=time_all, ymax=upper_b, ymin=lower_b), fill="grey60") +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("e    PK-Sim") +
  #ggtitle("PK-Sim") +
  scale_color_manual(values=c('grey60')) +
  th6



# Unused code that plots the full range of predictions

# out_PT_mean <- out_PT %>% 
#   group_by(time) %>%
#   dplyr::summarise(avg=mean(Cplasma))
# pred_PT <- ggplot() +
#   geom_line(data=out_PT,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
#  # geom_line(data=out_PT_mean,aes(x=time, y=avg)) +
#   geom_point(data=df, aes(x=time, y=conc), size=2.5) +
#   geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
#   xlim(-0.1,10.1) +
#   ylim(-0.1,8) +
#   xlab("Time (h)") +
#   ylab("Plasma concentration (mg/L)") +
#   ggtitle("a    PT") +
#   scale_color_manual(values=c('grey60')) +
#   th6 
# 
# out_Berez_mean <- out_Berez %>% 
#   group_by(time) %>%
#   dplyr::summarise(avg=mean(Cplasma))
# pred_Berez <- ggplot() +
#   geom_line(data=out_Berez,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
#  # geom_line(data=out_Berez_mean,aes(x=time, y=avg)) +
#   geom_point(data=df, aes(x=time, y=conc), size=2.5) +
#   geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
#   xlim(-0.1,10.1) +
#   ylim(-0.1,8) +
#   xlab("Time (h)") +
#   ylab("Plasma concentration (mg/L)") +
#   ggtitle("b    Berez") +
#   scale_color_manual(values=c('grey60')) +
#   th6 
# 
# out_RR_mean <- out_RR %>% 
#   group_by(time) %>%
#   dplyr::summarise(avg=mean(Cplasma))
# pred_RR <- ggplot() +
#   geom_line(data=out_RR,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
#  # geom_line(data=out_RR_mean,aes(x=time, y=avg)) +
#   geom_point(data=df, aes(x=time, y=conc), size=2.5) +
#   geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
#   xlim(-0.1,10.1) +
#   ylim(-0.1,8) +
#   xlab("Time (h)") +
#   ylab("Plasma concentration (mg/L)") +
#   ggtitle("c    RR") +
#   scale_color_manual(values=c('grey60')) +
#   th6
# 
# out_Schmitt_mean <- out_Schmitt %>% 
#   group_by(time) %>%
#   dplyr::summarise(avg=mean(Cplasma))
# pred_Schmitt <- ggplot() +
#   geom_line(data=out_Schmitt,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
#  # geom_line(data=out_Schmitt_mean,aes(x=time, y=avg)) +
#   geom_point(data=df, aes(x=time, y=conc), size=2.5) +
#   geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
#   xlim(-0.1,10.1) +
#   ylim(-0.1,8) +
#   xlab("Time (h)") +
#   ylab("Plasma concentration (mg/L)") +
#   ggtitle("d    Schmitt") + 
#   scale_color_manual(values=c('grey60')) +
#   th6
# 
# out_pksim_mean <- out_pksim %>% 
#   group_by(time) %>%
#   dplyr::summarise(avg=mean(Cplasma))
# pred_pksim <- ggplot() +
#   geom_line(data=out_pksim,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
#  # geom_line(data=out_pksim_mean,aes(x=time, y=avg)) +
#   geom_point(data=df, aes(x=time, y=conc), size=2.5) +
#   geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
#   xlim(-0.1,10.1) +
#   ylim(-0.1,8) +
#   xlab("Time (h)") +
#   ylab("Plasma concentration (mg/L)") +
#   ggtitle("e    PK-Sim") + 
#   scale_color_manual(values=c('grey60')) +
#   th6


# Note: certain combinations of parameters result is very large Kp estimates from Rodgers 
# and Rowland, Schmitt, and PK-Sim standard. This results in the oddly shaped lower 
# bound of the prediction interval for these methods.

fig8 <- grid.arrange(pred_PT, pred_Berez, pred_RR, pred_Schmitt, pred_pksim, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig8_95_interval.jpg", fig8, width=8, height=12)
