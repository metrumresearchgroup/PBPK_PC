######################################### Voriconazole - weak base ###################################
# Reference: 
# Compile model
mod <- mread("../model/voriPBPK_Adult")

# Load data (digitized observations and then calculated mean and sd)
df <- read.csv("../data/obs_voriconazole.csv")

# Unified tissue composition with coefficient of variation
dat_uni_cv <- read.csv("../data/unified_tissue_comp_error.csv")  
dat_uni <- read.csv("../data/unified_tissue_comp.csv")

# Number of simulated patients (use 1000 to make figure, but fewer during de-bugging)
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
    mutate(Method = "Poulin and Theil") %>%
    mutate(ID = j)
  out2 <- outFun(Kps2) %>%
    mutate(Method = "Berezhkovskiy") %>%
    mutate(ID = j)
  out3 <- outFun(Kps3) %>%
    mutate(Method = "Rodgers and Rowland") %>%
    mutate(ID = j)
  out4 <- outFun(Kps4) %>%
    mutate(Method = "Schmitt") %>%
    mutate(ID = j)
  out5 <- outFun(Kps5) %>%
    mutate(Method = "PK-Sim Standard") %>%
    mutate(ID = j)
  
  # Row-bind the predictions for each patient using each method
  out_PT <- rbind(out_PT, out1)
  out_Berez <- rbind(out_Berez, out2)
  out_RR <- rbind(out_RR, out3)
  out_Schmitt <- rbind(out_Schmitt, out4)
  out_pksim <- rbind(out_pksim, out5)
}

th6 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15),
               legend.position="none",
               plot.title = element_text(face="bold", size=15))

out_PT_mean <- out_PT %>% 
  group_by(time) %>%
  dplyr::summarise(avg=mean(Cplasma))
pred_PT <- ggplot() +
  geom_line(data=out_PT,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
 # geom_line(data=out_PT_mean,aes(x=time, y=avg)) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("a    Poulin and Theil") +
  scale_color_manual(values=c('grey60')) +
  th6 

out_Berez_mean <- out_Berez %>% 
  group_by(time) %>%
  dplyr::summarise(avg=mean(Cplasma))
pred_Berez <- ggplot() +
  geom_line(data=out_Berez,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
 # geom_line(data=out_Berez_mean,aes(x=time, y=avg)) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("b    Berezhkovskiy") +
  scale_color_manual(values=c('grey60')) +
  th6 

out_RR_mean <- out_RR %>% 
  group_by(time) %>%
  dplyr::summarise(avg=mean(Cplasma))
pred_RR <- ggplot() +
  geom_line(data=out_RR,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
 # geom_line(data=out_RR_mean,aes(x=time, y=avg)) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("c    Rodgers and Rowland") +
  scale_color_manual(values=c('grey60')) +
  th6

out_Schmitt_mean <- out_Schmitt %>% 
  group_by(time) %>%
  dplyr::summarise(avg=mean(Cplasma))
pred_Schmitt <- ggplot() +
  geom_line(data=out_Schmitt,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
 # geom_line(data=out_Schmitt_mean,aes(x=time, y=avg)) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("d    Schmitt") + 
  scale_color_manual(values=c('grey60')) +
  th6

out_pksim_mean <- out_pksim %>% 
  group_by(time) %>%
  dplyr::summarise(avg=mean(Cplasma))
pred_pksim <- ggplot() +
  geom_line(data=out_pksim,aes(x=time, y=Cplasma, linetype=Method, color=Method)) +
 # geom_line(data=out_pksim_mean,aes(x=time, y=avg)) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  xlim(-0.1,10.1) +
  ylim(-0.1,8) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("e    PK-Sim Standard") + 
  scale_color_manual(values=c('grey60')) +
  th6


# Note: certain combinations of parameters result is very large Kp estimates from Rodgers 
# and Rowland, Schmitt, and PK-Sim standard. This results in the oddly shaped lower 
# bound of the prediction interval for these methods.

fig6 <- grid.arrange(pred_PT, pred_Berez, pred_RR, pred_Schmitt, pred_pksim, ncol=2, nrow=3)
#ggsave(file="../deliv/figure/fig6_test_no_mean.jpg", fig6, width=8, height=12)
