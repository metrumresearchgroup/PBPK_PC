######################################### Voriconazole - weak base ###################################
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
source("getPartitionCoeff.R")

# Set functions 
filter <- dplyr::filter
mutate <- dplyr::mutate
select <- dplyr::select

# Compile model
mod <- mread("../model/voriPBPK_Adult")

# Load data (digitized observations and then calculated mean and sd)
df <- read.csv("../data/obs_voriconazole.csv")
sample <- df$time 

# Unified tissue composition with coefficient of variation
dat_uni_cv <- read.csv("../data/unified_tissue_comp_cv.csv")  
dat_uni <- read.csv("../data/unified_tissue_comp.csv")

# Number of simulated individuals (500)
num_ind <- 15

# Number of replicates 500)
num_replicates <- 15

# Total number of simulations
num_sim <- num_ind*num_replicates

# Track replication number and ID number for each simulation
rep_number <- rep(c(1:num_replicates), each=num_ind)
ind_number <- rep(c(1:num_ind), times=num_replicates)

set.seed(2891)

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

# Store PK parameters (AUC and Cmax)
PK_PT <- matrix(0,num_sim, 2)
PK_Berez <- PK_PT
PK_RR <- PK_PT
PK_Schmitt <- PK_PT
PK_pksim <- PK_PT

# Sample drug parameters from a uniform distribution with 100% increase from reported value as max and 
# 50% decrease from reported value as min.
# Modify to handle negative logP values
if (logP<0){
  logP_samp <- runif(num_replicates, min=logP+logP/2, max=logP/2)
}else{
  logP_samp <- runif(num_replicates, min=logP/2, max=logP+logP/2)
}
pKa_samp <- runif(num_replicates, min=pKa/2, max=pKa+pKa/2) # restrict to be non-negative, note: this will need to be adjusted if simulating zwitterions
fup_samp <- runif(num_replicates, min=fup/2, max=fup+fup/2)
BP_samp <- runif(num_replicates, min=BP/2, max=BP+BP/2)


# Create empty data frames to store simulation results
out_PT <- data.frame(ID=double(),RP=double(),time=double(),LUNG=double(),ADIPOSE=double(),BONE=double(),BRAIN=double(),
                     HEART=double(),KIDNEY=double(),MUSCLE=double(),SKIN=double(),LIVER=double(),SPLEEN=double(),
                     GUT=double(),ART=double(),VEN=double(),D=double(),PANCREAS=double(),REST=double(),
                     Cplasma=double(),Cvenous=double(),Method=character())
out_Berez <- out_PT
out_RR <- out_PT  
out_Schmitt <- out_PT
out_pksim <- out_PT



for (k in 1:num_sim){
  
   print(k) #track progress through the simulation
  
  # Identify the replicate and individual for each simulation
  rep <- rep_number[k]; 
  ind <- ind_number[k];
  
    
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

    # Calculate Kp values for sampled parameters
    Kps1 <- pcoeffs(logP=logP_samp[rep], pKa=pKa_samp[rep], fup=fup_samp[rep], BP=BP_samp[rep], type=type, pred="P&T", dat=dat_uni_samp)
    Kps2 <- pcoeffs(logP=logP_samp[rep], pKa=pKa_samp[rep], fup=fup_samp[rep], BP=BP_samp[rep], type=type, pred="Berez", dat=dat_uni_samp)
    Kps3 <- pcoeffs(logP=logP_samp[rep], pKa=pKa_samp[rep], fup=fup_samp[rep], BP=BP_samp[rep], type=type, pred="R&R", dat=dat_uni_samp)
    Kps4 <- pcoeffs(logP=logP_samp[rep], pKa=pKa_samp[rep], fup=fup_samp[rep], BP=BP_samp[rep], pred="Schmitt", dat=dat_uni_samp)
    Kps5 <- pcoeffs(logP=logP_samp[rep], pKa=pKa_samp[rep], fup=fup_samp[rep], BP=BP_samp[rep], pred="pksim", dat=dat_uni_samp)
  
    # Run PBPK model to generate predictions for each method
    out1 <- outFun(Kps1) %>%
      mutate(Method = "PT") %>%
      mutate(ID = ind) %>%
      mutate(RP = rep)
    out2 <- outFun(Kps2) %>%
      mutate(Method = "Berez") %>%
      mutate(ID = ind) %>%
      mutate(RP = rep)
    out3 <- outFun(Kps3) %>%
      mutate(Method = "RR") %>%
      mutate(ID = ind) %>%
      mutate(RP = rep)
    out4 <- outFun(Kps4) %>%
      mutate(Method = "Schmitt") %>%
      mutate(ID = ind) %>%
      mutate(RP = rep)
    out5 <- outFun(Kps5) %>%
      mutate(Method = "PK-Sim") %>%
      mutate(ID = ind) %>%
      mutate(RP = rep)
    
    # Row-bind the predictions for each patient using each method
    out_PT <- rbind(out_PT, out1)
    out_Berez <- rbind(out_Berez, out2)
    out_RR <- rbind(out_RR, out3)
    out_Schmitt <- rbind(out_Schmitt, out4)
    out_pksim <- rbind(out_pksim, out5)
    
}

# Store output of simulation
#vori_all <- rbind(out_PT,out_Berez,out_RR,out_Schmitt,out_pksim)
#vori_reduc <- vori_all %>%
#  select(ID, time, Cplasma, Method, RP)
#write.csv(vori_reduc, file = "../data/vori_reduc.csv")

# Read in output of simulation (50x50 scheme)
#vori_reduc <- read.csv("../data/vori_reduc.csv")
# out_PT <- vori_reduc %>%
#   filter(Method == "PT")
# out_Berez <- vori_reduc %>%
#   filter(Method == "Berez")
# out_RR <- vori_reduc %>%
#   filter(Method == "RR")
# out_Schmitt <- vori_reduc %>%
#   filter(Method == "Schmitt")
# out_pksim <- vori_reduc %>%
#   filter(Method == "PK-Sim")


# Save the outputs from each method separately (file size issues)
# out_PT <- out_PT %>%
#   select(ID, time, Cplasma, Method, RP)
# out_Berez <- out_Berez %>%
#   select(ID, time, Cplasma, Method, RP)
# out_RR <- out_RR %>%
#   select(ID, time, Cplasma, Method, RP)
# out_Schmitt <- out_Schmitt %>%
#   select(ID, time, Cplasma, Method, RP)
# out_pksim <- out_pksim %>%
#   select(ID, time, Cplasma, Method, RP)
# write.csv(out_PT, file = "../data/vori_PT.csv")
# write.csv(out_Berez, file = "../data/vori_Berez.csv")
# write.csv(out_RR, file = "../data/vori_RR.csv")
# write.csv(out_Schmitt, file = "../data/vori_Schmitt.csv")
# write.csv(out_pksim, file = "../data/vori_pksim.csv")

# Read in output of simulation
out_PT <- read.csv("../data/vori_PT.csv")
out_Berez <- read.csv("../data/vori_Berez.csv")
out_RR <- read.csv("../data/vori_RR.csv")
out_Schmitt <- read.csv("../data/vori_Schmitt.csv")
out_pksim <- read.csv("../data/vori_pksim.csv")


# Define function to calculate the 10th/90th percentiles and 95% prediction intervals
med <- function(x) as.numeric(quantile(x, 0.5))
lo1 <- function(x) as.numeric(quantile(x, 0.025))
hi1 <- function(x) as.numeric(quantile(x, 0.975))
lo2 <- function(x) as.numeric(quantile(x, 0.1))
hi2 <- function(x) as.numeric(quantile(x, 0.9))


calc_stat <- function(output){
  out_stat <- output %>%
    # summarize first across individuals, get average conc for each subject at each timepoint
    group_by(RP, time) %>%
    mutate(med = med(Cplasma),
           lo = lo2(Cplasma),
           hi = hi2(Cplasma)) %>%
    ungroup() %>%
    # now get summary across replicates
    group_by(time) %>%
    mutate(medMed = med(med),
           loMed = lo1(med),
           hiMed = hi1(med),
           medLo = med(lo),
           loLo = lo1(lo),
           hiLo = hi1(lo),
           medHi = med(hi),
           loHi = lo1(hi),
           hiHi = hi1(hi)) %>%
    ungroup()
  return(out_stat)
}

PT_stat <- calc_stat(out_PT)
Berez_stat <- calc_stat(out_Berez)
RR_stat <- calc_stat(out_RR)
Schmitt_stat <- calc_stat(out_Schmitt)
pksim_stat <- calc_stat(out_pksim)


# Calcuate and record the AUC and Cmax values for each method
PK_PT <- out_PT %>%
  group_by(ID, RP) %>%
  mutate(AUC = pk.calc.auc(Cplasma, time,interval=c(sample[1],last(sample))),
         Cmax = max(Cplasma))
PK_PT <- cbind(PK_PT[seq(1, length(out_PT$time), 1001), 7], PK_PT[seq(1, length(out_PT$time), 1001), 8])

PK_Berez <- out_Berez %>%
  group_by(ID, RP) %>%
  mutate(AUC = pk.calc.auc(Cplasma, time,interval=c(sample[1],last(sample))),
         Cmax = max(Cplasma))
PK_Berez <- cbind(PK_Berez[seq(1, length(out_PT$time), 1001), 7], PK_Berez[seq(1, length(out_PT$time), 1001), 8])

PK_RR <- out_RR %>%
  group_by(ID, RP) %>%
  mutate(AUC = pk.calc.auc(Cplasma, time,interval=c(sample[1],last(sample))),
         Cmax = max(Cplasma))
PK_RR <- cbind(PK_RR[seq(1, length(out_PT$time), 1001), 7], PK_RR[seq(1, length(out_PT$time), 1001), 8])

PK_Schmitt <- out_Schmitt %>%
  group_by(ID, RP) %>%
  mutate(AUC = pk.calc.auc(Cplasma, time,interval=c(sample[1],last(sample))),
         Cmax = max(Cplasma))
PK_Schmitt <- cbind(PK_Schmitt[seq(1, length(out_PT$time), 1001), 7], PK_Schmitt[seq(1, length(out_PT$time), 1001), 8])

PK_pksim <- out_pksim %>%
  group_by(ID, RP) %>%
  mutate(AUC = pk.calc.auc(Cplasma, time,interval=c(sample[1],last(sample))),
         Cmax = max(Cplasma))
PK_pksim <- cbind(PK_pksim[seq(1, length(out_PT$time), 1001), 7], PK_pksim[seq(1, length(out_PT$time), 1001), 8])


# Calculate the percent coefficient of variation for AUC and Cmax values for each method
CV_calc <- function(x){
  CV <- sig((sd(x) / mean(x))*100)
  return(CV)
}
CV_PT <- c(CV_calc(PK_PT$AUC), CV_calc(PK_PT$Cmax))
CV_Berez <- c(CV_calc(PK_Berez$AUC), CV_calc(PK_Berez$Cmax))
CV_RR <- c(CV_calc(PK_RR$AUC), CV_calc(PK_RR$Cmax))
CV_Schmitt <- c(CV_calc(PK_Schmitt$AUC), CV_calc(PK_Schmitt$Cmax))
CV_pksim <- c(CV_calc(PK_pksim$AUC), CV_calc(PK_pksim$Cmax))


CV_all <- rbind(CV_PT, CV_Berez, CV_RR, CV_Schmitt, CV_pksim) %>%
  'colnames<-' (c("AUC CV", "Cmax CV")) %>%
  'rownames<-' (c("PT", "Berez", "RR", "Schmitt", "pksim"))

CV_table <- kable(CV_all,"html") %>%
  kable_styling(full_width=F) 



##
th6 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 15), plot.title = element_text(face="bold", size=15), 
             legend.title = element_blank(), legend.justification=c(1,1), 
             legend.position=c(1,1))
cols <- c("Median"="black","10%"="red","90%"="blue");
pred_PT <-  
  ggplot() +
  geom_ribbon(data=PT_stat[1:1001,], aes(x=time, ymax=hiLo, ymin=loLo, fill="10%"), alpha=0.4) +
  geom_ribbon(data=PT_stat[1:1001,], aes(x=time, ymax=hiHi, ymin=loHi, fill="90%"), alpha=0.4) +
  geom_ribbon(data=PT_stat[1:1001,], aes(x=time, ymax=hiMed, ymin=loMed, fill="Median"), alpha=0.4) +
  geom_line(data=PT_stat[1:1001,], aes(x=time, y=medMed, colour="Median")) +
  geom_line(data=PT_stat[1:1001,], aes(x=time, y=medLo, colour="10%")) +
  geom_line(data=PT_stat[1:1001,], aes(x=time, y=medHi, colour="90%")) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  scale_colour_manual("",values=cols) +
  scale_fill_manual("",values=c("Median"="grey60","10%"="red","90%"="blue")) +
  #guides(fill = guide_legend(override.aes = list(colour="grey60"))) +
  guides(fill = guide_legend(override.aes= list(alpha = 0.2))) +
  xlim(-0.1,10.1) +
  ylim(-0.1,9) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("a    PT") +
  th6 

pred_Berez <- 
  ggplot() +
  geom_ribbon(data=Berez_stat[1:1001,], aes(x=time, ymax=hiLo, ymin=loLo, fill="10%"), alpha=0.4) +
  geom_ribbon(data=Berez_stat[1:1001,], aes(x=time, ymax=hiHi, ymin=loHi, fill="90%"), alpha=0.4) +
  geom_ribbon(data=Berez_stat[1:1001,], aes(x=time, ymax=hiMed, ymin=loMed, fill="Median"), alpha=0.4) +
  geom_line(data=Berez_stat[1:1001,], aes(x=time, y=medMed, colour="Median")) +
  geom_line(data=Berez_stat[1:1001,], aes(x=time, y=medLo, colour="10%")) +
  geom_line(data=Berez_stat[1:1001,], aes(x=time, y=medHi, colour="90%")) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  scale_colour_manual("",values=cols) +
  scale_fill_manual("",values=c("Median"="grey60","10%"="red","90%"="blue")) +
  #guides(fill = guide_legend(override.aes = list(colour="grey60"))) +
  guides(fill = guide_legend(override.aes= list(alpha = 0.2))) +
  xlim(-0.1,10.1) +
  ylim(-0.1,9) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("b    Berez") +
  th6 

pred_RR <- 
  ggplot() +
  geom_ribbon(data=RR_stat[1:1001,], aes(x=time, ymax=hiLo, ymin=loLo, fill="10%"), alpha=0.4) +
  geom_ribbon(data=RR_stat[1:1001,], aes(x=time, ymax=hiHi, ymin=loHi, fill="90%"), alpha=0.4) +
  geom_ribbon(data=RR_stat[1:1001,], aes(x=time, ymax=hiMed, ymin=loMed, fill="Median"), alpha=0.4) +
  geom_line(data=RR_stat[1:1001,], aes(x=time, y=medMed, colour="Median")) +
  geom_line(data=RR_stat[1:1001,], aes(x=time, y=medLo, colour="10%")) +
  geom_line(data=RR_stat[1:1001,], aes(x=time, y=medHi, colour="90%")) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  scale_colour_manual("",values=cols) +
  scale_fill_manual("",values=c("Median"="grey60","10%"="red","90%"="blue")) +
  #guides(fill = guide_legend(override.aes = list(colour="grey60"))) +
  guides(fill = guide_legend(override.aes= list(alpha = 0.2))) +
  xlim(-0.1,10.1) +
  ylim(-0.1,9) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("c    RR") +
  th6 

pred_Schmitt <- 
  ggplot() +
  geom_ribbon(data=Schmitt_stat[1:1001,], aes(x=time, ymax=hiLo, ymin=loLo, fill="10%"), alpha=0.4) +
  geom_ribbon(data=Schmitt_stat[1:1001,], aes(x=time, ymax=hiHi, ymin=loHi, fill="90%"), alpha=0.4) +
  geom_ribbon(data=Schmitt_stat[1:1001,], aes(x=time, ymax=hiMed, ymin=loMed, fill="Median"), alpha=0.4) +
  geom_line(data=Schmitt_stat[1:1001,], aes(x=time, y=medMed, colour="Median")) +
  geom_line(data=Schmitt_stat[1:1001,], aes(x=time, y=medLo, colour="10%")) +
  geom_line(data=Schmitt_stat[1:1001,], aes(x=time, y=medHi, colour="90%")) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  scale_colour_manual("",values=cols) +
  scale_fill_manual("",values=c("Median"="grey60","10%"="red","90%"="blue")) +
  #guides(fill = guide_legend(override.aes = list(colour="grey60"))) +
  guides(fill = guide_legend(override.aes= list(alpha = 0.2))) +
  xlim(-0.1,10.1) +
  ylim(-0.1,9) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("d    Schmitt") +
  th6 

pred_pksim <- 
  ggplot() +
  geom_ribbon(data=pksim_stat[1:1001,], aes(x=time, ymax=hiLo, ymin=loLo, fill="10%"), alpha=0.4) +
  geom_ribbon(data=pksim_stat[1:1001,], aes(x=time, ymax=hiHi, ymin=loHi, fill="90%"), alpha=0.4) +
  geom_ribbon(data=pksim_stat[1:1001,], aes(x=time, ymax=hiMed, ymin=loMed, fill="Median"), alpha=0.4) +
  geom_line(data=pksim_stat[1:1001,], aes(x=time, y=medMed, colour="Median")) +
  geom_line(data=pksim_stat[1:1001,], aes(x=time, y=medLo, colour="10%")) +
  geom_line(data=pksim_stat[1:1001,], aes(x=time, y=medHi, colour="90%")) +
  geom_point(data=df, aes(x=time, y=conc), size=2.5) +
  geom_errorbar(data = df, aes(x = time, ymin=conc-sd, ymax=conc+sd), width=0.1) +
  scale_colour_manual("",values=cols) +
  scale_fill_manual("",values=c("Median"="grey60","10%"="red","90%"="blue")) +
  #guides(fill = guide_legend(override.aes = list(colour="grey60"))) +
  guides(fill = guide_legend(override.aes= list(alpha = 0.2))) +
  xlim(-0.1,10.1) +
  ylim(-0.1,9) +
  xlab("Time (h)") +
  ylab("Plasma concentration (mg/L)") +
  ggtitle("e   PK-Sim") +
  th6 

fig8_new <- grid.arrange(pred_PT, pred_Berez, pred_RR, pred_Schmitt, pred_pksim, ncol=2, nrow=3)
ggsave(file="../deliv/figure/fig8_new.png", fig8_new, width=8, height=12)



#########################################################################################################

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
# CV_calc <- function(x){
#   CV <- sig((sd(x) / mean(x))*100)
#   return(CV)
# }
# CV_PT <- c(CV_calc(PK_PT[,1]), CV_calc(PK_PT[,2]))
# CV_Berez <- c(CV_calc(PK_Berez[,1]), CV_calc(PK_Berez[,2]))
# CV_RR <- c(CV_calc(PK_RR[,1]), CV_calc(PK_RR[,2]))
# CV_Schmitt <- c(CV_calc(PK_Schmitt[,1]), CV_calc(PK_Schmitt[,2]))
# CV_pksim <- c(CV_calc(PK_pksim[,1]), CV_calc(PK_pksim[,2]))
# 
# CV_all <- rbind(CV_PT, CV_Berez, CV_RR, CV_Schmitt, CV_pksim) %>%
#   'colnames<-' (c("AUC CV", "Cmax CV")) %>%
#   'rownames<-' (c("PT", "Berez", "RR", "Schmitt", "pksim"))
# 
# CV_table <- kable(CV_all,"html") %>%
#   kable_styling(full_width=F) 


# Note: certain combinations of parameters result is very large Kp estimates from Rodgers 
# and Rowland, Schmitt, and PK-Sim standard. This results in the oddly shaped lower 
# bound of the prediction interval for these methods.
