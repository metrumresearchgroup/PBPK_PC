#This function generates a list of the desired partition coefficients
.libPaths("lib")
library(dplyr)

pcoeffs <- function(logP, pKa, fup, BP, type=3, pred="P&T", dat){
  if(pred=="P&T"){
    source("CalcKp_P&T.R")
    pcoeff <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
  }else if(pred=="Berez"){  #Berezhkovskiy 
    source("CalcKp_Berez.R")
    pcoeff <- calcKp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
  }else if(pred == "pksim"){  #standard PK-Sim, Willmann et al. 2008
    source("CalcKp_pksim.R")  
    pcoeff <- calcKp_pksim(logP=logP, fup=fup, dat=dat)
  }else if(pred == "Schmitt"){  #Schmitt, Walter 2008
    source("CalcKp_Schmitt.R")
    pcoeff <- calcKp_Schmitt(logP=logP, pKa=pKa, fup=fup, type=type, dat=dat)
  }else{
    source("CalcKp_R&R.R")  #Rodgers and Rowland 2006
    pcoeff <- calcKp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
  }

  return(pcoeff)
}
