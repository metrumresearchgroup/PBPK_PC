# Impact of Partition Coefficient Prediction Methods on PBPK Output Manuscript
## Layout:
- PBPK model files are under `model`:
  - Alfentanil model `alfentanilPBPK_Adult.cpp`
  - Artemether model `artemetherPBPK_adult.cpp`
  - Caffeine model `caffeinePBPK_adult.cpp`
  - Digoxin model `digoxinPBPK_Adult.cpp`
  - Metoprolol model `metoprololPBPK_Adult.cpp`
  - Midazolam model `midazolamPBPK_Adult.cpp`
  - Nevirapine model `nevirapinePBPK_Adult.cpp`
  - Nifedipine model `nifedipinePBPK_Adult.cpp`
  - Ofloxacin model `ofloxacinPBPK_Adult.cpp`
  - Thiopental model `thiopentalPBPK_Adult.cpp`
  - Voriconazole model `voriPBPK.cpp`
- Digitized observed data are under `data`:
  - Observed alfentanil data `obs_alfentanil.csv`
  - Observed artemether data `obs_artemether.csv`  
  - Observed caffeine data `obs_caffeine.csv`
  - Observed digoxin data `obs_iv_digoxin.csv`
  - Observed metoprolol data `obs_metoprolol.csv`
  - Observed midazolam data `obs_midazolam.csv`
  - Observed nevirapine data `obs_nevirapine.csv`
  - Observed nifedipine data `obs_nifedipine.csv`
  - Observed ofloxacin data `obs_iv_ofloxacin.csv`
  - Observed S-thiopental data `obs_thiopental.csv`
  - Observed R-thiopental data `obs_thiopental_R.csv`
  - Observed voriconazole data `obs_voriconazole.csv`
- Tissue composition data are under `data`:
  - Poulin and Theil `tissue_comp_P&T.csv`
  - Rodgers and Rowland `tissue_comp_R&R.csv`
  - Schmitt `tissue_comp_Schmitt.csv`
  - PK-Sim default `PKSim_tissue_comp_pksim.csv`
  - Standardized tissue composition `unified_tissue_comp.csv` 
- Scripts for calculating tissue:plasma partition coefficients are under `script`:
  - Poulin and Theil `CalcKp_P&T.R`
  - Poulin and Theil Vss (used in verification) `CalcVss_P&T.R`
  - Berezhkovskiy `CalcKp_Berez.R`
  - Rodgers and Rowland `CalcKp_R&R.R`
  - Rodgers and Rowland Kpu (used in verification) `CalcKpu_R&R.R`
  - Schmitt `CalcKp_Schmitt.R`
  - PK-Sim default `CalcKp_pksim.R`
- Scripts for running simulations are under `script`:   
  - Alfentanil simulation `PBPK_sim_alfentanil.R`
  - Artemether simulation `PBPK_sim_artemether.R`
  - Caffeine simulation `PBPK_sim_caffeine.R`
  - Digoxin simulation `PBPK_sim_digoxin.R`
  - Metoprolol simulation `PBPK_sim_metoprolol.R`
  - Midazolam simulation `PBPK_sim_midazolam.R`
  - Nevirapine simulation `PBPK_sim_nevirapine.R`
  - Nifedipine simulation `PBPK_sim_nifedipine.R`
  - Ofloxacin simulation `PBPK_sim_ofloxacin.R`
  - Thiopental simulation `PBPK_sim_thiopental.R`
  - Voriconazole simulation `PBPK_sim_voriconazole.R`
- Script for reproducting manuscript figures is under `script`

## To reproduce manuscript figures:
1. Make the folder: `script` your working directory
2. Run `pkgSetup.R` to install required packages
3. Run each section of `script.R` to reproduce **Figures 2-7** and **Table 2** 

## Update (April 25, 2022):
`CalcKp_Schmitt.R` was updated to correct a typo in the calculation of K_n_l for diprotic bases.

## Update (January 26, 2023):
`CalcKp_R&R.R` was updated to correct a typo in the calculation of Kpu_bc, which impacts the calculation of Kp values for moderate to strong bases. 