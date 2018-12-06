## This script installs all required packages and creates a "lib" directory
## Before you proceed make the main repo directory your working directory with the function setwd()
author <- c("ahmede")

pkgs <- c("dplyr", "ggplot2", "tidyr", "mrgsolve", "magrittr", "FME", "sensitivity", "nloptr", "gridExtra", "PKPDmisc",
          "kableExtra", "PKNCA","msm")


pkgRoot <- "pkg"
pkgDir <- file.path(pkgRoot, "src", "contrib")
pkgDir <- normalizePath(pkgDir, mustWork=FALSE)
libDir <- "lib"

if(!dir.exists(pkgDir)) dir.create(pkgDir, recursive = TRUE)
if(!dir.exists(libDir)) dir.create(libDir)

.libPaths(libDir)

user <- Sys.info()["user"]

fromCRAN <- user %in% author

local_repos <- paste0("file://",file.path(getwd(),pkgRoot))
metrum_repos <- "https://cran.rstudio.com/"
cran_repos <- "https://cran.rstudio.com/"
repos <- c(mrg = metrum_repos, cran = cran_repos, local = local_repos)

if(file.exists(file.path(pkgDir,"PACKAGES"))){
  available <- available.packages(repos = repos["local"])[,"Package"]
} else{
  available <- NULL
  file.create(file.path(pkgDir,"PACKAGES"))
  tools::write_PACKAGES(pkgDir)
}

if(fromCRAN){
  
  newpkgs <- setdiff(pkgs, available)
  
  tools::write_PACKAGES(pkgDir)
  
  if(length(newpkgs) > 0){
    ## These packages are installed either from mrg or cran
    install.packages(newpkgs,
                     lib=libDir,
                     repos = repos[c("mrg", "cran")],
                     destdir=pkgDir,
                     type="source",
                     dependencies =  c("Depends", "Imports", "LinkingTo"),
                     INSTALL_opts="--no-multiarch")
    
    tools::write_PACKAGES(pkgDir)
  }
  
  ## If multiple authors qcing each other, a package could be available
  ## but uninstalled.  Install from local.
  uninstalled <- setdiff(pkgs, installed.packages(libDir))
  
  if(length(uninstalled)>0){
    install.packages(uninstalled,
                     lib = libDir,
                     repos = repos["local"],
                     type = "source",
                     INSTALL_opts="--no-multiarch")
  }    
}


if(!fromCRAN){
  installed <- row.names(installed.packages(libDir))
  newpkgs <- setdiff(pkgs, installed)
  if(length(newpkgs)>0){
    install.packages(newpkgs,
                     lib = libDir,
                     repos = repos["local"],
                     type = "source",
                     INSTALL_opts="--no-multiarch")
    
  }
}
