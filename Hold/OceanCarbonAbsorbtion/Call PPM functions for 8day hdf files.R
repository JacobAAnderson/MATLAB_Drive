rm(list = ls())
library(RNetCDF)
library(R.matlab)
library(ncdf4)
library(reshape2)
library(raster)
library(tidyr)
library(raster)
library(ncdf4)
library(RNetCDF)
library(ggplot2)
library(dplyr)
library(plyr)

source('~/OSU/Models/NPP_models/PPM_function.R')

#Read in Data
setwd("~/OSU/Models/")

files.bbp  <- list.files(path = "Sat_data/BBP_8day_9km/", pattern = "*.hdf", full.names=T)
files.k490 <- list.files(path = "Sat_data/KD_8day_9km/", pattern = "*.hdf", full.names=T)
files.chl  <- list.files(path = "Sat_data/CHL_8day_9km/", pattern = "*.hdf", full.names=T)
files.irr  <- list.files(path = "Sat_data/PAR_8day_9km/", pattern = "*.hdf", full.names=T)
files.mld  <- list.files(path = "Sat_data/MLD_8day_9km/", pattern = "*.hdf", full.names=T)

files.order <- c(1:39)
mn          <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12","13",
                 "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25","26",
                 "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38","39")
doy         <- seq(15, 350, length.out=12 )


for (i in 1:39){
 

  bbp <- read_hdf_files(files.bbp,files.order) 
  MLD <- read_hdf_files(files.mld,files.order) 
  chl <- read_hdf_files(files.chl,files.order) 
  k490 <- read_hdf_files(files.k490,files.order) 
  irr <- read_hdf_files(files.irr,files.order) 
  
  
  #1) Regular Model, Data automatically saved
  PPM(bbp, chl, irr, k490, MLD,  doy[i], mn[i]) 
  
  #3) Return PAR at depth z
 #z       <- MLD_stretch_larger(MLD)
#PAR.MLD <- CpBM.PARz(bbp, chl, irr, K490, MLD,  doy[i], mn[i], z)
#  save(PAR.MLD, file = paste(getwd(), "/RData/PARz/PAR ", mn[i],".RData", sep="")) 
  
#  z        <- MLD_stretch(ZNO3)
#  PAR.ZNO3 <- CpBM.PARz(bbp, chl, irr, K490, MLD, ZNO3,  doy[i], mn[i], z)
#  save(PAR.ZNO3, file = paste(getwd(), "/RData/PARz/ZNO3 ", mn[i],".RData", sep=""))   
}

