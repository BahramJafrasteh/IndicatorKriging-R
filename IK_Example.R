# Indicator kriging example file

# Input data should be comma delimited csv 
# where the last column is the target
# or Matlab .mat file (Xtr, Ytr, Xt, Yt)


cat("\014")
rm(list = ls())
# Required Libraries
library(caret)
library(lattice)
library(ggplot2)
library(R.matlab)
library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 

# Adding Source file
source("/media//Indicator_Krigign.R")

# Data Path
path <- "/media/Data";
setwd(path)
# warning off
options(warn=-1)
# CDF values of Thresholds (You are free to chnage them between (0, 1) )
CDFThs = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
# Output File Name
OutFileName = paste(path ,"/Out_Test", sep = "")
# Format of Input Data
DataFormat = "csv" # csv
# Input File Name
# csv
DataFileName = "Sarch_train_1.csv"
Indexes_Inputs = "Sarch_test_1.csv"
if ( isEqual(DataFormat, "Matlab") )
{
  # matlab
  DataFileName <- paste("data_divide_NN", toString(1),
                        ".mat", sep="")
  Indexes_Inputs = c(2,3,4)
}

# Indicator Kriging Algorithm
IK = IndicatorKriging ( DataFileName, Indexes_Inputs = Indexes_Inputs, CDFThs = CDFThs,
                  DataFormat = DataFormat, OutFileName = OutFileName )






