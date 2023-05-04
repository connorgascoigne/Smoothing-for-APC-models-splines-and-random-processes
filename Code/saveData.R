# packages ----

library(tidyverse)

# directories ---- 

# extract file location of this script
codePath <- rstudioapi::getActiveDocumentContext()$path
codePathSplitted <- strsplit(codePath, "/")[[1]]

# retrieve directories
homeDir <- paste(codePathSplitted[1: (length(codePathSplitted)-2)], collapse = "/")
codeDir <- paste0(homeDir, '/Code')

# need to make sure data is made
if(!dir.exists(paths = paste0(homeDir, '/Data'))) {
  dir.create(path = paste0(homeDir, '/Data'))
}

dataDir <- paste0(homeDir, '/Data')

# load ----

## functions ----

source(paste0(codeDir, '/functions.R'))

# create data sets ----

## arguments ---

nSims <- 100; age <- 10:84; period <- 2000:2020; N <- 750000; M <- 5

## data ---

allData <- list()
for(i in 1:nSims){
  allData[[i]] <- data.sim(age = age, period = period, N = N, M = M,
                           FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
  print(i)
}

## save ----

saveRDS(allData, file = paste0(dataDir, '/allSimulatedData.rds'))
