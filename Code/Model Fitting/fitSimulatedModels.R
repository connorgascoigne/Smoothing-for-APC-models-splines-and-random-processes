# packages ----

library(tidyverse)
library(mgcv)
library(INLA)

# directories ---- 

# extract file location of this script
codePath <- rstudioapi::getActiveDocumentContext()$path
codePathSplitted <- strsplit(codePath, "/")[[1]]

# retrieve directories
homeDir <- paste(codePathSplitted[1: (length(codePathSplitted)-3)], collapse = "/")
codeDir <- paste0(homeDir, '/Code')
dataDir <- paste0(homeDir, '/Data')

# need to make sure results is made

if(!dir.exists(paths = paste0(codeDir, '/Results'))) {
  dir.create(path = paste0(codeDir, '/Results'))
}

resultsDir <- paste0(codeDir, '/Results')

# saving format ----

height <- width <- 10
textSize <- 17.5

# load ----

## functions ----

source(paste0(codeDir, '/functions.R'))

## data ----

allData <- readRDS(file = paste0(dataDir, '/allSimulatedData.rds'))

# model fitting ----

## model fit ----

nObs <- nrow(allData[[1]])
nSims <- length(allData)

crSpline <-
  bsSpline <-
  psSpline <-
  rw1PC1 <-
  rw1PC2 <-
  rw2PC1 <-
  rw2PC2 <- 
  matrix(NA, nrow = nObs, ncol = nSims)

for(i in 1:nSims){
  
  # full data
  dataFull <- allData[[i]]
  
  # data for spline estimate function
  ## estimates
  dataEst_spline <- 
    dataFull %>% 
    dplyr::filter(period %in% 2000:2017)
  ## prediction
  dataPred_spline <-
    dataFull %>% 
    dplyr::select(age, period, cohort)
  
  # data for random walk estimate function
  dataPred_randomWalk <-
    dataFull %>% 
    dplyr::mutate(y = dplyr::if_else(period > 2017, NA, y))
  
  # estimates
  ## spline
  ### cr basis
  crSpline[,i] <- spline.estimates(dataEst = dataEst_spline, dataPred = dataPred_spline,
                                   mod = 'apc', slopeDrop = 'c', bs = 'cr',
                                   knots = list(age = 5, period = 5, cohort = 8),
                                   fixed = list(age = F, period = F, cohort = F))
  ### bs
  bsSpline[,i] <- spline.estimates(dataEst = dataEst_spline, dataPred = dataPred_spline,
                                   mod = 'apc', slopeDrop = 'c', bs = 'bs',
                                   knots = list(age = 5, period = 5, cohort = 8),
                                   fixed = list(age = F, period = F, cohort = F))
  ### ps
  psSpline[,i] <- spline.estimates(dataEst = dataEst_spline, dataPred = dataPred_spline,
                                   mod = 'apc', slopeDrop = 'c', bs = 'ps',
                                   knots = list(age = 5, period = 5, cohort = 8),
                                   fixed = list(age = F, period = F, cohort = F))
  ## random walk
  ### rw1
  #### pc1
  rw1PC1[,i] <- randomWalk.estimates(dataPred = dataPred_randomWalk,
                                     mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                                     pc.u = 1, pc.alpha = 0.01,
                                     control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                                     inla.mode = c('classic', 'twostage', 'experimental')[3],
                                     control.compute = list(config = TRUE), verbose = FALSE)
  #### pc2
  rw1PC2[,i] <- randomWalk.estimates(dataPred = dataPred_randomWalk,
                                     mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                                     pc.u = 3, pc.alpha = 0.01,
                                     control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                                     inla.mode = c('classic', 'twostage', 'experimental')[3],
                                     control.compute = list(config = TRUE), verbose = FALSE)
  ### rw2
  #### pc1
  rw2PC1[,i] <- randomWalk.estimates(dataPred = dataPred_randomWalk,
                                     mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                                     pc.u = 1, pc.alpha = 0.01,
                                     control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                                     inla.mode = c('classic', 'twostage', 'experimental')[3],
                                     control.compute = list(config = TRUE), verbose = FALSE)
  #### pc2
  rw2PC2[,i] <- randomWalk.estimates(dataPred = dataPred_randomWalk,
                                     mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                                     pc.u = 3, pc.alpha = 0.01,
                                     control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                                     inla.mode = c('classic', 'twostage', 'experimental')[3],
                                     control.compute = list(config = TRUE), verbose = FALSE)
  
  print(i)
  
}

# collect results arguments
CI = 0.95
allNames <- c('crSpline',
              'bsSpline',
              'psSpline',
              'rw1PC1',
              'rw1PC2',
              'rw2PC1',
              'rw2PC2')

final <- 
  rbind(collect.results(simulationResults = eval(parse(text = allNames[1])), allData = allData, CI = CI, name = allNames[1]),
        collect.results(simulationResults = eval(parse(text = allNames[2])), allData = allData, CI = CI, name = allNames[2]),
        collect.results(simulationResults = eval(parse(text = allNames[3])), allData = allData, CI = CI, name = allNames[3]),
        collect.results(simulationResults = eval(parse(text = allNames[4])), allData = allData, CI = CI, name = allNames[4]),
        collect.results(simulationResults = eval(parse(text = allNames[5])), allData = allData, CI = CI, name = allNames[5]),
        collect.results(simulationResults = eval(parse(text = allNames[6])), allData = allData, CI = CI, name = allNames[6]),
        collect.results(simulationResults = eval(parse(text = allNames[7])), allData = allData, CI = CI, name = allNames[7])) %>% 
  dplyr::mutate(model = model %>% factor(., levels = allNames))

ggplot2::ggplot(data = final %>% dplyr::filter(period > 2017), aes(x = age)) +
  ggplot2::geom_line(aes(y = median, group = interaction(model, period),  colour = model), linetype = 'dashed') +
  ggplot2::geom_line(aes(y = lower, group = interaction(model, period),  colour = model), linetype = 'dashed') +
  ggplot2::geom_line(aes(y = upper, group = interaction(model, period),  colour = model), linetype = 'dashed') +
  # ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper, group = interaction(model, period), fill = model),
  #                      alpha = 0.25, colour = NA) +
  ggplot2::geom_line(data = final %>% dplyr::filter(period > 2017) %>% dplyr::select(age, period, cohort, true) %>% dplyr::distinct(), 
                     aes(y = true), linetype = 'dotted') +
  ggplot2::facet_wrap(~ as.factor(period)) +
  my.theme(legend.title = element_blank())

biasBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = bias, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dashed') +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Bias') +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank())

mseBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = mse, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dashed') +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Bias') +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank())


ggplot2::ggsave(filename = paste0(resultsDir, '/biasBP.png'),
                plot = biasBP,
                width = width, height = height)

ggplot2::ggsave(filename = paste0(resultsDir, '/mseBP.png'),
                plot = mseBP,
                width = width, height = height)
