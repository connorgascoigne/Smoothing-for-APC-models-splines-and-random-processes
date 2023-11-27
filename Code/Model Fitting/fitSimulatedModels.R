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
textSize <- 25

# load ----

## functions ----

source(paste0(codeDir, '/functions.R'))

## data ----

allData <- readRDS(file = paste0(dataDir, '/allSimulatedData.rds'))

# model fitting ----

## fit models ----

nSims <- length(allData)

# set up cluster 
nCores <- parallel::detectCores() - 1
myCluster <- parallel::makeCluster(nCores, type = "PSOCK")
doParallel::registerDoParallel(cl = myCluster)
# timed parallel loop
ptm <- proc.time() # start timer
allResults <- foreach::foreach(
  i = 1:nSims,
  # list packages that are used in loop
  .packages = c('mgcv', 'INLA', 'tidyverse')
) %dopar% {
  # cr spline basis
  crSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'cr',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))
  # bs spline basis
  bsSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'bs',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))
  # tp spline basis
  tpSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'tp',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))
  # rw2 U = 1
  rw2PC1 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 1, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)
  # rw2 U = 3
  rw2PC2 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 3, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)
  
  # rw2 U = 6
  rw2PC3 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 6, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)
  
  list(crSpline = crSpline, 
       bsSpline = bsSpline, 
       tpSpline = tpSpline,
       rw2PC1 = rw2PC1, 
       rw2PC2 = rw2PC2, 
       rw2PC3 = rw2PC3)
  
}
proc.time() - ptm # stop timer
# close cluster
parallel::stopCluster(cl = myCluster)

## collect results ----

CI <- 0.95

for(i in 1:nSims){
  
  estimationResults <- collect.simulation.results(allBasisResults = allResults[[i]], trueData = allData[[i]], CI = CI, periods = 2000:2017)
  predictionResults <- collect.simulation.results(allBasisResults = allResults[[i]], trueData = allData[[i]], CI = CI, periods = 2018:2020)
  
  scoresTemp <- rbind(estimationResults, predictionResults)
  
  if(i == 1){
    # make sure we are starting from scratch
    scores <- NULL
    scores <- scoresTemp
  } else {
    scores <- cbind(scores, scoresTemp %>% dplyr::select(-periods, -model, -metric))
  }
  
}

colnames(scores) <- c('periods', 'model', 'metric', paste0('theta:', 1:nSims))
allScoresFinal <-
  scores  %>% 
  dplyr::mutate(model = model %>% factor(., levels = c('crSpline', 'bsSpline', 'tpSpline', 'rw2PC1', 'rw2PC2', 'rw2PC3'),
                                         labels = c('CRS', 'BS', 'TPRS', 'U = 1', 'U = 3', 'U = 6')),
                type = periods %>% factor(., levels = c('2000:2017', '2018:2020'), labels = c('Estimation', 'Prediction')),
                dplyr::select(., starts_with('theta:')) %>% 
                  apply(., 1, my.summary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .)) %>% 
  dplyr::select(-starts_with('theta:'))

# results ----

## plots ----

maeBP <- 
  ggplot2::ggplot(data = allScoresFinal %>% dplyr::filter(metric == 'mae'), aes(x = model, y = mean, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Absolute Error') +
  ggplot2::facet_grid(~ type) +
  ggplot2::scale_fill_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank()); maeBP
mseBP <- 
  ggplot2::ggplot(data = allScoresFinal %>% dplyr::filter(metric == 'mse'), aes(x = model, y = mean, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Square Error') +
  ggplot2::scale_fill_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank()); mseBP

ggplot2::ggsave(filename = paste0(resultsDir, '/maeBP.png'),
                plot = maeBP,
                width = width, height = height)

ggplot2::ggsave(filename = paste0(resultsDir, '/mseBP.png'),
                plot = mseBP,
                width = width, height = height)

## interval scores ----

allScoresTable <- 
  allScoresFinal %>% 
  dplyr::filter(metric %in% c('is', 'width')) %>% 
  dplyr::select(model, metric, type, mean) %>% 
  dplyr::mutate(mean = round(mean*100, digits = 2)) %>% 
  tidyr::pivot_wider(., names_from = 'type', values_from = 'mean') %>% 
  tidyr::pivot_wider(., names_from = 'metric', values_from = c('Estimation', 'Prediction')); allScoresTable


print(xtable::xtable(allScoresTable), 
      include.rownames = FALSE,
      file = paste0(resultsDir, '/allScoreTableSimulatedData.txt'))
