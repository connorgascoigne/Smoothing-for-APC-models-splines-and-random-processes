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

nSims <- length(allData)

# set up cluster 
nCores <- parallel::detectCores() - 1
myCluster <- parallel::makeCluster(nCores, type = "PSOCK")
doParallel::registerDoParallel(cl = myCluster)
# timed parallel loop
ptm <- proc.time() # start timer
allResults <- foreach::foreach(
  i = 1:nSims,
  # list packages that are used in process
  .packages = c('mgcv', 'INLA', 'tidyverse')
) %dopar% {
  # process in loop
  # estimates
  ## spline
  ### cr basis
  crSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'cr',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  ### bs
  bsSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'bs',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  ### ps
  psSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'ps',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  ## random walk
  # ### rw1
  # #### pc1
  # rw1PC1 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
  #                              mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
  #                              pc.u = 1, pc.alpha = 0.01,
  #                              control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
  #                              inla.mode = c('classic', 'twostage', 'experimental')[3],
  #                              control.compute = list(config = TRUE), verbose = FALSE)$yHat
  # #### pc2
  # rw1PC2 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
  #                              mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
  #                              pc.u = 3, pc.alpha = 0.01,
  #                              control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
  #                              inla.mode = c('classic', 'twostage', 'experimental')[3],
  #                              control.compute = list(config = TRUE), verbose = FALSE)$yHat
  ### rw2
  #### pc1
  rw2PC1 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 1, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)$yHat
  #### pc2
  rw2PC2 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 3, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)$yHat
  
  list(crSpline = crSpline, 
       bsSpline = bsSpline, 
       psSpline = psSpline,
       # rw1PC1 = rw1PC1, 
       # rw1PC2 = rw1PC2,
       rw2PC1 = rw2PC1, 
       rw2PC2 = rw2PC2)
  
}
proc.time() - ptm # stop timer
# close cluster
parallel::stopCluster(cl = myCluster)

CI <- 0.95

final <- 
  rbind(collect.results(simulationResults = sapply(allResults, `[[`, 'crSpline'), allData = allData, CI = CI, name = 'crSpline'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'bsSpline'), allData = allData, CI = CI, name = 'bsSpline'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'psSpline'), allData = allData, CI = CI, name = 'psSpline'),
        # collect.results(simulationResults = sapply(allResults, `[[`, 'rw1PC1'), allData = allData, CI = CI, name = 'rw1PC1'),
        # collect.results(simulationResults = sapply(allResults, `[[`, 'rw1PC2'), allData = allData, CI = CI, name = 'rw1PC2'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'rw2PC1'), allData = allData, CI = CI, name = 'rw2PC1'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'rw2PC2'), allData = allData, CI = CI, name = 'rw2PC2')) %>% 
  dplyr::mutate(model = model %>% factor(., levels = c('crSpline', 'bsSpline', 'psSpline', 'rw1PC1', 'rw1PC2', 'rw2PC1', 'rw2PC2')),
                type = type %>% factor(., levels = c('estimate', 'prediction'), labels = c('Estimation', 'In-sample prediction')))

# results ----

## plots ----

# ggplot2::ggplot(data = final %>% dplyr::filter(period > 2017), aes(x = age)) +
#   ggplot2::geom_line(aes(y = median, group = interaction(model, period),  colour = model), linetype = 'dashed') +
#   ggplot2::geom_line(aes(y = lower, group = interaction(model, period),  colour = model), linetype = 'dashed') +
#   ggplot2::geom_line(aes(y = upper, group = interaction(model, period),  colour = model), linetype = 'dashed') +
#   # ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper, group = interaction(model, period), fill = model),
#   #                      alpha = 0.25, colour = NA) +
#   ggplot2::geom_line(data = final %>% dplyr::filter(period > 2017) %>% dplyr::select(age, period, cohort, true) %>% dplyr::distinct(), 
#                      aes(y = true), linetype = 'dotted') +
#   ggplot2::facet_wrap(~ as.factor(period)) +
#   my.theme(legend.title = element_blank())

maeBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = mae, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dashed') +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Absolute Error') +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank()); maeBP

mseBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = mse, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dashed') +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Square Error') +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank()); mseBP


ggplot2::ggsave(filename = paste0(resultsDir, '/biasBP.png'),
                plot = biasBP,
                width = width, height = height)

ggplot2::ggsave(filename = paste0(resultsDir, '/mseBP.png'),
                plot = mseBP,
                width = width, height = height)

## interval scoring metric ----

# average log-rate for data generated
true_logRate <- log(sapply(allData, `[[`, 'y')/sapply(allData, `[[`, 'N')) %>% rowMeans()
# confidence level
alpha <- 0.05
# upper and lower for each model
## cr spline
crSpline_lower <- final %>% dplyr::filter(model == 'crSpline') %>% dplyr::pull(lower)
crSpline_upper <- final %>% dplyr::filter(model == 'crSpline') %>% dplyr::pull(upper)
## bs spline
bsSpline_lower <- final %>% dplyr::filter(model == 'bsSpline') %>% dplyr::pull(lower)
bsSpline_upper <- final %>% dplyr::filter(model == 'bsSpline') %>% dplyr::pull(upper)
## ps spline
psSpline_lower <- final %>% dplyr::filter(model == 'psSpline') %>% dplyr::pull(lower)
psSpline_upper <- final %>% dplyr::filter(model == 'psSpline') %>% dplyr::pull(upper)
# ## rw1
# rw1PC1_lower <- final %>% dplyr::filter(model == 'rw1PC1') %>% dplyr::pull(lower)
# rw1PC1_upper <- final %>% dplyr::filter(model == 'rw1PC1') %>% dplyr::pull(upper)
# rw1PC2_lower <- final %>% dplyr::filter(model == 'rw1PC2') %>% dplyr::pull(lower)
# rw1PC2_upper <- final %>% dplyr::filter(model == 'rw1PC2') %>% dplyr::pull(upper)
## rw2
rw2PC1_lower <- final %>% dplyr::filter(model == 'rw2PC1') %>% dplyr::pull(lower)
rw2PC1_upper <- final %>% dplyr::filter(model == 'rw2PC1') %>% dplyr::pull(upper)
rw2PC2_lower <- final %>% dplyr::filter(model == 'rw2PC2') %>% dplyr::pull(lower)
rw2PC2_upper <- final %>% dplyr::filter(model == 'rw2PC2') %>% dplyr::pull(upper)

crSplineScore <- interval.score(lower = crSpline_lower, upper = crSpline_upper, true = true_logRate, alpha = 0.05)
bsSplineScore <- interval.score(lower = bsSpline_lower, upper = bsSpline_upper, true = true_logRate, alpha = 0.05)
psSplineScore <- interval.score(lower = psSpline_lower, upper = psSpline_upper, true = true_logRate, alpha = 0.05)
# rw1PC1Score <- interval.score(lower = rw1PC1_lower, upper = rw1PC1_upper, true = true_logRate, alpha = 0.05)
# rw1PC2Score <- interval.score(lower = rw1PC2_lower, upper = rw1PC2_upper, true = true_logRate, alpha = 0.05)
rw2PC1Score <- interval.score(lower = rw2PC1_lower, upper = rw2PC1_upper, true = true_logRate, alpha = 0.05)
rw2PC2Score <- interval.score(lower = rw2PC2_lower, upper = rw2PC2_upper, true = true_logRate, alpha = 0.05)

allScores <- 
  data.frame(model = c('crSpline', 'bsSpline', 'psSpline', 'rw2PC1', 'rw2PC2'),
             averageScore = c(crSplineScore$averageScore,
                              bsSplineScore$averageScore,
                              psSplineScore$averageScore,
                              rw2PC1Score$averageScore,
                              rw2PC2Score$averageScore),
             averageWidth = c(crSplineScore$averageWidth,
                              bsSplineScore$averageWidth,
                              psSplineScore$averageWidth,
                              rw2PC1Score$averageWidth,
                              rw2PC2Score$averageWidth),
             coverage = c(crSplineScore$coverage,
                          bsSplineScore$coverage,
                          psSplineScore$coverage,
                          rw2PC1Score$coverage,
                          rw2PC2Score$coverage))

print(xtable::xtable(allScores), include.rownames = FALSE)
