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
  # list packages that are used in loop
  .packages = c('mgcv', 'INLA', 'tidyverse')
) %dopar% {
  # cr spline basis
  crSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'cr',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  # bs spline basis
  bsSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'bs',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  # tp spline basis
  tpSpline <- spline.fit(data = allData[[i]], predictFrom = 2017,
                         mod = 'apc', slopeDrop = 'c', bs = 'tp',
                         knots = list(age = 10, period = 10, cohort = 12),
                         fixed = list(age = F, period = F, cohort = F))$yHat
  # rw2 U = 1
  rw2PC1 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 1, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)$yHat
  # rw2 U = 3
  rw2PC2 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 3, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)$yHat
  
  # rw2 U = 6
  rw2PC3 <- randomWalk.fit(data = allData[[i]], predictFrom = 2017,
                           mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                           pc.u = 6, pc.alpha = 0.01,
                           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                           inla.mode = c('classic', 'twostage', 'experimental')[3],
                           control.compute = list(config = TRUE), verbose = FALSE)$yHat
  
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

CI <- 0.95

final <- 
  rbind(collect.results(simulationResults = sapply(allResults, `[[`, 'crSpline'), allData = allData, CI = CI, name = 'crSpline'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'bsSpline'), allData = allData, CI = CI, name = 'bsSpline'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'tpSpline'), allData = allData, CI = CI, name = 'tpSpline'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'rw2PC1'), allData = allData, CI = CI, name = 'rw2PC1'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'rw2PC2'), allData = allData, CI = CI, name = 'rw2PC2'),
        collect.results(simulationResults = sapply(allResults, `[[`, 'rw2PC3'), allData = allData, CI = CI, name = 'rw2PC3')) %>% 
  dplyr::mutate(model = model %>% factor(., levels = c('crSpline', 'bsSpline', 'tpSpline', 'rw2PC1', 'rw2PC2', 'rw2PC3')),
                type = type %>% factor(., levels = c('estimate', 'prediction'), labels = c('Estimation', 'In-sample prediction')))

# results ----

## plots ----

maeBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = mae, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Absolute Error') +
  ggplot2::facet_grid(~ type) +
  my.theme(text = element_text(size = textSize),
           legend.title = element_blank(),
           axis.text.x = element_blank()); maeBP

mseBP <- 
  ggplot2::ggplot(data = final, aes(x = model, y = mse, col = model, fill = model)) +
  ggplot2::geom_hline(aes(yintercept = 0), colour = 'black', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_boxplot(alpha = 0.2) +
  ggplot2::labs(x = '', y = 'Mean Square Error') +
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

## interval scoring metric ----

# average log-rate for data generated
trueData <- 
  expand.grid(age = seq(from = 12.5, to = 82.5, by = 5), 
              period = 2000:2020) %>%
  dplyr::mutate(cohort = period - age,
                mean = age.fun(age - mean(age)) + 
                  per.fun(period - mean(period)) + 
                  coh.fun(cohort - mean(cohort))) %>% 
  # need to arrange to match order of the data!!!!
  dplyr::arrange(age, period)
true_logRate_estimate <- trueData %>% dplyr::filter(period <= 2017) %>% dplyr::pull(mean)
true_logRate_inSamplePredcition <- trueData %>% dplyr::filter(period > 2017) %>% dplyr::pull(mean)

# estimation
## cr spline
crSpline_lower_estimate <- final %>% dplyr::filter(model == 'crSpline', type == 'Estimation') %>% dplyr::pull(lower)
crSpline_upper_estimate <- final %>% dplyr::filter(model == 'crSpline', type == 'Estimation') %>% dplyr::pull(upper)
## bs spline
bsSpline_lower_estimate <- final %>% dplyr::filter(model == 'bsSpline', type == 'Estimation') %>% dplyr::pull(lower)
bsSpline_upper_estimate <- final %>% dplyr::filter(model == 'bsSpline', type == 'Estimation') %>% dplyr::pull(upper)
## tp spline
tpSpline_lower_estimate <- final %>% dplyr::filter(model == 'tpSpline', type == 'Estimation') %>% dplyr::pull(lower)
tpSpline_upper_estimate <- final %>% dplyr::filter(model == 'tpSpline', type == 'Estimation') %>% dplyr::pull(upper)
## rw2 U = 1
rw2PC1_lower_estimate <- final %>% dplyr::filter(model == 'rw2PC1', type == 'Estimation') %>% dplyr::pull(lower)
rw2PC1_upper_estimate <- final %>% dplyr::filter(model == 'rw2PC1', type == 'Estimation') %>% dplyr::pull(upper)
## rw2 U = 3
rw2PC2_lower_estimate <- final %>% dplyr::filter(model == 'rw2PC2', type == 'Estimation') %>% dplyr::pull(lower)
rw2PC2_upper_estimate <- final %>% dplyr::filter(model == 'rw2PC2', type == 'Estimation') %>% dplyr::pull(upper)
## rw2 U = 6
rw2PC3_lower_estimate <- final %>% dplyr::filter(model == 'rw2PC3', type == 'Estimation') %>% dplyr::pull(lower)
rw2PC3_upper_estimate <- final %>% dplyr::filter(model == 'rw2PC3', type == 'Estimation') %>% dplyr::pull(upper)
## scores
crSplineScore_estimate <- interval.score(lower = crSpline_lower_estimate, upper = crSpline_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))
bsSplineScore_estimate <- interval.score(lower = bsSpline_lower_estimate, upper = bsSpline_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))
tpSplineScore_estimate <- interval.score(lower = tpSpline_lower_estimate, upper = tpSpline_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))
rw2PC1Score_estimate <- interval.score(lower = rw2PC1_lower_estimate, upper = rw2PC1_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))
rw2PC2Score_estimate <- interval.score(lower = rw2PC2_lower_estimate, upper = rw2PC2_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))
rw2PC3Score_estimate <- interval.score(lower = rw2PC3_lower_estimate, upper = rw2PC3_upper_estimate, true = true_logRate_estimate, alpha = (1-CI))

# in-sample prediction
## cr spline
crSpline_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'crSpline', type == 'In-sample prediction') %>% dplyr::pull(lower)
crSpline_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'crSpline', type == 'In-sample prediction') %>% dplyr::pull(upper)
## bs spline
bsSpline_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'bsSpline', type == 'In-sample prediction') %>% dplyr::pull(lower)
bsSpline_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'bsSpline', type == 'In-sample prediction') %>% dplyr::pull(upper)
## tp spline
tpSpline_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'tpSpline', type == 'In-sample prediction') %>% dplyr::pull(lower)
tpSpline_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'tpSpline', type == 'In-sample prediction') %>% dplyr::pull(upper)
## rw2 U = 1
rw2PC1_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC1', type == 'In-sample prediction') %>% dplyr::pull(lower)
rw2PC1_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC1', type == 'In-sample prediction') %>% dplyr::pull(upper)
## rw2 U = 3
rw2PC2_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC2', type == 'In-sample prediction') %>% dplyr::pull(lower)
rw2PC2_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC2', type == 'In-sample prediction') %>% dplyr::pull(upper)
## rw2 U = 6
rw2PC3_lower_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC3', type == 'In-sample prediction') %>% dplyr::pull(lower)
rw2PC3_upper_inSamplePredcition <- final %>% dplyr::filter(model == 'rw2PC3', type == 'In-sample prediction') %>% dplyr::pull(upper)
## scores
crSplineScore_inSamplePredcition <- interval.score(lower = crSpline_lower_inSamplePredcition, upper = crSpline_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))
bsSplineScore_inSamplePredcition <- interval.score(lower = bsSpline_lower_inSamplePredcition, upper = bsSpline_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))
tpSplineScore_inSamplePredcition <- interval.score(lower = tpSpline_lower_inSamplePredcition, upper = tpSpline_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))
rw2PC1Score_inSamplePredcition <- interval.score(lower = rw2PC1_lower_inSamplePredcition, upper = rw2PC1_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))
rw2PC2Score_inSamplePredcition <- interval.score(lower = rw2PC2_lower_inSamplePredcition, upper = rw2PC2_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))
rw2PC3Score_inSamplePredcition <- interval.score(lower = rw2PC3_lower_inSamplePredcition, upper = rw2PC3_upper_inSamplePredcition, true = true_logRate_inSamplePredcition, alpha = (1-CI))



estimateScores <- 
  data.frame(model = c('crSpline', 'bsSpline', 'tpSpline', 'rw2PC1', 'rw2PC2', 'rw2PC3'),
             averageScore = c(crSplineScore_estimate$averageScore,
                              bsSplineScore_estimate$averageScore,
                              tpSplineScore_estimate$averageScore,
                              rw2PC1Score_estimate$averageScore,
                              rw2PC2Score_estimate$averageScore,
                              rw2PC3Score_estimate$averageScore),
             averageWidth = c(crSplineScore_estimate$averageWidth,
                              bsSplineScore_estimate$averageWidth,
                              tpSplineScore_estimate$averageWidth,
                              rw2PC1Score_estimate$averageWidth,
                              rw2PC2Score_estimate$averageWidth,
                              rw2PC3Score_estimate$averageWidth),
             coverage = c(crSplineScore_estimate$coverage,
                          bsSplineScore_estimate$coverage,
                          tpSplineScore_estimate$coverage,
                          rw2PC1Score_estimate$coverage,
                          rw2PC2Score_estimate$coverage,
                          rw2PC3Score_estimate$coverage)) %>% 
  dplyr::mutate(across(-model, ~ round(.x * 100, digits = 2)))


inSamplePredictionScores <- 
  data.frame(model = c('crSpline', 'bsSpline', 'tpSpline', 'rw2PC1', 'rw2PC2', 'rw2PC3'),
             averageScore = c(crSplineScore_inSamplePredcition$averageScore,
                              bsSplineScore_inSamplePredcition$averageScore,
                              tpSplineScore_inSamplePredcition$averageScore,
                              rw2PC1Score_inSamplePredcition$averageScore,
                              rw2PC2Score_inSamplePredcition$averageScore,
                              rw2PC3Score_inSamplePredcition$averageScore),
             averageWidth = c(crSplineScore_inSamplePredcition$averageWidth,
                              bsSplineScore_inSamplePredcition$averageWidth,
                              tpSplineScore_inSamplePredcition$averageWidth,
                              rw2PC1Score_inSamplePredcition$averageWidth,
                              rw2PC2Score_inSamplePredcition$averageWidth,
                              rw2PC3Score_inSamplePredcition$averageWidth),
             coverage = c(crSplineScore_inSamplePredcition$coverage,
                          bsSplineScore_inSamplePredcition$coverage,
                          tpSplineScore_inSamplePredcition$coverage,
                          rw2PC1Score_inSamplePredcition$coverage,
                          rw2PC2Score_inSamplePredcition$coverage,
                          rw2PC3Score_inSamplePredcition$coverage)) %>% 
  dplyr::mutate(across(-model, ~ round(.x * 100, digits = 2)))


allScores <- 
  cbind(estimateScores %>% dplyr::mutate(type = 'Estimation'),
        inSamplePredictionScores %>% dplyr::mutate(type = 'In-sample prediction'))
allScores[,c(1:4,7:9)]
print(xtable::xtable(allScores[,c(1:4,7:9)]), include.rownames = FALSE)
