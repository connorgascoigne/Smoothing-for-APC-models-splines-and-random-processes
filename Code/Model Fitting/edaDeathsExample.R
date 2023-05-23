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

data <- 
  read.csv(paste0(dataDir, '/mortalityExamplesData.csv')) %>% 
  dplyr::mutate(Age = age_lower + 2.5,
                alcohol_rate = log((alcohol_deaths + 0.5)/population),
                selfHarm_rate = log((suicide_deaths+0.5)/population))


# data exploration ----

## heat maps ----

### alcohol ----

alcoholObservedHeatmap <-
  ggplot2::ggplot(data, aes(y = Age, x = Year, fill = alcohol_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholObservedHeatmap

alcoholObservedHeatmap_25plus <-
  ggplot2::ggplot(data %>% dplyr::filter(Age > 24), aes(y = Age, x = Year, fill = alcohol_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholObservedHeatmap_25plus

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholObservedHeatmap.png'),
                plot = alcoholObservedHeatmap,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholObservedHeatmap_25plus.png'),
                plot = alcoholObservedHeatmap_25plus,
                height = height, width = width)

### self harm ----

selfHarmObservedHeatmap <-
  ggplot2::ggplot(data, aes(y = Age, x = Year, fill = selfHarm_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmObservedHeatmap

selfHarmObservedHeatmap_25plus <-
  ggplot2::ggplot(data %>% dplyr::filter(Age > 24), aes(y = Age, x = Year, fill = selfHarm_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmObservedHeatmap_25plus

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmObservedHeatmap.png'),
                plot = selfHarmObservedHeatmap,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmObservedHeatmap_25plus.png'),
                plot = selfHarmObservedHeatmap_25plus,
                height = height, width = width)


# model fits ----

## alcohol ----

### data ----

alcoholData <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = alcohol_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y) %>% 
  dplyr::filter(age >= 25)

### model fits ----

#### spline ----

splineFit <- 
  spline.fit(data = alcoholData,
             predictFrom = 2018,
             mod = 'apc', slopeDrop = 'c', bs = 'tp', 
             knots = list(age = 10, period = 10, cohort = 12),
             fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw2Fit <- 
  randomWalk.fit(data = alcoholData,
                 predictFrom = 2018,
                 mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                 pc.u = 1, pc.alpha = 0.01,
                 control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                 inla.mode = c('classic', 'twostage', 'experimental')[3],
                 control.compute = list(config = TRUE), verbose = FALSE)


### results collecting ----

alcoholResults <-
  rbind(splineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'Spline'),
        rw2Fit %>% dplyr::select(-y, -N, -Age_Group) %>% dplyr::mutate(model = 'RW2')) %>%
  dplyr::mutate(model = model %>% factor(., levels = c('Spline', 'RW2'))) %>%
  dplyr::left_join(., alcoholData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')


### plots ----  

#### heat map ----

alcoholPredictedHeatmap_spline <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholPredictedHeatmap_spline

alcoholPredictedHeatmap_rw2 <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholPredictedHeatmap_rw2

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholPredictedHeatmap_spline.png'),
                plot = alcoholPredictedHeatmap_spline,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholPredictedHeatmap_rw2.png'),
                plot = alcoholPredictedHeatmap_rw2,
                height = height, width = width)


#### direct comp ----

alcoholPredicted_splineVsRW2 <-
  ggplot2::ggplot(alcoholResults %>% 
                    tidyr::pivot_wider(., names_from = model, values_from = c('yHat', 'lower', 'upper')) %>% 
                    dplyr::mutate(type = dplyr::if_else(period <= (max(period) -3), 'Estimation', 'In-sample prediction')),
                  aes(x = yHat_Spline, y = yHat_RW2, shape = type, colour = type)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_point() +
  ggplot2::scale_shape_manual(values = c(1, 2)) +
  ggplot2::scale_color_manual(values = c('green4', 'blue3')) +
  ggplot2::labs(x = 'Spline', y = 'RW2') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           legend.position = c(0.85, 0.1)); alcoholPredicted_splineVsRW2

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholPredicted_splineVsRW2.png'),
                plot = alcoholPredicted_splineVsRW2,
                height = height, width = width)

#### line plots ----

alcoholPredictedLineplot <-
  ggplot2::ggplot(data = 
                    alcoholResults %>% 
                    dplyr::mutate(Age_Group = Age_Group %>% stringr::str_replace(., 'Aged ', '') %>% stringr::str_replace(., ' years', '')), 
                  aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = alcoholData %>% 
                        dplyr::mutate(Age_Group = Age_Group %>% stringr::str_replace(., 'Aged ', '') %>% stringr::str_replace(., ' years', '')), 
                      aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Year') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_blank()); alcoholPredictedLineplot

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholPredictedLineplot.png'),
                plot = alcoholPredictedLineplot,
                height = height, width = width)


## self harm ----

### data ----

selfHarmData <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = suicide_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y) %>% 
  dplyr::filter(age >= 25)

### model fits ----

#### spline ----

splineFit <- 
  spline.fit(data = selfHarmData,
             predictFrom = 2018,
             mod = 'apc', slopeDrop = 'c', bs = 'tp', 
             knots = list(age = 10, period = 10, cohort = 12),
             fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw2Fit <- 
  randomWalk.fit(data = selfHarmData,
                 predictFrom = 2018,
                 mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                 pc.u = 1, pc.alpha = 0.01,
                 control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                 inla.mode = c('classic', 'twostage', 'experimental')[3],
                 control.compute = list(config = TRUE), verbose = FALSE)


### results collecting ----

selfHarmResults <-
  rbind(splineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'Spline'),
        rw2Fit %>% dplyr::select(-y, -N, -Age_Group) %>% dplyr::mutate(model = 'RW2')) %>%
  dplyr::mutate(model = model %>% factor(., levels = c('Spline', 'RW2'))) %>%
  dplyr::left_join(., selfHarmData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')

### plots ----  

#### heat map ----

selfHarmPredictedHeatmap_spline <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmPredictedHeatmap_spline

selfHarmPredictedHeatmap_rw2 <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmPredictedHeatmap_rw2

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmPredictedHeatmap_spline.png'),
                plot = selfHarmPredictedHeatmap_spline,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmPredictedHeatmap_rw2.png'),
                plot = selfHarmPredictedHeatmap_rw2,
                height = height, width = width)


#### direct comp ----

selfHarmPredicted_splineVsRW2 <-
  ggplot2::ggplot(selfHarmResults %>% 
                    tidyr::pivot_wider(., names_from = model, values_from = c('yHat', 'lower', 'upper')) %>% 
                    dplyr::mutate(type = dplyr::if_else(period <= (max(period) -3), 'Estimation', 'In-sample prediction')),
                  aes(x = yHat_Spline, y = yHat_RW2, shape = type, colour = type)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_point() +
  ggplot2::scale_shape_manual(values = c(1, 2)) +
  ggplot2::scale_color_manual(values = c('green4', 'blue3')) +
  ggplot2::labs(x = 'Spline', y = 'RW2') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           legend.position = c(0.85, 0.1)); alcoholPredicted_splineVsRW2

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmPredicted_splineVsRW2.png'),
                plot = selfHarmPredicted_splineVsRW2,
                height = height, width = width)

#### line plots ----

selfHarmPredictedLineplot <-
  ggplot2::ggplot(data = 
                    selfHarmResults %>% 
                    dplyr::mutate(Age_Group = Age_Group %>% stringr::str_replace(., 'Aged ', '') %>% stringr::str_replace(., ' years', '')), 
                  aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted', linewidth = 1) +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = 
                        selfHarmData %>% 
                        dplyr::mutate(Age_Group = Age_Group %>% stringr::str_replace(., 'Aged ', '') %>% stringr::str_replace(., ' years', '')), 
                      aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Year') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_blank()); selfHarmPredictedLineplot

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmPredictedLineplot.png'),
                plot = selfHarmPredictedLineplot,
                height = height, width = width)

# models scores ----

alcoholSplineScore <- find.score(true = alcoholData, results = alcoholResults, model = 'Spline', predictFrom = 2017, CI = 0.95)
alcoholRW2Score <- find.score(true = alcoholData, results = alcoholResults, model = 'RW2', predictFrom = 2017, CI = 0.95)

selfHarmSplineScore <- find.score(true = selfHarmData, results = selfHarmResults, model = 'Spline', predictFrom = 2017, CI = 0.95)
selfHarmRW2Score <- find.score(true = selfHarmData, results = selfHarmResults, model = 'RW2', predictFrom = 2017, CI = 0.95)


estimateScores <- 
  data.frame(data = rep(c('Alcohol', 'Self Harm'), each = 2),
             model = rep(c('Spline', 'RW2'), times = 2),
             averageScore = c(alcoholSplineScore$scoreEstimate$averageScore,
                              alcoholRW2Score$scoreEstimate$averageScore,
                              selfHarmSplineScore$scoreEstimate$averageScore,
                              selfHarmRW2Score$scoreEstimate$averageScore),
             averageWidth = c(alcoholSplineScore$scoreEstimate$averageWidth,
                              alcoholRW2Score$scoreEstimate$averageWidth,
                              selfHarmSplineScore$scoreEstimate$averageWidth,
                              selfHarmRW2Score$scoreEstimate$averageWidth),
             coverage = c(alcoholSplineScore$scoreEstimate$coverage,
                          alcoholRW2Score$scoreEstimate$coverage,
                          selfHarmSplineScore$scoreEstimate$coverage,
                          selfHarmRW2Score$scoreEstimate$coverage)) %>% 
  dplyr::mutate(across(-c('model', 'data'), ~ round(.x*100, digits = 2))); estimateScores


predictScores <- 
  data.frame(data = rep(c('Alcohol', 'Self Harm'), each = 2),
             model = rep(c('Spline', 'RW2'), times = 2),
             averageScore = c(alcoholSplineScore$scorePredict$averageScore,
                              alcoholRW2Score$scorePredict$averageScore,
                              selfHarmSplineScore$scorePredict$averageScore,
                              selfHarmRW2Score$scorePredict$averageScore),
             averageWidth = c(alcoholSplineScore$scorePredict$averageWidth,
                              alcoholRW2Score$scorePredict$averageWidth,
                              selfHarmSplineScore$scorePredict$averageWidth,
                              selfHarmRW2Score$scorePredict$averageWidth),
             coverage = c(alcoholSplineScore$scorePredict$coverage,
                          alcoholRW2Score$scorePredict$coverage,
                          selfHarmSplineScore$scorePredict$coverage,
                          selfHarmRW2Score$scorePredict$coverage)) %>% 
  dplyr::mutate(across(-c('model', 'data'), ~ round(.x*100, digits = 2))); predictScores


allScores <- 
  cbind(estimateScores %>% dplyr::mutate(type = 'Estimation'),
        predictScores %>% dplyr::mutate(type = 'In-sample prediction'))
allScores[,c(1:5,9:11)]
print(xtable::xtable(allScores[,c(1:5,9:11)]), include.rownames = FALSE)





