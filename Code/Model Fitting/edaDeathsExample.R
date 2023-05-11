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

alcoholDeathsPlot_all <-
  ggplot2::ggplot(data, aes(y = Age, x = Year, fill = alcohol_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholDeathsPlot_all

alcoholDeathsPlot_zoom <-
  ggplot2::ggplot(data %>% dplyr::filter(Age > 24), aes(y = Age, x = Year, fill = alcohol_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); alcoholDeathsPlot_zoom

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholDeathsPlot_all.png'),
                plot = alcoholDeathsPlot_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/alcoholDeathsPlot_zoom.png'),
                plot = alcoholDeathsPlot_zoom,
                height = height, width = width)

### self harm ----

selfHarmDeathsPlot_all <-
  ggplot2::ggplot(data, aes(y = Age, x = Year, fill = selfHarm_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmDeathsPlot_all

selfHarmDeathsPlot_zoom <-
  ggplot2::ggplot(data %>% dplyr::filter(Age > 24), aes(y = Age, x = Year, fill = selfHarm_rate)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); selfHarmDeathsPlot_zoom

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmDeathsPlot_all.png'),
                plot = selfHarmDeathsPlot_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/selfHarmDeathsPlot_zoom.png'),
                plot = selfHarmDeathsPlot_zoom,
                height = height, width = width)


# model fits ----

## alcohol ----

### data ----

#### estimation ----

alcoholDataAll <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = alcohol_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y)

alcoholDataReduced <-
  alcoholDataAll %>% 
  dplyr::filter(period %in% min(period):(max(period)-3))

#### prediction ----

alcoholData.splinePredict <- 
  alcoholDataAll %>% 
  dplyr::select(age, period, cohort)

alcoholData.randomWalkPredict <-
  alcoholDataAll %>% 
  dplyr::select(-Age_Group) %>% 
  dplyr::mutate(y = dplyr::if_else(period > (max(period)-3), NA, y))

### model fits ----

#### spline ----

crSplineFit <- 
  spline.real.fit(dataEst = alcoholDataReduced, dataPred = alcoholData.splinePredict,
                   mod = 'apc', slopeDrop = 'c', bs = 'cr', 
                   knots = list(age = 10, period = 10, cohort = 12),
                   fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

bsSplineFit <- 
  spline.real.fit(dataEst = alcoholDataReduced, dataPred = alcoholData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'bs', 
                  knots = list(age = 10, period = 10, cohort = 12),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

psSplineFit <- 
  spline.real.fit(dataEst = alcoholDataReduced, dataPred = alcoholData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'ps', 
                  knots = list(age = 10, period = 10, cohort = 12),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw1FitPC1 <- 
  randomWalk.real.fit(data = alcoholData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                      pc.u = 1, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw2FitPC1 <- 
  randomWalk.real.fit(data = alcoholData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                      pc.u = 1, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw1FitPC2 <- 
  randomWalk.real.fit(data = alcoholData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                      pc.u = 3, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw2FitPC2 <- 
  randomWalk.real.fit(data = alcoholData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                      pc.u = 3, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)


### results collecting ----

# alcoholResults <- 
#   rbind(crSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'CR spline'), 
#         bsSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'BS spline'),
#         psSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'PS spline'),
#         rw1FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1: PC1'),
#         rw2FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2: PC1'),
#         rw1FitPC2 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1: PC2'),
#         rw2FitPC2 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2: PC2')) %>% 
#   dplyr::mutate(model = model %>% factor(., levels = c('CR spline', 'BS spline', 'PS spline',
#                                                        'RW1: PC1', 'RW2: PC1', 'RW1: PC2', 'RW2: PC2'))) %>% 
#   dplyr::left_join(., alcoholData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')

alcoholResults <-
  rbind(crSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'Spline'),
        rw1FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1'),
        rw2FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2')) %>%
  dplyr::mutate(model = model %>% factor(., levels = c('Spline', 'RW1', 'RW2'))) %>%
  dplyr::left_join(., alcoholDataAll %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')

alcoholResultsWide <-
  alcoholResults %>% 
  tidyr::pivot_wider(., names_from = model, values_from = c('yHat', 'lower', 'upper'))


### plots ----  

#### heat map ----

predictedAlcoholDeathsHeatmap_spline_allYears <-
  ggplot2::ggplot(alcoholResults %>% filter(model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_spline_allYears

predictedAlcoholDeathsHeatmap_rw1_allYears <-
  ggplot2::ggplot(alcoholResults %>% filter(model == 'RW1'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_rw1_allYears

predictedAlcoholDeathsHeatmap_rw2_allYears <-
  ggplot2::ggplot(alcoholResults %>% filter(model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_rw2_allYears

predictedAlcoholDeathsHeatmap_spline_25plus <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(age > 24, model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_spline_25plus

predictedAlcoholDeathsHeatmap_rw1_25plus <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(age > 24, model == 'RW1'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_rw1_25plus

predictedAlcoholDeathsHeatmap_rw2_25plus <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(age > 24, model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_rw2_25plus

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_spline_allYears.png'),
                plot = predictedAlcoholDeathsHeatmap_spline_allYears,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_rw1_allYears.png'),
                plot = predictedAlcoholDeathsHeatmap_rw1_allYears,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_rw2_allYears.png'),
                plot = predictedAlcoholDeathsHeatmap_rw2_allYears,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_spline_25plus.png'),
                plot = predictedAlcoholDeathsHeatmap_spline_25plus,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_rw1_25plus.png'),
                plot = predictedAlcoholDeathsHeatmap_rw1_25plus,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_rw2_25plus.png'),
                plot = predictedAlcoholDeathsHeatmap_rw2_25plus,
                height = height, width = width)


#### direct comp ----

splineVsRW1 <-
  ggplot2::ggplot(alcoholResultsWide, aes(x = yHat_Spline, y = yHat_RW1)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); splineVsRW1

splineVsRW2 <-
  ggplot2::ggplot(alcoholResultsWide, aes(x = yHat_Spline, y = yHat_RW2)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); splineVsRW2

RW2VsRW1 <-
  ggplot2::ggplot(alcoholResultsWide, aes(x = yHat_RW2, y = yHat_RW1)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); RW2VsRW1

#### line plots ----

predictedAlcoholDeathsLine_all <-
  ggplot2::ggplot(data = alcoholResults, aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = alcoholDataAll, aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsLine_all

predictedAlcoholDeathsLine_last4 <-
  ggplot2::ggplot(data = alcoholResults %>% dplyr::filter(age > 64), aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = alcoholDataAll %>% dplyr::filter(age > 64), aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsLine_last4

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsLine_all.png'),
                plot = predictedAlcoholDeathsLine_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsLine_last4.png'),
                plot = predictedAlcoholDeathsLine_last4,
                height = height, width = width)



## self harm ----

### data ----

#### estimation ----

selfHarmDataAll <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = suicide_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y)

selfHarmDataReduced <-
  selfHarmDataAll %>% 
  dplyr::filter(period %in% min(period):(max(period)-3))

#### prediction ----

selfHarmData.splinePredict <- 
  selfHarmDataAll %>% 
  dplyr::select(age, period, cohort)

selfHarmData.randomWalkPredict <-
  selfHarmDataAll %>% 
  dplyr::select(-Age_Group) %>% 
  dplyr::mutate(y = dplyr::if_else(period > (max(period)-3), NA, y))

### model fits ----

#### spline ----

crSplineFit <- 
  spline.real.fit(dataEst = selfHarmDataReduced, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'cr', 
                  knots = list(age = 10, period = 10, cohort = 12),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

bsSplineFit <- 
  spline.real.fit(dataEst = selfHarmDataReduced, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'bs', 
                  knots = list(age = 10, period = 10, cohort = 12),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

psSplineFit <- 
  spline.real.fit(dataEst = selfHarmDataReduced, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'ps', 
                  knots = list(age = 10, period = 10, cohort = 12),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw1FitPC1 <- 
  randomWalk.real.fit(data = selfHarmData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                      pc.u = 1, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw2FitPC1 <- 
  randomWalk.real.fit(data = selfHarmData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                      pc.u = 1, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw1FitPC2 <- 
  randomWalk.real.fit(data = selfHarmData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw1',
                      pc.u = 3, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)

rw2FitPC2 <- 
  randomWalk.real.fit(data = selfHarmData.randomWalkPredict,
                      mod = 'apc', slopeDrop = 'c', randomWalk = 'rw2',
                      pc.u = 3, pc.alpha = 0.01,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE)


### results collecting ----

# selfHarmResults <- 
#   rbind(crSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'CR spline'), 
#         bsSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'BS spline'),
#         psSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'PS spline'),
#         rw1FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1: PC1'),
#         rw2FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2: PC1'),
#         rw1FitPC2 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1: PC2'),
#         rw2FitPC2 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2: PC2')) %>% 
#   dplyr::mutate(model = model %>% factor(., levels = c('CR spline', 'BS spline', 'PS spline',
#                                                        'RW1: PC1', 'RW2: PC1', 'RW1: PC2', 'RW2: PC2'))) %>% 
#   dplyr::left_join(., selfHarmData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')

selfHarmResults <-
  rbind(crSplineFit %>% dplyr::select(-se) %>%  dplyr::mutate(model = 'Spline'),
        rw1FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW1'),
        rw2FitPC1 %>% dplyr::select(-y, -N) %>% dplyr::mutate(model = 'RW2')) %>%
  dplyr::mutate(model = model %>% factor(., levels = c('Spline', 'RW1', 'RW2'))) %>%
  dplyr::left_join(., selfHarmDataAll %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')

selfHarmResultsWide <-
  selfHarmResults %>% 
  tidyr::pivot_wider(., names_from = model, values_from = c('yHat', 'lower', 'upper'))

### plots ----  

#### heat map ----

predictedSelfHarmDeathsHeatmap_spline_allYears <-
  ggplot2::ggplot(selfHarmResults %>% filter(model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_spline_allYears

predictedSelfHarmDeathsHeatmap_rw1_allYears <-
  ggplot2::ggplot(selfHarmResults %>% filter(model == 'RW1'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_rw1_allYears

predictedSelfHarmDeathsHeatmap_rw2_allYears <-
  ggplot2::ggplot(selfHarmResults %>% filter(model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_rw2_allYears

predictedSelfHarmDeathsHeatmap_spline_25plus <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(age > 24, model == 'Spline'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_spline_25plus

predictedSelfHarmDeathsHeatmap_rw1_25plus <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(age > 24, model == 'RW1'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_rw1_25plus

predictedSelfHarmDeathsHeatmap_rw2_25plus <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(age > 24, model == 'RW2'), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_rw2_25plus

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_spline_allYears.png'),
                plot = predictedSelfHarmDeathsHeatmap_spline_allYears,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_rw1_allYears.png'),
                plot = predictedSelfHarmDeathsHeatmap_rw1_allYears,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_rw2_allYears.png'),
                plot = predictedSelfHarmDeathsHeatmap_rw2_allYears,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_spline_25plus.png'),
                plot = predictedSelfHarmDeathsHeatmap_spline_25plus,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_rw1_25plus.png'),
                plot = predictedSelfHarmDeathsHeatmap_rw1_25plus,
                height = height, width = width)
ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_rw2_25plus.png'),
                plot = predictedSelfHarmDeathsHeatmap_rw2_25plus,
                height = height, width = width)

#### direct comp ----

splineVsRW1 <-
  ggplot2::ggplot(selfHarmResultsWide, aes(x = yHat_Spline, y = yHat_RW1)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); splineVsRW1

splineVsRW2 <-
  ggplot2::ggplot(selfHarmResultsWide, aes(x = yHat_Spline, y = yHat_RW2)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); splineVsRW2

RW2VsRW1 <-
  ggplot2::ggplot(selfHarmResultsWide, aes(x = yHat_RW2, y = yHat_RW1)) +
  ggplot2::geom_abline(aes(intercept = 0, slope = 1), colour = 'red3', linetype = 'dotted') +
  ggplot2::geom_point() +
  my.theme(); RW2VsRW1


#### line plots ----

predictedSelfHarmDeathsLine_all <-
  ggplot2::ggplot(data = selfHarmResults, aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = selfHarmDataAll, aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsLine_all

predictedSelfHarmDeathsLine_last4 <-
  ggplot2::ggplot(data = selfHarmResults %>% dplyr::filter(age > 64), aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = max(period)-3), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = selfHarmDataAll %>% dplyr::filter(age > 64), aes(y = log(y/N))) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsLine_last4

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsLine_all.png'),
                plot = predictedSelfHarmDeathsLine_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsLine_last4.png'),
                plot = predictedSelfHarmDeathsLine_last4,
                height = height, width = width)






