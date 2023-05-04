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

alcoholData <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = alcohol_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y)

#### prediction ----

alcoholData.splinePredict <- 
  expand.grid(age = alcoholData$age %>% unique() %>% sort(),
              period = min(alcoholData$period):(max(alcoholData$period)+3)) %>% 
  dplyr::mutate(cohort = period - age)

alcoholData.randomWalkPredict <-
  expand.grid(age = alcoholData$age %>% unique() %>% sort(),
              period = min(alcoholData$period):(max(alcoholData$period)+3)) %>%
  dplyr::mutate(cohort = period - age) %>% 
  dplyr::left_join(., alcoholData %>% dplyr::select(age, period, cohort, y, N), by = c('age', 'period', 'cohort'))

### model fits ----

#### spline ----

crSplineFit <- 
  spline.real.fit(dataEst = alcoholData, dataPred = alcoholData.splinePredict,
                   mod = 'apc', slopeDrop = 'c', bs = 'cr', 
                   knots = list(age = 5, period = 5, cohort = 8),
                   fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

bsSplineFit <- 
  spline.real.fit(dataEst = alcoholData, dataPred = alcoholData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'bs', 
                  knots = list(age = 5, period = 5, cohort = 8),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

psSplineFit <- 
  spline.real.fit(dataEst = alcoholData, dataPred = alcoholData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'ps', 
                  knots = list(age = 5, period = 5, cohort = 8),
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
  dplyr::left_join(., alcoholData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')


### plots ----  

#### heat map ----

predictedAlcoholDeathsHeatmap_all <-
  ggplot2::ggplot(alcoholResults, aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  ggplot2::facet_wrap(~ model) +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_all

predictedAlcoholDeathsHeatmap_zoom <-
  ggplot2::ggplot(alcoholResults %>% dplyr::filter(age > 24), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  ggplot2::facet_wrap(~ model) +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsHeatmap_zoom

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_all.png'),
                plot = predictedAlcoholDeathsHeatmap_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedAlcoholDeathsHeatmap_zoom.png'),
                plot = predictedAlcoholDeathsHeatmap_zoom,
                height = height, width = width)


#### line plots ----

predictedAlcoholDeathsLine_all <-
  ggplot2::ggplot(data = alcoholResults, aes(x = period, group = interaction(Age_Group, model), color = model)) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat)) + 
  ggplot2:: geom_line(aes(y = lower), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper), linetype = 'dashed') +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedAlcoholDeathsLine_all

predictedAlcoholDeathsLine_last4 <-
  ggplot2::ggplot(data = alcoholResults %>% dplyr::filter(age > 64), aes(x = period, group = interaction(Age_Group, model), color = model)) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat)) + 
  ggplot2:: geom_line(aes(y = lower), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper), linetype = 'dashed') +
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

selfHarmData <-
  data %>% 
  dplyr::mutate(age = Age, 
                period = Year,
                cohort = period - age,
                N = population,
                y = suicide_deaths) %>% 
  dplyr::select(Age_Group, age, period, cohort, N, y)

#### prediction ----

selfHarmData.splinePredict <- 
  expand.grid(age = selfHarmData$age %>% unique() %>% sort(),
              period = min(selfHarmData$period):(max(selfHarmData$period)+3)) %>% 
  dplyr::mutate(cohort = period - age)

selfHarmData.randomWalkPredict <-
  expand.grid(age = selfHarmData$age %>% unique() %>% sort(),
              period = min(selfHarmData$period):(max(selfHarmData$period)+3)) %>%
  dplyr::mutate(cohort = period - age) %>% 
  dplyr::left_join(., selfHarmData %>% dplyr::select(age, period, cohort, y, N), by = c('age', 'period', 'cohort'))

### model fits ----

#### spline ----

crSplineFit <- 
  spline.real.fit(dataEst = selfHarmData, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'cr', 
                  knots = list(age = 5, period = 5, cohort = 8),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

bsSplineFit <- 
  spline.real.fit(dataEst = selfHarmData, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'bs', 
                  knots = list(age = 5, period = 5, cohort = 8),
                  fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

psSplineFit <- 
  spline.real.fit(dataEst = selfHarmData, dataPred = selfHarmData.splinePredict,
                  mod = 'apc', slopeDrop = 'c', bs = 'ps', 
                  knots = list(age = 5, period = 5, cohort = 8),
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
  dplyr::left_join(., selfHarmData %>% dplyr::select(Age_Group, age) %>% dplyr::distinct(), by = 'age')


### plots ----  

#### heat map ----

predictedSelfHarmDeathsHeatmap_all <-
  ggplot2::ggplot(selfHarmResults, aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(10, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  ggplot2::facet_wrap(~ model) +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_all

predictedSelfHarmDeathsHeatmap_zoom <-
  ggplot2::ggplot(selfHarmResults %>% dplyr::filter(age > 24), aes(y = age, x = period, fill = yHat)) + 
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_viridis_c('Log rate', option = 'G', direction = 1) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::labs(x = 'Year', y = 'Age') +
  ggplot2::facet_wrap(~ model) +
  my.theme(text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsHeatmap_zoom

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_all.png'),
                plot = predictedSelfHarmDeathsHeatmap_all,
                height = height, width = width)

ggplot2::ggsave(filename = paste0(resultsDir, '/predictedSelfHarmDeathsHeatmap_zoom.png'),
                plot = predictedSelfHarmDeathsHeatmap_zoom,
                height = height, width = width)


#### line plots ----

predictedSelfHarmDeathsLine_all <-
  ggplot2::ggplot(data = selfHarmResults, aes(x = period, group = interaction(Age_Group, model), color = model)) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat)) + 
  ggplot2:: geom_line(aes(y = lower), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper), linetype = 'dashed') +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(25, 85, 5)) +
  ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3')) + 
  ggplot2::labs(y = 'Log-Rate', x = 'Period') +
  ggplot2::facet_wrap(~ Age_Group, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)); predictedSelfHarmDeathsLine_all

predictedSelfHarmDeathsLine_last4 <-
  ggplot2::ggplot(data = selfHarmResults %>% dplyr::filter(age > 64), aes(x = period, group = interaction(Age_Group, model), color = model)) +
  ggplot2::geom_vline(aes(xintercept = 2021), color = 'red3', linetype = 'dotted') +
  ggplot2::geom_line(aes(y = yHat)) + 
  ggplot2:: geom_line(aes(y = lower), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper), linetype = 'dashed') +
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






