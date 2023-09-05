# packages ----

library(tidyverse)
library(mgcv)

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

# data sim arguments
age <- 10:84; N <- 750000; M <- 5
ageGroups <- seq(min(age), max(age) + M, by = M)

set.seed(2468)

data <- 
  data.frame(age1 = age, N = N) %>% 
  dplyr::mutate(
    mean = age.fun(age1 - mean(age1)),
    rate = N * exp(mean)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(y = rpois(n = 1, lambda = rate)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ageInts = findInterval(age1, ageGroups, all.inside = T),
                age = (ageGroups[ageInts] + ageGroups[ageInts + 1]) / 2) %>% 
  dplyr::group_by(age) %>% 
  dplyr::summarise(y = sum(y),
                   N = sum(N)) %>% 
  dplyr::ungroup()

# model fitting ----

## max values ----

A <- data$age %>% unique %>% length()

## models ----

lowKnots_yesPenalty <- 
  gam(y ~ s(age, bs = 'tp', k = 3, fx = FALSE),
      offset = log(N), family = 'poisson', data = data, method = 'REML')

lowKnots_noPenalty <- 
  gam(y ~ s(age, bs = 'tp', k = 3, fx = TRUE),
      offset = log(N), family = 'poisson', data = data, method = 'REML')

highKnots_yesPenalty <- 
  gam(y ~ s(age, bs = 'tp', k = A, fx = FALSE),
      offset = log(N), family = 'poisson', data = data, method = 'REML')

highKnots_noPenalty <- 
  gam(y ~ s(age, bs = 'tp', k = A, fx = TRUE),
      offset = log(N), family = 'poisson', data = data, method = 'REML')

## collect results ----

results <- 
  rbind(data, data, data, data) %>% 
  dplyr::mutate(estimates = c(lowKnots_yesPenalty$linear.predictors, 
                              lowKnots_noPenalty$linear.predictors, 
                              highKnots_yesPenalty$linear.predictors, 
                              highKnots_noPenalty$linear.predictors),
                type = 
                  rep(c('Low Knots & Penalisation', 
                        'Low Knots', 
                        'High Knots & Penalisation', 
                        'High Knots'), each = nrow(data)) %>% 
                  factor(., levels = c('Low Knots & Penalisation', 
                                       'Low Knots', 
                                       'High Knots & Penalisation', 
                                       'High Knots')))

## plot ----

p1 <-
  ggplot2::ggplot(results, aes(x = age)) +
  ggplot2::geom_point(aes(y = log(y)), color = 'black') +
  ggplot2::geom_line(aes(y = estimates, group = type, color = type)) +
  ggplot2::labs(y = 'Linear Predictor', x = 'Age') +
  ggplot2::scale_color_manual(values = c('red3', 'green4', 'blue3', 'purple3')) +
  my.theme(legend.title = element_blank(),
           legend.position = c(0.7, 0.5),
           text = element_text(size = textSize)); p1

## save ----

ggplot2::ggsave(filename = paste0(resultsDir, '/whyIncludePenalty.png'),
                plot = p1,
                height = height, width = width)
