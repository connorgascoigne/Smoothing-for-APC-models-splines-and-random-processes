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

# A <- data$age %>% unique %>% length()
A <- 6

## models ----

tprs.fit <- 
  gam(y ~ - 1 + s(age, bs = 'tp', k = A, fx = FALSE),
      offset = log(N), family = 'poisson', data = data, method = 'REML', fit = FALSE)

crs.fit <- 
  gam(y ~ - 1 + s(age, bs = 'cr', k = A, fx = FALSE),
      offset = log(N), family = 'poisson', data = data, method = 'REML', fit = FALSE)

bs.fit <- 
  gam(y ~ - 1 + s(age, bs = 'bs', k = A, fx = FALSE),
      offset = log(N), family = 'poisson', data = data, method = 'REML', fit = FALSE)

## collect results ----

tprs.X <- tprs.fit$X %>% as.data.frame()
crs.X <- crs.fit$X %>% as.data.frame()
bs.X <- bs.fit$X %>% as.data.frame()

colnames(tprs.X) <- colnames(crs.X) <- colnames(bs.X) <- paste0('b:', 1:(A-1))

all.X <-
  dplyr::bind_rows(tprs.X %>% 
                     dplyr::mutate(type = 'TPRS',
                                   id = 1:n()) %>% 
                     tidyr::pivot_longer(cols = -c('type', 'id')),
                   crs.X %>% 
                     dplyr::mutate(type = 'CRS',
                                   id = 1:n()) %>% 
                     tidyr::pivot_longer(cols = -c('type', 'id')),
                   bs.X %>% 
                     dplyr::mutate(type = 'BS',
                                   id = 1:n()) %>% 
                     tidyr::pivot_longer(cols = -c('type', 'id'))) %>% 
  dplyr::mutate(type = type %>% factor(., levels = c('CRS', 'BS', 'TPRS')),
                name = name %>% factor(., levels = paste0('b:', 1:(A-1))))

latex.labels <- 
  c(latex2exp::TeX(r'($b_1$)'),
  latex2exp::TeX(r'($b_2$)'),
  latex2exp::TeX(r'($b_3$)'),
  latex2exp::TeX(r'($b_4$)'),
  latex2exp::TeX(r'($b_5$)'))

## plot ----

p1 <-
  ggplot2::ggplot(all.X, aes(x = id, y = value, group = interaction(name, type), colour = name, linetype = name)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(~ type) +
  ggplot2::labs(x = '', y = '') +
  ggplot2::scale_color_manual(values = c('red3', 'green4', 'blue3', 'purple3', 'orange2'),
                              labels = latex.labels) +
  ggplot2::scale_linetype(labels = latex.labels) + 
  my.theme(legend.title = element_blank(),
           axis.text = element_blank(),
           text = element_text(size = textSize)); p1

## save ----

ggplot2::ggsave(filename = paste0(resultsDir, '/differentBasisExamples.png'),
                plot = p1,
                height = height, width = 3*width)
