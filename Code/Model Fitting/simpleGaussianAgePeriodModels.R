# packages ----

library(tidyverse)
library(mgcv)
library(INLA)

# directories ---- 

# # extract file location of this script
# codePath <- rstudioapi::getActiveDocumentContext()$path
# codePathSplitted <- strsplit(codePath, "/")[[1]]
# 
# # retrieve directories
# homeDir <- paste(codePathSplitted[1: (length(codePathSplitted)-3)], collapse = "/")
# codeDir <- paste0(homeDir, '/Code')
# dataDir <- paste0(homeDir, '/Data')
# 
# # need to make sure results is made
# 
# if(!dir.exists(paths = paste0(codeDir, '/Results'))) {
#   dir.create(path = paste0(codeDir, '/Results'))
# }
# 
# resultsDir <- paste0(codeDir, '/Results')

# saving format ----

height <- width <- 10
textSize <- 25

# load ----

## functions ----

# source(paste0(codeDir, '/functions.R'))

# theme for plots
my.theme<-function(...){
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black'),
        # legend.title=element_blank(),
        legend.text.align=0,
        legend.key=element_rect(fill=NA),
        ...)
}

# interval score metrics
interval.score <- function(lower, upper, true, alpha){
  
  # lower = rw2PC1_lower; upper = rw2PC1_upper; true = true_logRate; alpha = 0.05;
  
  # each part seperately
  dispersion <- (upper - lower)
  overprediction <- 2/alpha * (lower - true) * (true < lower)
  underprediction <- 2/alpha * (true - upper) * (true > upper)
  
  averageScore <- (dispersion + underprediction + overprediction) %>% mean()
  averageWidth <- dispersion %>% mean()
  coverage <- ((lower < true & true < upper) %>% sum())/length(true)
  
  return(list(averageScore = averageScore, averageWidth = averageWidth, coverage = coverage))
}

### rewrite data generating function for normal ----

age.fun <- function(a){
  (0.3*a - 0.01*(a^2))
}

per.fun <- function(p){
  (-0.04*p + 0.02*(p^2))
}

## data ----

age1 <- 10:84
M <- 5
ageGroups <- seq(min(age1), max(age1) + M, by = M)
ageInts <- findInterval(age1, ageGroups, all.inside = T)
age <- ((ageGroups[ageInts] + ageGroups[ageInts + 1]) / 2) %>% unique()
period <- 2000:2020

# define data
data <- 
  expand.grid(age = age, 
              period = period) %>%
  dplyr::mutate(mean = age.fun(age - mean(age)) + per.fun(period - mean(period))) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(y = rnorm(n = 1, mean = mean, sd = 1)) %>% 
  dplyr::ungroup()

# models ----

## mgcv ----

dataMGCV.est <- 
  data %>% 
  dplyr::filter(period %in% 2000:2017) %>%
  dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                period_id = period %>% as.factor() %>% as.numeric())

dataMGCV.pred <-
  data %>% 
  dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                period_id = period %>% as.factor() %>% as.numeric())

simpleMGCV <-
  mgcv::gam(formula = 
              y ~ 
              s(age_id, bs = 'cr', k = 12, fx = FALSE) + 
              s(period_id, bs = 'cr', k = 15, fx = FALSE),
            family = 'gaussian', 
            data = dataMGCV.est, 
            method = 'REML',
            control = list(nthreads = parallel::detectCores()))

reparamMGCV <-
  mgcv::gam(formula = 
              y ~ 
              age_id +
              period_id +
              s(age_id, bs = 'cr', k = 12, fx = FALSE) + 
              s(period_id, bs = 'cr', k = 15, fx = FALSE),
            family = 'gaussian', 
            data = dataMGCV.est, 
            method = 'REML',
            control = list(nthreads = parallel::detectCores()))

simpleMGCV.pred <- predict(object = simpleMGCV, newdata = dataMGCV.pred, type = 'link', se.fit = TRUE)
reparamMGCV.pred <- predict(object = reparamMGCV, newdata = dataMGCV.pred, type = 'link', se.fit = TRUE)

simpleAgeLambda <- simpleMGCV$sp[1] %>% as.numeric()
simplePeriodLambda <- simpleMGCV$sp[2] %>% as.numeric()

reparamAgeLambda <- reparamMGCV$sp[1] %>% as.numeric()
reparamPeriodLambda <- reparamMGCV$sp[2] %>% as.numeric()

## inla ----

dataINLA <-
  data %>% 
  dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                period_id = period %>% as.factor() %>% as.numeric(),
                age2_id = age_id,
                period2_id = period_id) %>% 
  dplyr::mutate(y = dplyr::if_else(period > 2017, NA, y)) 

simpleINLA <- 
  INLA::inla(formula = 
               y ~ 
               f(age2_id, model = 'rw2', 
                 hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.05),
                                          initial = log(simpleAgeLambda), fixed = TRUE)),
                 # scale.model = TRUE,
                 constr = TRUE) +
               f(period2_id, model = 'rw2', 
                 hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.05),
                                          initial = log(simplePeriodLambda), fixed = TRUE)),
                 # scale.model = TRUE,
                 constr = TRUE),
             family = 'normal', 
             data = dataINLA,
             control.compute = list(config = TRUE),
             control.predictor = list(compute = FALSE, link = 1),
             control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
             inla.mode = 'experimental')

reparamINLA <- 
  INLA::inla(formula = 
               y ~ 
               age_id + 
               period_id + 
               f(age2_id, model = 'rw2', 
                 hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.05),
                                          initial = log(reparamAgeLambda), fixed = TRUE)),
                 # scale.model = TRUE, 
                 constr = TRUE) +
               f(period2_id, model = 'rw2', 
                 hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.05),
                                          initial = log(reparamPeriodLambda), fixed = TRUE)),
                 # scale.model = TRUE, 
                 constr = TRUE),
             family = 'normal', 
             data = dataINLA,
             control.compute = list(config = TRUE),
             control.predictor = list(compute = FALSE, link = 1),
             control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
             inla.mode = 'experimental')


# results ----

## collect ----

simpleINLA.results <-
  data %>%
  dplyr::mutate(yHat = simpleINLA$summary.linear.predictor$`0.5quant`,
                lower = simpleINLA$summary.linear.predictor$`0.025quant`,
                upper = simpleINLA$summary.linear.predictor$`0.975quant`)

reparamINLA.results <-
  data %>%
  dplyr::mutate(yHat = reparamINLA$summary.linear.predictor$`0.5quant`,
                lower = reparamINLA$summary.linear.predictor$`0.025quant`,
                upper = reparamINLA$summary.linear.predictor$`0.975quant`)

simpleMGCV.results <-
  data %>%
  dplyr::mutate(yHat = simpleMGCV.pred$fit %>% as.vector(),
                se = simpleMGCV.pred$se.fit %>% as.vector()) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lower = yHat + qnorm(p = 0.025, mean = 0, sd = 1) * se,
                upper = yHat + qnorm(p = 0.975, mean = 0, sd = 1) * se) %>%
  dplyr::ungroup()

reparamMGCV.results <-
  data %>%
  dplyr::mutate(yHat = reparamMGCV.pred$fit %>% as.vector(),
                se = reparamMGCV.pred$se.fit %>% as.vector()) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lower = yHat + qnorm(p = 0.025, mean = 0, sd = 1) * se,
                upper = yHat + qnorm(p = 0.975, mean = 0, sd = 1) * se) %>%
  dplyr::ungroup()

results <-
  rbind(simpleINLA.results %>% dplyr::select(-mean, -y) %>% dplyr::mutate(model = 'INLA w/o linear trend'),
        reparamINLA.results %>% dplyr::select(-mean, -y) %>% dplyr::mutate(model = 'INLA w linear trend'),
        simpleMGCV.results %>% dplyr::select(-se, -mean, -y) %>%  dplyr::mutate(model = 'MGCV w/o linear trend'),
        reparamMGCV.results %>% dplyr::select(-se, -mean, -y) %>%  dplyr::mutate(model = 'MGCV w linear trend')) %>%
  dplyr::mutate(model = model %>% factor(., levels = c('INLA w/o linear trend', 
                                                       'INLA w linear trend',
                                                       'MGCV w/o linear trend', 
                                                       'MGCV w linear trend')))

## plot ----

p1 <- 
  ggplot2::ggplot(data = results, aes(x = period)) +
  ggplot2::geom_vline(aes(xintercept = 2017), color = 'red3', linetype = 'dashed', linewidth = 1) +
  ggplot2::geom_line(aes(y = yHat, color = model)) + 
  ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
  ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
  ggplot2::geom_point(data = data, aes(y = mean)) +
  ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
  # ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3', 'pink2')) + 
  ggplot2::labs(y = '', x = 'Year') +
  ggplot2::facet_wrap(~ age, scales = 'free') +
  my.theme(legend.title = element_blank(),
           text = element_text(size = textSize),
           axis.text.x = element_blank(),
           legend.position = 'top'); p1


# ggplot2::ggplot(data = results %>% filter(age%in%c(17.5, 37.5, 57.5, 77.5)),
#                 aes(x = period)) +
#   ggplot2::geom_vline(aes(xintercept = 2017), color = 'red3', linetype = 'dashed', linewidth = 1) +
#   ggplot2::geom_line(aes(y = yHat, color = model)) + 
#   ggplot2:: geom_line(aes(y = lower, color = model), linetype = 'dashed') + 
#   ggplot2:: geom_line(aes(y = upper, color = model), linetype = 'dashed') +
#   ggplot2::geom_point(data = data %>% filter(age%in%c(17.5, 37.5, 57.5, 77.5)), aes(y = mean)) +
#   ggplot2::scale_x_continuous(breaks = seq(2005, 2021, 1)) +
#   # ggplot2::scale_colour_manual(values = c('green4', 'blue3', 'purple3', 'pink2')) + 
#   ggplot2::labs(y = '', x = 'Year') +
#   ggplot2::facet_wrap(~ age, scales = 'free') +
#   my.theme(legend.title = element_blank(),
#            text = element_text(size = textSize),
#            axis.text.x = element_blank(),
#            legend.position = 'top')

ggplot2::ggsave(
  plot = p1,
  filename = 'simpleGaussianAgePeriodModel.png',
  width = 10, height = 10
  )

## scores ----

reparamMGCV.allScore <- interval.score(lower = reparamMGCV.results$lower, upper = reparamMGCV.results$upper, true = data$y, alpha = 0.05)
simpleMGCV.allScore <- interval.score(lower = simpleMGCV.results$lower, upper = simpleMGCV.results$upper, true = data$y, alpha = 0.05)

reparamINLA.allScore <- interval.score(lower = reparamINLA.results$lower, upper = reparamINLA.results$upper, true = data$y, alpha = 0.05)
simpleINLA.allScore <- interval.score(lower = simpleINLA.results$lower, upper = simpleINLA.results$upper, true = data$y, alpha = 0.05)

reparamMGCV.allScore$averageScore
simpleMGCV.allScore$averageScore
reparamINLA.allScore$averageScore
simpleINLA.allScore$averageScore
