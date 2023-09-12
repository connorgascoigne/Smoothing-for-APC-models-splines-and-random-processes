#functions ----

## spline ----

spline.fit <-
  function(data, predictFrom, mod = c('apc','ac','ap','pc','a','p','c')[1],
           bs = NULL, knots = NULL, fixed = NULL){
    
    dataEst <- 
      data %>% 
      dplyr::filter(period %in% min(data$period):predictFrom)
    
    # if data needs additional changes
    ## center data
    dataEst2 <-
      dataEst %>%
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric())
    
    ak <- knots$age; pk <- knots$period; ck <- knots$cohort
    afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
    
    formula <-
      y ~
      s(age_id, bs = bs, k = ak, fx = afx) +
      s(period_id, bs = bs, k = pk, fx = pfx) +
      s(cohort_id, bs = bs, k = ck, fx = cfx)
    
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in
    # # re-parameterisation
    if (mod == 'ac'){
      formula <- update(formula, ~.- s(period_id, bs = bs, k = pk, fx = pfx))
    } else if (mod == 'pc'){
      formula <- update(formula, ~. - s(age_id, bs = bs, k = ak, fx = afx))
    } else if (mod == 'ap'){
      formula <- update(formula, ~. - s(cohort_id, bs = bs, k = ck, fx = cfx))
    } else if (mod == 'c'){
      formula <- update(formula, ~.- s(age_id, bs = bs, k = ak, fx = afx) - s(period_id, bs = bs, k = pk, fx = pfx))
    } else if (mod == 'a'){
      formula <- update(formula, ~.- s(period_id, bs = bs, k = pk, fx = pfx) - s(cohort_id, bs = bs, k = ck, fx = cfx))
    } else if (mod == 'p'){
      formula <- update(formula, ~.- s(age_id, bs = bs, k = ak, fx = afx) - s(cohort_id, bs = bs, k = ck, fx = cfx))
    }
    
    # use all cores avaliable
    ctrl <- list(nthreads = parallel::detectCores())
    
    # model fit
    fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataEst2, method = 'REML', control = ctrl)
    
    # need to match the naming convention in spline.fit()
    dataPred2 <-
      data %>% 
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric())
    
    prediction <- predict(object = fit, newdata = dataPred2, type = 'link', se.fit = TRUE)
    
    results <-
      data %>%
      dplyr::select(age, period, cohort) %>% 
      dplyr::mutate(yHat = prediction$fit %>% as.vector(),
                    se = prediction$se.fit %>% as.vector()) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(lower = yHat + qnorm(p = 0.025, mean = 0, sd = 1) * se,
                    upper = yHat + qnorm(p = 0.975, mean = 0, sd = 1) * se) %>%
      dplyr::ungroup()
    
    return(results)
    
  }

## random walk ----

randomWalk.fit <-
  function(data, predictFrom, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
           randomWalk = c('rw1', 'rw2')[2],
           pc.u = 1, pc.alpha = 0.01,
           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
           inla.mode = c('classic', 'twostage', 'experimental')[3],
           control.compute = list(config = TRUE), verbose = FALSE, ...){
    
    # data augmentation
    data2 <-
      data %>% 
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric(),
                    age2_id = age_id,
                    period2_id = period_id,
                    cohort2_id = cohort_id) %>% 
      dplyr::mutate(y = dplyr::if_else(period > predictFrom, NA, y))
    
    # define the formula for INLA
    # pc hyper priors
    ## temproal
    hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
    
    # full temporal formula
    formula <-
      y ~ 
      f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE)
    
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in
    # # re-parameterisation
    if (mod == 'ac'){
      formula <- update(formula, ~.
                        - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'pc'){
      formula <- update(formula, ~.
                        - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'ap'){
      formula <- update(formula, ~.
                        - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'c'){
      formula <- update(formula, ~.
                        - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'a'){
      formula <- update(formula, ~.
                        - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if(mod == 'p'){
      formula <- update(formula, ~.
                        - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    }
    
    # model fit
    fit <-
      INLA::inla(formula, family = 'poisson',
                 data = data2, offset = log(data2$N),
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 safe = TRUE,
                 ...)
    
    results <-
      data %>%
      dplyr::mutate(yHat = fit$summary.linear.predictor$`0.5quant`,
                    lower = fit$summary.linear.predictor$`0.025quant`,
                    upper = fit$summary.linear.predictor$`0.975quant`)
    
    return(results)
  }


## data ----

data <- 
  read.csv(paste0(dataDir, '/mortalityExamplesData.csv')) %>% 
  dplyr::mutate(Age = age_lower + 2.5,
                alcohol_rate = log((alcohol_deaths + 0.5)/population),
                selfHarm_rate = log((suicide_deaths+0.5)/population))


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
             mod = 'apc', bs = 'tp', 
             knots = list(age = 10, period = 10, cohort = 12),
             fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw2Fit <- 
  randomWalk.fit(data = alcoholData,
                 predictFrom = 2018,
                 mod = 'apc', randomWalk = 'rw2',
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
             mod = 'apc', bs = 'tp', 
             knots = list(age = 10, period = 10, cohort = 12),
             fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### random walk ----

rw2Fit <- 
  randomWalk.fit(data = selfHarmData,
                 predictFrom = 2018,
                 mod = 'apc', randomWalk = 'rw2',
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






