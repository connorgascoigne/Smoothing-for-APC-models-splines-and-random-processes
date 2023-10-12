# data generating ----

## temporal functions ----

age.fun <- function(a){
  -5*log(10) + (1/50)*(0.3*a - 0.01*(a^2))
}

per.fun <- function(p){
  (1/50)*(-0.04*p + 0.02*(p^2))
}

coh.fun <- function(c){
  (1/50)*(0.35*c - 0.0015*(c^2))
}

## data simulation ----

data.sim <- function(age = NULL, period = NULL, 
                     M, N, 
                     FUNage, FUNperiod, FUNcohort){
  
  # # inputs to run all parts of the function
  # age = 10:84; period = 2000:2020; N = 25;
  # M = 5;
  # FUNage = age.fun;  FUNperiod = per.fun; FUNcohort = coh.fun
  
  # suppress summarise message
  options(dplyr.summarise.inform = FALSE)
  
  # grouped age
  ageGroups <- seq(min(age), max(age) + M, by = M)
  
  # define data
  data <- 
    expand.grid(age1 = age, period = period) %>%
    dplyr::mutate(
      cohort = period - age1,
      mean = 
        FUNage(age1 - mean(age1)) + FUNperiod(period - mean(period)) + FUNcohort(cohort - mean(cohort)),
      N = N,
      rate = N * exp(mean)) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(y = rpois(n = 1, lambda = rate)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(ageInts = findInterval(age1, ageGroups, all.inside = T),
                  age = (ageGroups[ageInts] + ageGroups[ageInts + 1]) / 2) %>% 
    dplyr::group_by(age, period) %>% 
    dplyr::summarise(y = sum(y),
                     N = sum(N)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(cohort = period - age) %>% 
    dplyr::relocate(cohort, .after = period)
  
  return(data)
  
}

# model fits ----


## spline ----

spline.fit <-
  function(data, predictFrom, mod = c('apc','ac','ap','pc','a','p','c')[1],
           bs = NULL, knots = NULL, fixed = NULL, slopeDrop = NULL){
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2020; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # predictFrom = 2017
    # mod = 'apc'; slopeDrop = 'c'; bs = 'cr'
    # knots = list(age = 5, period = 5, cohort = 8)
    # fixed = list(age = F,period = F,cohort = F)
    
    if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
    if(is.null(bs)) stop('Warning: Need a basis function.')
    if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
    if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
    if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
    
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
      age_id + period_id + cohort_id +
      s(age_id, bs = bs, k = ak, fx = afx) +
      s(period_id, bs = bs, k = pk, fx = pfx) +
      s(cohort_id, bs = bs, k = ck, fx = cfx)
    
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in
    # # re-parameterisation
    if (mod == 'apc'){
      if(slopeDrop == 'a'){
        formula <- update(formula, ~. - age_id)
      } else if (slopeDrop == 'p'){
        formula <- update(formula, ~. - period_id)
      } else if (slopeDrop == 'c'){
        formula <- update(formula, ~. - cohort_id)
      }
    } else if (mod == 'ac'){
      formula <- update(formula, ~.- period_id - s(period_id, bs = bs, k = pk, fx = pfx))
    } else if (mod == 'pc'){
      formula <- update(formula, ~. - age_id - s(age_id, bs = bs, k = ak, fx = afx))
    } else if (mod == 'ap'){
      formula <- update(formula, ~. - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
    } else if (mod == 'c'){
      formula <- update(formula, ~.- age_id - s(age_id, bs = bs, k = ak, fx = afx) - period_id - s(period_id, bs = bs, k = pk, fx = pfx))
    } else if (mod == 'a'){
      formula <- update(formula, ~.- period_id - s(period_id, bs = bs, k = pk, fx = pfx) - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
    } else if (mod == 'p'){
      formula <- update(formula, ~.- age_id - s(age_id, bs = bs, k = ak, fx = afx) - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
    }
    
    # use all cores avaliable
    ctrl <- list(nthreads = parallel::detectCores())
    
    # model fit
    fit <- mgcv::gam(formula, 
                     offset = log(N), 
                     family = 'poisson', 
                     data = dataEst2, 
                     method = 'REML', 
                     control = ctrl)
    
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
           slopeDrop = NULL,
           pc.u = 1, pc.alpha = 0.01,
           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
           inla.mode = c('classic', 'twostage', 'experimental')[3],
           control.compute = list(config = TRUE), verbose = FALSE, ...){
    
    
    # # inputs to run all parts of the function
    # data = alcoholData.randomWalkPredict
    # mod = 'apc'; slopeDrop = 'c'; randomWalk = 'rw1';
    # pc.u = 1; pc.alpha = 0.01;
    # control.inla = list(strategy = 'adaptive', int.strategy = 'auto');
    # inla.mode = c('classic', 'twostage', 'experimental')[3];
    # control.compute = list(config = TRUE); verbose = FALSE
    
    if (!"config" %in% names(control.compute)) {
      message("config = TRUE is added to control.compute so that posterior draws can be taken.")
      control.compute$config <- TRUE
    }
    if(mod == 'apc' && is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
    if(mod == 'apc' &&! (slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
    
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
      y ~ age_id + period_id + cohort_id +
      f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE)
    
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in
    # # re-parameterisation
    if (mod == 'apc'){
      if(slopeDrop == 'a'){
        formula <- update(formula, ~. - age_id)
      } else if (slopeDrop == 'p'){
        formula <- update(formula, ~. - period_id)
      } else if (slopeDrop == 'c'){
        formula <- update(formula, ~. - cohort_id)
      }
    } else if (mod == 'ac'){
      formula <- update(formula, ~.
                        - period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'pc'){
      formula <- update(formula, ~.
                        - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'ap'){
      formula <- update(formula, ~.
                        - cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'c'){
      formula <- update(formula, ~.
                        - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if (mod == 'a'){
      formula <- update(formula, ~.
                        - period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    } else if(mod == 'p'){
      formula <- update(formula, ~.
                        - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
                          cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
    }
    
    # model fit
    fit <-
      INLA::inla(formula, family = 'poisson',
                 data = data2, 
                 offset = log(data2$N),
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

# additional functions ----

# summary for estimates
my.summary <- function(x, CI = 0.95) {
  
  # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
  
  # lower and upper CI calculators
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}

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
  coverage <- ((lower < true & true < upper) %>% sum)/length(true)
  
  return(list(averageScore = averageScore, averageWidth = averageWidth, coverage = coverage))
}

# mae, mse, is, width and coverage for one set of data
model.score.summary <- function(results, trueData, CI, periods, model) {
  
  # results = splineResults; trueData = trueData2; CI = CI; periods = periods; model = 'Spline'
  
  res <-
    results %>% 
    dplyr::left_join(., trueData, by = c('age', 'period', 'cohort')) %>% 
    dplyr::filter(period %in% periods)
  
  periodLabel <- paste0(min(periods),':',max(periods))
  
  truth <- res$true
  estimate <- res$yHat
  lower <- res$lower
  upper <- res$upper
  
  
  scores <- interval.score(lower = lower, upper = upper, true = truth, alpha = (1-CI))
  
  mae <- data.frame(periods = periodLabel, model = model, mae = abs(estimate - truth))
  mse <- data.frame(periods = periodLabel, model = model, mse = (estimate - truth)^2)
  is <- data.frame(periods = periodLabel, model = model, is = scores$averageScore)
  width <- data.frame(periods = periodLabel, model = model, width = scores$averageWidth)
  coverage <- data.frame(periods = periodLabel, model = model, coverage = scores$coverage)
  
  return(list(mae = mae, mse = mse, is = is, width = width, coverage = coverage))
  
}

# scores for all simulation data models at once
collect.simulation.results <- function(allBasisResults, trueData, CI, periods){
  
  # allBasisResults = allResults[[i]]; trueData = trueData; CI = CI; periods = 2000:2017
  
  crSplineResults <- model.score.summary(results = allBasisResults$crSpline, trueData = trueData, CI = CI, periods = periods, model = 'crSpline')
  bsSplineResults <- model.score.summary(results = allBasisResults$bsSpline, trueData = trueData, CI = CI, periods = periods, model = 'bsSpline')
  tprsSplineResults <- model.score.summary(results = allBasisResults$tpSpline, trueData = trueData, CI = CI, periods = periods, model = 'tpSpline')
  rw2PC1Results <- model.score.summary(results = allBasisResults$rw2PC1, trueData = trueData, CI = CI, periods = periods, model = 'rw2PC1')
  rw2PC2Results <- model.score.summary(results = allBasisResults$rw2PC2, trueData = trueData, CI = CI, periods = periods, model = 'rw2PC2')
  rw2PC3Results <- model.score.summary(results = allBasisResults$rw2PC3, trueData = trueData, CI = CI, periods = periods, model = 'rw2PC3')
  
  scores <- 
    rbind(crSplineResults$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          bsSplineResults$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          tprsSplineResults$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          rw2PC1Results$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          rw2PC2Results$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          rw2PC3Results$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae') %>% relocate(metric, .after = model),
          crSplineResults$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          bsSplineResults$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          tprsSplineResults$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          rw2PC1Results$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          rw2PC2Results$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          rw2PC3Results$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse') %>% relocate(metric, .after = model),
          crSplineResults$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          bsSplineResults$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          tprsSplineResults$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          rw2PC1Results$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          rw2PC2Results$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          rw2PC3Results$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is') %>% relocate(metric, .after = model),
          crSplineResults$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          bsSplineResults$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          tprsSplineResults$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          rw2PC1Results$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          rw2PC2Results$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          rw2PC3Results$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width') %>% relocate(metric, .after = model),
          crSplineResults$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model),
          bsSplineResults$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model),
          tprsSplineResults$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model),
          rw2PC1Results$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model),
          rw2PC2Results$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model),
          rw2PC3Results$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage') %>% relocate(metric, .after = model))
  
  scores
  
}

# scores for all suicide data models at once
collect.suicide.results <- function(allModelResults, trueData, CI, periods){
  
  # allModelResults = alcoholResults; trueData = alcoholResults; CI = 0.95; periods = 2000:2017
  
  trueData2 <- 
    trueData %>% 
    dplyr::mutate(true = log_rate) %>% 
    dplyr::select(age, period, cohort, true)
  
  splineResults <-
    allModelResults %>% 
    dplyr::filter(model == 'Spline')
  
  rw2Results <-
    allModelResults %>% 
    dplyr::filter(model == 'RW2')
  
  splineScores <- model.score.summary(results = splineResults, trueData = trueData2, CI = CI, periods = periods, model = 'Spline')
  rw2Scores <- model.score.summary(results = rw2Results, trueData = trueData2, CI = CI, periods = periods, model = 'RW2')
  
  
  scores <- 
    rbind(splineScores$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae'),
          rw2Scores$mae %>% dplyr::rename(score = mae) %>% dplyr::mutate(metric = 'mae'),
          splineScores$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse'),
          rw2Scores$mse %>% dplyr::rename(score = mse) %>% dplyr::mutate(metric = 'mse'),
          splineScores$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is'),
          rw2Scores$is %>% dplyr::rename(score = is) %>% dplyr::mutate(metric = 'is'),
          splineScores$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width'),
          rw2Scores$width %>% dplyr::rename(score = width) %>% dplyr::mutate(metric = 'width'),
          splineScores$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage'),
          rw2Scores$coverage %>% dplyr::rename(score = coverage) %>% dplyr::mutate(metric = 'coverage'))
  
  
  scores
  
}
