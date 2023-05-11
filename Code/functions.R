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

## spline model ----

spline.fit <- 
  function(data, mod = c('apc','ac','ap','pc','a','p','c')[1],
           bs = NULL, knots = NULL, fixed = NULL, slopeDrop = NULL){
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2020; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # mod = 'apc'; slopeDrop = 'c'; bs = 'cr'
    # knots = list(age = 5, period = 5, cohort = 8)
    # fixed = list(age = F,period = F,cohort = F)
    
    if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
    if(is.null(bs)) stop('Warning: Need a basis function.')
    if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
    if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
    if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
    
    # if data needs additional changes
    ## center data
    data2 <- 
      data %>%
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric(),
                    index_id = interaction(age, period) %>% as.numeric())
    
    ak <- knots$age; pk <- knots$period; ck <- knots$cohort
    afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
    
    formula <- 
      as.formula(y ~  
                   age_id + period_id + cohort_id + 
                   s(age_id, bs = bs, k = ak, fx = afx) + 
                   s(period_id, bs = bs, k = pk, fx = pfx) + 
                   s(cohort_id, bs = bs, k = ck, fx = cfx) +
                   s(index_id, bs = 're'))
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
    fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = data2, method = 'REML', control = ctrl)
    
    return(fit)
  }

## random walk models ----

randomWalk.fit <-
  function(data, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
           randomWalk = c('rw1', 'rw2')[2],
           slopeDrop = NULL, 
           pc.u = 1, pc.alpha = 0.01,
           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
           inla.mode = c('classic', 'twostage', 'experimental')[3],
           control.compute = list(config = TRUE), verbose = FALSE, ...){
    
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2020; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # mod = 'apc'; slopeDrop = 'c'; randomWalk <- c('rw1', 'rw2')[2];
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
                    cohort2_id = cohort_id,
                    index_id = interaction(age, period) %>% as.numeric())
    
    # define the formula for INLA
    # pc hyper priors
    ## temproal
    hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
    
    # full temporal formula
    formula <- 
      y ~ age_id + period_id + cohort_id +
      f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(index_id, model = 'iid', hyper = hyper_pc_time)
    
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
                 data = data2, offset = log(data2$N),
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 safe = TRUE,
                 ...)
    
    return(fit)
  }

# model estimates ----

## spline model ----

spline.estimates <-
  function(dataEst, dataPred,
           mod = c('apc','ac','ap','pc','a','p','c')[1], 
           bs = NULL, knots = NULL, fixed = NULL, slopeDrop = NULL){
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2020; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # mod = 'apc'; slopeDrop = 'c'; bs = 'cr'
    # knots = list(age = 5, period = 5, cohort = 8)
    # fixed = list(age = F,period = F,cohort = F)
    
    
    # the reparameterised fit
    model <- spline.fit(data = dataEst, mod = mod, bs = bs, knots = knots, fixed = fixed, slopeDrop = slopeDrop)
    
    # model fit
    fit <- model
    
    # need to match the naming convention in spline.fit()
    dataPred2 <- 
      dataPred %>% 
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric(),
                    index_id = interaction(age, period) %>% as.numeric())
    
    yhat <- predict(object = fit, newdata = dataPred2, type = 'link')
    
    return(yhat)
    
  }

## random walk model ----

randomWalk.estimates <- 
  function(dataPred, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
           randomWalk = c('rw1', 'rw2')[2],
           slopeDrop = NULL, 
           pc.u = 1, pc.alpha = 0.01,
           control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
           inla.mode = c('classic', 'twostage', 'experimental')[3],
           control.compute = list(config = TRUE), verbose = FALSE, ...){
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2017; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # mod = 'apc'; slopeDrop = 'c'; randomWalk <- c('rw1', 'rw2')[2];
    # pc.u = 1; pc.alpha = 0.01;
    # control.inla = list(strategy = 'adaptive', int.strategy = 'auto');
    # inla.mode = c('classic', 'twostage', 'experimental')[3];
    # control.compute = list(config = TRUE); verbose = FALSE
    
    
    model <- randomWalk.fit(data = dataPred, mod = mod, randomWalk = randomWalk,
                            slopeDrop = slopeDrop, pc.u = pc.u, pc.alpha = pc.alpha,
                            inla.mode = inla.mode, control.compute = control.compute,
                            verbose = verbose)
    
    yhat <- model$summary.linear.predictor$`0.5quant`
    
    return(yhat)
    
  }

# real data ----

## spline model ----

spline.real.fit <- 
  function(dataEst, dataPred, mod = c('apc','ac','ap','pc','a','p','c')[1],
           bs = NULL, knots = NULL, fixed = NULL, slopeDrop = NULL){
    
    # # inputs to run all parts of the function
    # age = 10:99; period = 2000:2020; N = 25;
    # M = 5;
    # data = data.sim(age = age, period = period, N = N, M = M,
    #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
    # mod = 'apc'; slopeDrop = 'c'; bs = 'cr'
    # knots = list(age = 5, period = 5, cohort = 8)
    # fixed = list(age = F,period = F,cohort = F)
    
    if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
    if(is.null(bs)) stop('Warning: Need a basis function.')
    if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
    if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
    if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
    
    # if data needs additional changes
    ## center data
    dataEst2 <- 
      dataEst %>%
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric(),
                    index_id = interaction(age, period) %>% as.numeric())
    
    ak <- knots$age; pk <- knots$period; ck <- knots$cohort
    afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
    
    formula <- 
      as.formula(y ~  
                   age_id + period_id + cohort_id + 
                   s(age_id, bs = bs, k = ak, fx = afx) + 
                   s(period_id, bs = bs, k = pk, fx = pfx) + 
                   s(cohort_id, bs = bs, k = ck, fx = cfx) +
                   s(index_id, bs = 're'))
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
    fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataEst2, method = 'REML', control = ctrl)
    
    # need to match the naming convention in spline.fit()
    dataPred2 <- 
      dataPred %>% 
      dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
                    period_id = period %>% as.factor() %>% as.numeric(),
                    cohort_id = cohort %>% as.factor() %>% as.numeric(),
                    index_id = interaction(age, period) %>% as.numeric())
    
    prediction <- predict(object = fit, newdata = dataPred2, type = 'link', se.fit = TRUE, exclude = 's(index_id)')
    
    results <-
      dataPred %>% 
      dplyr::mutate(yHat = prediction$fit %>% as.vector(),
                    se = prediction$se.fit %>% as.vector()) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(lower = yHat + qnorm(p = 0.025, mean = 0, sd = 1) * se,
                    upper = yHat + qnorm(p = 0.975, mean = 0, sd = 1) * se) %>% 
      dplyr::ungroup()
    
    return(results)
    
  }

## random walk model ----

randomWalk.real.fit <-
  function(data, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
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
                    cohort2_id = cohort_id,
                    index_id = interaction(age, period) %>% as.numeric())
    
    # define the formula for INLA
    # pc hyper priors
    ## temproal
    hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
    
    # full temporal formula
    formula <- 
      y ~ age_id + period_id + cohort_id +
      f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
      f(index_id, model = 'iid', hyper = hyper_pc_time)
    
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
    # fit <- 
    #   INLA::inla(formula, family = 'poisson', 
    #              data = data2, offset = log(data2$N),
    #              control.compute = control.compute,
    #              control.predictor = list(compute = FALSE, link = 1),
    #              control.inla = control.inla,
    #              lincomb = NULL,
    #              inla.mode = inla.mode,
    #              verbose = verbose,
    #              safe = TRUE,
    #              ...)
    fit <- 
      INLA::inla(formula, family = 'poisson', 
                 data = data2, offset = log(data2$N),
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 safe = TRUE)
    
    nSims <- 100
    CI = 0.95; 
    
    # lower and upper CI calculators
    lowerCI <- (1 - CI)/2
    upperCI <- 1 - lowerCI
    
    allSamples <- INLA::inla.posterior.sample(n = nSims, result = fit, intern = TRUE, verbose = FALSE) # <--- remove selection from havard
    
    lineaPredFun <- function(newdata = NA, 
                             mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
                             slopeDrop = NULL){
      
      newdata <- newdata %>% as.data.frame()
      
      if (mod == 'apc'){
        if(slopeDrop == 'a'){
          result <-
            # intercept
            Intercept +
            # linear terms
            period_id * newdata$period_id +
            cohort_id * newdata$cohort_id +
            # random effects
            age2_id[newdata$age2_id] +
            period2_id[newdata$period2_id] +
            cohort2_id[newdata$cohort2_id]
        } else if (slopeDrop == 'p'){
          result <-
            Intercept +
            # linear terms
            age_id * newdata$age_id +
            cohort_id * newdata$cohort_id +
            # random effects
            age2_id[newdata$age2_id] +
            period2_id[newdata$period2_id] +
            cohort2_id[newdata$cohort2_id]
        } else if (slopeDrop == 'c'){
          result <-
            # intercept
            Intercept +
            # linear terms
            age_id * newdata$age_id +
            period_id * newdata$period_id +
            # random effects
            age2_id[newdata$age2_id] +
            period2_id[newdata$period2_id] +
            cohort2_id[newdata$cohort2_id]
        }
      } else if (mod == 'ac'){
        result <-
          # intercept
          Intercept +
          # linear terms
          age_id * newdata$age_id +
          cohort_id * newdata$cohort_id +
          # random effects
          age2_id[newdata$age2_id] +
          cohort2_id[newdata$cohort2_id]
      } else if (mod == 'pc'){
        result <-
          # intercept
          Intercept +
          # linear terms
          period_id * newdata$period_id +
          cohort_id * newdata$cohort_id +
          # random effects
          period2_id[newdata$period2_id] +
          cohort2_id[newdata$cohort2_id]
      } else if (mod == 'ap'){
        result <-
          # intercept
          Intercept +
          # linear terms
          age_id * newdata$age_id +
          period_id * newdata$period_id +
          # random effects
          age2_id[newdata$age2_id] +
          period2_id[newdata$period2_id]
      } else if (mod == 'c'){
        result <-
          # intercept
          Intercept +
          # linear terms
          cohort_id * newdata$cohort_id +
          # random effects
          cohort2_id[newdata$cohort2_id]
      } else if (mod == 'a'){
        result <-
          # intercept
          Intercept +
          # linear terms
          age_id * newdata$age_id +
          # random effects
          age2_id[newdata$age2_id]
      } else if(mod == 'p'){
        result <-
          # intercept
          Intercept +
          # linear terms
          period_id * newdata$period_id +
          # random effects
          period2_id[newdata$period2_id]
      }
      
      return(result)
      
    }
    
    theta <- array(NA, c(nrow(data2), nSims))
    for(i in 1:nrow(data2)){
      theta[i,] <- INLA::inla.posterior.sample.eval(fun = lineaPredFun, samples = allSamples, mod = mod, slopeDrop = slopeDrop, newdata = data2[i,])
    }
    colnames(theta) <- paste0('theta:', 1:nSims)
    theta2 <- cbind(data2, theta)
    
    
    my.summary <- function(x, lowerCI, upperCI) {
      
      # x = theta[,1]; lowerCI = lowerCI; upperCI = upperCI
      
      qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
      data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
      
    }
    
    
    results <-
      theta2 %>% 
      dplyr::mutate(dplyr::select(., starts_with('theta:')) %>%
                      apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>% 
                      lapply(., data.frame) %>%
                      do.call(rbind, .),
                    yHat = median) %>%
      dplyr::select(age, period, cohort, N, y, yHat, lower, upper)
    
    return(results)
  }

# additional functions ----

# summary for estimates
my.summary <- function(x, lowerCI, upperCI) {
  
  # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}

# collect simulation results
collect.results <- function(simulationResults, allData, CI, name){
  
  # setting up
  ## true
  true <- 
    allData[[1]] %>%
    dplyr::select(age, period, cohort) %>%
    dplyr::mutate(true = age.fun(age - mean(age)) + per.fun(period - mean(period)) + coh.fun(cohort - mean(cohort)))
  ## number of sims
  nSims <- length(allData)
  
  # lower and upper CI calculators
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  # colnames
  colnames(simulationResults) <- paste0('theta:', 1:nSims) 
  
  # results
  results <-
    cbind(true , simulationResults) %>% 
    dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                    apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>% 
                    lapply(., data.frame) %>%
                    do.call(rbind, .),
                  bias = 
                    dplyr::across(dplyr::starts_with('theta:'), ~  . - true) %>% 
                    rowMeans(),
                  mse = 
                    dplyr::across(dplyr::starts_with('theta:'), ~  (. - true)^2) %>% 
                    rowMeans(),
                  model = name,
                  type = dplyr::if_else(period > 2017, 'prediction', 'estimate')) %>%
    dplyr::select(-starts_with('theta:'))
  
  results
  
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

# get the plots from simulation runs
get.final.plot <- function(simRun, ...){
  
  # simRun = allResults
  
  # extracting the data frames from the simulation output
  simEst <- simRun$simEstData
  simCurv <- simRun$simCurvData
  truEff <- simRun$trueEffData
  trucurv <- simRun$trueCurvData
  
  # binding the estimate and curvature data frames together
  ## simulation
  simRes <- rbind(simEst %>% mutate(Type = 'est'),
                  simCurv %>% mutate(Type = 'curv')) %>%
    relocate(Type, .before = 's00001')
  
  ## true values
  true <- rbind(truEff %>% mutate(Type = 'est'),
                trucurv %>% mutate(Type = 'curv'))
  
  trueTemp <-
    true %>%
    mutate(truth = Estimate) %>%
    select(Effect, Group, Type, truth)
  
  summary <- 
    simRes %>%
    left_join(., trueTemp, by = c('Effect', 'Group', 'Type')) %>%
    mutate(Estimate = select(., starts_with('s0')) %>% rowMeans,
           Bias = (select(., starts_with('s0')) - truth) %>% rowMeans,
           MSE = ((select(., starts_with('s0')) - truth)^2) %>% rowMeans) %>%
    select(-starts_with('s0'), -truth)
  
  # combining the simulation results and the true results together for specific scenario
  ## simulation results in long format
  simLong <- 
    summary %>% 
    pivot_longer(cols = c('Estimate', 'Bias', 'MSE'),
                 names_to = 'Statistic',
                 values_to = 'Value') %>% 
    filter(!(Type == 'est' & Statistic %in% c('Bias', 'MSE')))
  ## true results in long format
  trueLong <- 
    true %>% 
    pivot_longer(cols = 'Estimate',
                 names_to = 'Statistic',
                 values_to = 'Value')
  # ## combinging long simulation and true data frames
  finalDat <- 
    rbind(simLong, trueLong) %>% 
    mutate(int = factor(interaction(Type, Statistic),
                        levels = c('est.Estimate', 'curv.Estimate', 'curv.Bias', 'curv.MSE'),
                        labels = c('Effect', 'Curvature', 'Curvature Bias', 'Curvature MSE')))
  
  # creating the plots
  ## plots of continuous variables: Effects and curvatures
  p1 <- 
    ggplot()+
    geom_point(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
               aes(x = Group, y = Value, col = Model, shape = Model))+
    geom_line(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
              aes(x = Group, y = Value, col = Model, linetype = Model))+
    facet_grid(int~Effect, scales = 'free')+
    labs(x = 'Relative Years', y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             ...)
  ## plots of discrete variables: curvature bias and MSE box plots
  p2 <- 
    ggplot()+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('Bias')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    geom_jitter(data = finalDat %>% filter(Statistic %in% c('Bias')),
                aes(x = Model, y = Value, col = Model, shape = Model),
                alpha = 0.5)+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('MSE')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    geom_jitter(data = finalDat %>% filter(Statistic %in% c('MSE')),
                aes(x = Model, y = Value, col = Model, shape = Model),
                alpha = 0.5)+
    facet_grid(int~Effect, scales = 'free')+
    labs(y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
             ...)
  ## combining plots vertically
  p <- p1/p2
  ## return combined plot
  # p
  
  # have plots for each model with their simulations over
  ## get the simulation data frame in order
  simFitsLong <- 
    simRes %>% 
    pivot_longer(cols = starts_with('s0'),
                 names_to = 'Simulation',
                 values_to = 'Estimate') %>% 
    mutate(int = factor(Type,
                        levels = c('est', 'curv'),
                        labels = c('Effect', 'Curvature'))) %>% 
    select(-Type)
  # levels(simFitsLong$Type) <- c('Effect', 'Curvature')
  
  # set up loop to creat plot for each model and store as a list
  simPlots <- list()
  models <- unique(simFitsLong$Model)
  for(i in 1:length(models)){
    dataSim <- simFitsLong %>% filter(Model == models[i])
    dataFinal <- finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE')), Model %in% c('True', models[i]))
    p3 <- 
      ggplot()+
      geom_line(data = dataSim, alpha = 0.05,
                aes(x = Group, y = Estimate, col = Model, linetype = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataSim, alpha = 0.05, 
                 aes(x = Group, y = Estimate, col = Model, shape = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataFinal,
                 aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = dataFinal,
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               ...)
    simPlots[[i]] <- p3
  }
  
  individualEffectPlots <- list()
  effects <- c('Age', 'Period', 'Cohort')
  for(i in 1:length(effects)){
    individualEffectDat <- 
      finalDat %>% 
      filter(Effect == effects[i])
    
    # creating the plots
    ## plots of continuous variables: Effects and curvatures
    p1.a <- 
      ggplot()+
      geom_point(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
                 aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               ...)
    ## plots of discrete variables: curvature bias and MSE box plots
    p2.a <- 
      ggplot()+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
                  aes(x = Model, y = Value, col = Model, shape = Model),
                  alpha = 0.5)+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
                  aes(x = Model, y = Value, col = Model, shape = Model),
                  alpha = 0.5)+
      facet_grid(int~Effect, scales = 'free')+
      labs(y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
               ...)
    individualEffectPlots[[i]] <- p1.a/p2.a
  }
  
  return(list(p = p, simPlots = simPlots, individualEffectPlots = individualEffectPlots))
  
}

# get the plots from sensitivity simulation runs
get.sensitivity.final.plot <- function(simRun, ...){
  
  # simRun = allAgeResults
  
  # extracting the data frames from the simulation output
  simEst <- simRun$simEstData
  simCurv <- simRun$simCurvData
  truEff <- simRun$trueEffData
  trucurv <- simRun$trueCurvData
  
  # binding the estimate and curvature data frames together
  ## simulation
  simRes <- rbind(simEst %>% mutate(Type = 'est'),
                  simCurv %>% mutate(Type = 'curv')) %>%
    relocate(Type, .before = 's00001')
  
  ## true values
  true <- rbind(truEff %>% mutate(Type = 'est'),
                trucurv %>% mutate(Type = 'curv'))
  
  trueTemp <-
    true %>%
    mutate(truth = Estimate) %>%
    select(Effect, Group, Type, truth)
  
  summary <- 
    simRes %>%
    left_join(., trueTemp, by = c('Effect', 'Group', 'Type')) %>%
    mutate(Estimate = select(., starts_with('s0')) %>% rowMeans,
           Bias = (select(., starts_with('s0')) - truth) %>% rowMeans,
           MSE = ((select(., starts_with('s0')) - truth)^2) %>% rowMeans) %>%
    select(-starts_with('s0'), -truth)
  
  # combining the simulation results and the true results together for specific scenario
  ## simulation results in long format
  simLong <- 
    summary %>% 
    pivot_longer(cols = c('Estimate', 'Bias', 'MSE'),
                 names_to = 'Statistic',
                 values_to = 'Value') %>% 
    filter(!(Type == 'est' & Statistic %in% c('Bias', 'MSE')))
  ## true results in long format
  trueLong <- 
    true %>% 
    pivot_longer(cols = 'Estimate',
                 names_to = 'Statistic',
                 values_to = 'Value')
  # ## combinging long simulation and true data frames
  finalDat <- 
    rbind(simLong, trueLong) %>% 
    mutate(Model= factor(Model,
                         levels = c('Default', 'PC1', 'PC2', 'PC3', 'Ga1', 'Ga2', 'Ga3', 'True'),
                         labels = c('Default', 'PC1', 'PC2', 'PC3', 'Ga1', 'Ga2', 'Ga3', 'True')),
           int = factor(interaction(Type, Statistic),
                        levels = c('est.Estimate', 'curv.Estimate', 'curv.Bias', 'curv.MSE'),
                        labels = c('Effect', 'Curvature', 'Curvature Bias', 'Curvature MSE')))
  
  # creating the plots
  ## plots of continuous variables: Effects and curvatures
  p1 <- 
    ggplot()+
    # geom_point(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
    #            aes(x = Group, y = Value, col = Model, shape = Model))+
    geom_line(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
              aes(x = Group, y = Value, col = Model, linetype = Model))+
    facet_grid(int~Effect, scales = 'free')+
    labs(x = 'Relative Years', y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             ...)
  ## plots of discrete variables: curvature bias and MSE box plots
  p2 <- 
    ggplot()+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('Bias')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    # geom_jitter(data = finalDat %>% filter(Statistic %in% c('Bias')),
    #             aes(x = Model, y = Value, col = Model, shape = Model),
    #             alpha = 0.5)+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('MSE')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    # geom_jitter(data = finalDat %>% filter(Statistic %in% c('MSE')),
    #             aes(x = Model, y = Value, col = Model, shape = Model),
    #             alpha = 0.5)+
    facet_grid(int~Effect, scales = 'free')+
    labs(y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
             ...)
  ## combining plots vertically
  p <- p1/p2
  ## return combined plot
  # p
  
  # have plots for each model with their simulations over
  ## get the simulation data frame in order
  simFitsLong <- 
    simRes %>% 
    pivot_longer(cols = starts_with('s0'),
                 names_to = 'Simulation',
                 values_to = 'Estimate') %>% 
    mutate(Model = factor(Model,
                          levels = c('Default', 'PC1', 'PC2', 'PC3', 'Ga1', 'Ga2', 'Ga3', 'True'),
                          labels = c('Default', 'PC1', 'PC2', 'PC3', 'Ga1', 'Ga2', 'Ga3', 'True')),
           int = factor(Type,
                        levels = c('est', 'curv'),
                        labels = c('Effect', 'Curvature'))) %>% 
    select(-Type)
  # levels(simFitsLong$Type) <- c('Effect', 'Curvature')
  
  # set up loop to creat plot for each model and store as a list
  simPlots <- list()
  models <- unique(simFitsLong$Model)
  for(i in 1:length(models)){
    dataSim <- simFitsLong %>% filter(Model == models[i])
    dataFinal <- finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE')), Model %in% c('True', models[i]))
    p3 <- 
      ggplot()+
      geom_line(data = dataSim, alpha = 0.05,
                aes(x = Group, y = Estimate, col = Model, linetype = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataSim, alpha = 0.05, 
                 aes(x = Group, y = Estimate, col = Model, shape = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataFinal,
                 aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = dataFinal,
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               ...)
    simPlots[[i]] <- p3
  }
  
  individualEffectPlots <- list()
  effects <- c('Age', 'Period', 'Cohort')
  for(i in 1:length(effects)){
    individualEffectDat <- 
      finalDat %>% 
      filter(Effect == effects[i])
    
    # creating the plots
    ## plots of continuous variables: Effects and curvatures
    p1.a <- 
      ggplot()+
      # geom_point(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
      #            aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               ...)
    ## plots of discrete variables: curvature bias and MSE box plots
    p2.a <- 
      ggplot()+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      # geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
      #             aes(x = Model, y = Value, col = Model, shape = Model),
      #             alpha = 0.5)+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      # geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
      #             aes(x = Model, y = Value, col = Model, shape = Model),
      #             alpha = 0.5)+
      facet_grid(int~Effect, scales = 'free')+
      labs(y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
               ...)
    individualEffectPlots[[i]] <- p1.a/p2.a
  }
  
  return(list(p = p, simPlots = simPlots, individualEffectPlots = individualEffectPlots))
  
}


# # OLD:::: real data ----
# 
# ## spline model ----
# 
# spline.real.fit <- 
#   function(dataEst, dataPred, mod = c('apc','ac','ap','pc','a','p','c')[1],
#            bs = NULL, knots = NULL, fixed = NULL, slopeDrop = NULL){
#     
#     # # inputs to run all parts of the function
#     # age = 10:99; period = 2000:2020; N = 25;
#     # M = 5;
#     # data = data.sim(age = age, period = period, N = N, M = M,
#     #                 FUNage = age.fun,  FUNperiod = per.fun, FUNcohort = coh.fun)
#     # mod = 'apc'; slopeDrop = 'c'; bs = 'cr'
#     # knots = list(age = 5, period = 5, cohort = 8)
#     # fixed = list(age = F,period = F,cohort = F)
#     
#     if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
#     if(is.null(bs)) stop('Warning: Need a basis function.')
#     if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
#     if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
#     if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
#     
#     # if data needs additional changes
#     ## center data
#     dataEst2 <- 
#       dataEst %>%
#       dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
#                     period_id = period %>% as.factor() %>% as.numeric(),
#                     cohort_id = cohort %>% as.factor() %>% as.numeric(),
#                     index_id = interaction(age, period) %>% as.numeric())
#     
#     ak <- knots$age; pk <- knots$period; ck <- knots$cohort
#     afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
#     
#     formula <- 
#       as.formula(y ~  
#                    age_id + period_id + cohort_id + 
#                    s(age_id, bs = bs, k = ak, fx = afx) + 
#                    s(period_id, bs = bs, k = pk, fx = pfx) + 
#                    s(cohort_id, bs = bs, k = ck, fx = cfx) +
#                    s(index_id, bs = 're'))
#     # update pre and post forumla for the temporal models
#     # # do not remove any linear terms from pre fit as they are all needed later in 
#     # # re-parameterisation
#     if (mod == 'apc'){
#       if(slopeDrop == 'a'){
#         formula <- update(formula, ~. - age_id)
#       } else if (slopeDrop == 'p'){
#         formula <- update(formula, ~. - period_id)
#       } else if (slopeDrop == 'c'){
#         formula <- update(formula, ~. - cohort_id)
#       }
#     } else if (mod == 'ac'){
#       formula <- update(formula, ~.- period_id - s(period_id, bs = bs, k = pk, fx = pfx))
#     } else if (mod == 'pc'){
#       formula <- update(formula, ~. - age_id - s(age_id, bs = bs, k = ak, fx = afx))
#     } else if (mod == 'ap'){
#       formula <- update(formula, ~. - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
#     } else if (mod == 'c'){
#       formula <- update(formula, ~.- age_id - s(age_id, bs = bs, k = ak, fx = afx) - period_id - s(period_id, bs = bs, k = pk, fx = pfx))
#     } else if (mod == 'a'){
#       formula <- update(formula, ~.- period_id - s(period_id, bs = bs, k = pk, fx = pfx) - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
#     } else if (mod == 'p'){
#       formula <- update(formula, ~.- age_id - s(age_id, bs = bs, k = ak, fx = afx) - cohort_id - s(cohort_id, bs = bs, k = ck, fx = cfx))
#     }
#     
#     # use all cores avaliable
#     ctrl <- list(nthreads = parallel::detectCores())
#     
#     # model fit
#     fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataEst2, method = 'REML', control = ctrl)
#     
#     # need to match the naming convention in spline.fit()
#     dataPred2 <- 
#       dataPred %>% 
#       dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
#                     period_id = period %>% as.factor() %>% as.numeric(),
#                     cohort_id = cohort %>% as.factor() %>% as.numeric(),
#                     index_id = interaction(age, period) %>% as.numeric())
#     
#     prediction <- predict(object = fit, newdata = dataPred2, type = 'link', se.fit = TRUE)
#     
#     results <-
#       dataPred %>% 
#       dplyr::mutate(yHat = prediction$fit %>% as.vector(),
#                     se = prediction$se.fit %>% as.vector()) %>% 
#       dplyr::rowwise() %>% 
#       dplyr::mutate(lower = yHat + qnorm(p = 0.025, mean = 0, sd = 1) * se,
#                     upper = yHat + qnorm(p = 0.975, mean = 0, sd = 1) * se) %>% 
#       dplyr::ungroup()
#     
#     return(results)
#     
#   }
# 
# ## random walk model ----
# 
# randomWalk.real.fit <-
#   function(data, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1],
#            randomWalk = c('rw1', 'rw2')[2],
#            slopeDrop = NULL, 
#            pc.u = 1, pc.alpha = 0.01,
#            control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
#            inla.mode = c('classic', 'twostage', 'experimental')[3],
#            control.compute = list(config = TRUE), verbose = FALSE, ...){
#     
#     
#     # # inputs to run all parts of the function
#     # data = alcoholData.randomWalkPredict
#     # mod = 'apc'; slopeDrop = 'c'; randomWalk = 'rw1';
#     # pc.u = 1; pc.alpha = 0.01;
#     # control.inla = list(strategy = 'adaptive', int.strategy = 'auto');
#     # inla.mode = c('classic', 'twostage', 'experimental')[3];
#     # control.compute = list(config = TRUE); verbose = FALSE
#     
#     if (!"config" %in% names(control.compute)) {
#       message("config = TRUE is added to control.compute so that posterior draws can be taken.")
#       control.compute$config <- TRUE
#     }
#     if(mod == 'apc' && is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
#     if(mod == 'apc' &&! (slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
#     
#     # data augmentation
#     data2 <- 
#       data %>%
#       dplyr::mutate(age_id = age %>% as.factor() %>% as.numeric(),
#                     period_id = period %>% as.factor() %>% as.numeric(),
#                     cohort_id = cohort %>% as.factor() %>% as.numeric(),
#                     age2_id = age_id,
#                     period2_id = period_id,
#                     cohort2_id = cohort_id,
#                     index_id = interaction(age, period) %>% as.numeric())
#     
#     # define the formula for INLA
#     # pc hyper priors
#     ## temproal
#     hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
#     
#     # full temporal formula
#     formula <- 
#       y ~ age_id + period_id + cohort_id +
#       f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
#       f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
#       f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) +
#       f(index_id, model = 'iid', hyper = hyper_pc_time)
#     
#     # update pre and post forumla for the temporal models
#     # # do not remove any linear terms from pre fit as they are all needed later in 
#     # # re-parameterisation
#     if (mod == 'apc'){
#       if(slopeDrop == 'a'){
#         formula <- update(formula, ~. - age_id)
#       } else if (slopeDrop == 'p'){
#         formula <- update(formula, ~. - period_id)
#       } else if (slopeDrop == 'c'){
#         formula <- update(formula, ~. - cohort_id)
#       }
#     } else if (mod == 'ac'){
#       formula <- update(formula, ~.
#                         - period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     } else if (mod == 'pc'){
#       formula <- update(formula, ~.
#                         - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     } else if (mod == 'ap'){
#       formula <- update(formula, ~.
#                         - cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     } else if (mod == 'c'){
#       formula <- update(formula, ~.
#                         - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
#                           period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     } else if (mod == 'a'){
#       formula <- update(formula, ~.
#                         - period_id - f(period2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
#                           cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     } else if(mod == 'p'){
#       formula <- update(formula, ~.
#                         - age_id - f(age2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE) -
#                           cohort_id - f(cohort2_id, model = randomWalk, hyper = hyper_pc_time, constr = TRUE))
#     }
#     
#     # model fit
#     fit <- 
#       INLA::inla(formula, family = 'poisson', 
#                  data = data2, offset = log(data2$N),
#                  control.compute = control.compute,
#                  control.predictor = list(compute = FALSE, link = 1),
#                  control.inla = control.inla,
#                  lincomb = NULL,
#                  inla.mode = inla.mode,
#                  verbose = verbose,
#                  safe = TRUE,
#                  ...)
#     
#     results <-
#       data %>% 
#       dplyr::mutate(yHat = fit$summary.linear.predictor$`0.5quant`,
#                     lower = fit$summary.linear.predictor$`0.025quant`,
#                     upper = fit$summary.linear.predictor$`0.975quant`)
#     
#     return(results)
#   }
