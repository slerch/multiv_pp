# (univariate) post-processing of ensemble forecast object "ensfc"

# apply post-processing separately in each dimension, and return the parameters of the normal distribution in out of sample evaluation period

# requirements:
#   installation of the 'scoringRules' package with version >= 0.9.2

# input:
#   fcmodel: character string; indicating which model wass used to generate ensfc
#   ensfc, ensfc_init: ensemble forecasts to be used in training and computation of parameters, 
#       should be output of generate_ensfc.R
#   obs, obs_init: observations to be used in training,
#       should be output of generate_observations.R
#   train: string indicating type of training: 
#       "init" = use initial training period,
#       "moving" = use last trainlength instances
#   trainlength: number of forecast instances to be used:
#       if train = "init", use first trainlength instances
#       if train = "moving", this parameter is required
#       if no value is assigned, trainlength is set to the maximum possible number
#   emos_plus: logical indicating whether variance coefficients should be restricted to positive values
#       use this if postprocessing produces NA or NaN parameter values 
#       (may happen for example with short training periods or extreme parameter choices)
#   lower/upper: only required for Setting 2 (vector of lower / upper bounds)

# output:
#   array of forecast distribution parameters, 
#       first dimension: forecast instance in evaulation period
#       second dimension: dimension in forecasting problem
#       third dimension: vector of parameter values


postproc <- function(fcmodel, ensfc, ensfc_init, obs, obs_init, train, trainlength = NULL, emos_plus = FALSE, lower = NULL, upper = NULL){
  require(scoringRules)
  
  # Setting 1
  if(fcmodel == 1){
    
    # extract internal CRPS function for normal distributions
    crps_norm <- scoringRules::crps_norm
    
    # objective function for minimum CRPS estimation of EMOS coefficients
    objective_fun <- function(par, ens_mean, ens_var, obs_train){
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      return(sum(crps_norm(y = obs_train, location = m, scale = s)))
    }
    
    #
    #
    # newer scoringRules versions now contain gradient function for CRPS of normal distribution, might be used instead!
    #
    #
    
    # CRPS gradient
    crpsgrad <- function(y, location, scale){
      z <- y
      if(!identical(location, 0) | !identical(scale, 1)){
        z <- (y - location) / scale
      }
      dmu <- 1-2*pnorm(z)
      dsig <- 2*dnorm(z)-1/sqrt(pi)
      return(cbind(dmu, dsig))
    }
    
    # wrapper for gradient function to use in optim
    gradfun_wrapper <- function(par, obs_train, ens_mean, ens_var){
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      dcrps_dtheta <- crpsgrad(y = obs_train, location = m, scale = s)
      out1 <- dcrps_dtheta[,1] %*% cbind(1, ens_mean)
      out2 <- dcrps_dtheta[,2] %*% cbind(1/(2*sqrt(par[3]+par[4]*ens_var)), 
                                         ens_var/(2*sqrt(par[3]+par[4]*ens_var)))
      return(cbind(out1, out2))
    }
    
    # check dimensions of input ensfc and obs objects
    if(any(dim(obs_init) != dim(ensfc_init)[c(1,3)])){
      stop("dimensions of ensfc_init and obs_init do not match")
    }
    # if(train == "moving" & any(dim(obs_init) != dim(ensfc_init)[c(1,3)])){
    #   stop("dimensions of ensfc and obs do not match")
    # }
    if(any(dim(ensfc_init)[c(2,3)] != dim(ensfc)[c(2,3)])){
      stop("dimensions of ensfc_init and ensfc do not match")
    }
    
    # check training parameters
    if(!is.element(train, c("init", "moving"))){
      stop("'train' needs to be 'init' or 'moving'.")
    }
    if(train == "moving"){
      stop("currently only implemented for train = 'init'")
    }
    if(train == "init" & is.null(trainlength)){
      trainlength <- dim(ensfc_init)[1]
    }
    if(trainlength > dim(ensfc_init)[1]){
      trainlength <- dim(ensfc_init)[1]
      message("specified 'trainlength' is too large, using ninit instead")
    }
    
    
    # generate arrays to save parameter values to
    # dimensions: forecast instance in evaluation period; dimension; parameters 
    n <- dim(ensfc)[1]
    d <- dim(ensfc)[3]
    emos_param <- array(NA, dim = c(n, d, 2))
    
    # iterating over dimensions, estimate EMOS coefficients using minimum CRPS estimation
    emos_coefs <- matrix(NA, nrow = d, ncol = 4)
    for(dd in 1:d){
      # extract cases in training period
      ensfc_tr <- ensfc_init[,,dd]
      ensfc_tr_mean <- apply(ensfc_tr, 1, mean)
      ensfc_tr_var <- apply(ensfc_tr, 1, var)
      obs_tr <- obs_init[,dd]
      # estimate EMOS coefficients
      # try BFGS
      coefs <- NULL
      if(!emos_plus){
        try(coefs <- optim(par = c(1,1,1,1),
                           fn = objective_fun,
                           ens_mean = ensfc_tr_mean,
                           ens_var = ensfc_tr_var,
                           obs_train = obs_tr,
                           method = "BFGS",
                           gr = gradfun_wrapper)$par,
            silent = TRUE)
      }
      # if BFGS didn't work or emos_plus specified as TRUE, use L-BFGS-B instead
      if(emos_plus | class(coefs) == "try-error"){
        coefs <- optim(par = c(1,1,1,1),
                       fn = objective_fun,
                       ens_mean = ensfc_tr_mean,
                       ens_var = ensfc_tr_var,
                       obs_train = obs_tr,
                       method = "L-BFGS-B",
                       lower = c(-20, -20, 1e-5, 1e-5),
                       upper = rep(20, 4),
                       gr = gradfun_wrapper)$par
      }
      # save coefficients
      emos_coefs[dd,] <- coefs
    }
    
    # interating over days in the evaluation period and dimensions:
    #   compute EMOS parameters = parameters of (normal) forecast distributions 
    for(nn in 1:n){
      for(dd in 1:d){
        ensfc_nndd <- ensfc[nn,,dd]
        ensfc_nndd_mean <- mean(ensfc_nndd)
        ensfc_nndd_var <- var(ensfc_nndd)
        loc <- c(1, ensfc_nndd_mean) %*% emos_coefs[dd, 1:2]
        sc <- sqrt(c(1, ensfc_nndd_var) %*% emos_coefs[dd, 3:4])
        emos_param[nn,dd,] <- c(loc, sc)
      }
    }
    
    return(emos_param)
  }
  
  # Setting 2
  if(fcmodel == 2){
    
    # ensfc <- ens$ensfc
    # ensfc_init <- ens$ensfc_init
    # 
    # obs <- verobs$obs
    # obs_init <- verobs$obs_init
    
    
    # check whether the same bounds are used for input and output
    # if(!identical(ens$lower,verobs$lower)){
    #   stop("lower bounds for forecasts and observations do not match")
    # }
    # 
    # if(!identical(ens$upper,verobs$upper)){
    #   stop("upper bounds for forecasts and observations do not match")
    # }
    # 
    
    lowerB <- lower
    upperB <- upper
    
    # extract internal CRPS function for normal distributions
    crps_tnorm <- scoringRules:::crps_tnorm
    
    # objective function for minimum CRPS estimation of EMOS coefficients
    objective_fun <- function(par, ens_mean, ens_var, obs_train, lower_bound, upper_bound){
      par[3:4] <- par[3:4]^2 #force positive variance
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      return(sum(crps_tnorm(y = obs_train, location = m, scale = s, lower = lower_bound, upper = upper_bound)))
    }
    
    
    # CRPS gradient
    crpsgrad <- scoringRules:::gradcrps_tnorm
    
    
    # wrapper for gradient function to use in optim
    gradfun_wrapper <- function(par, obs_train, ens_mean, ens_var, lower_bound, upper_bound){
      par[3:4] <- par[3:4]^2
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      dcrps_dtheta <- crpsgrad(y = obs_train, location = m, scale = s, lower = lower_bound, upper = upper_bound)
      out1 <- dcrps_dtheta[,1] %*% cbind(1, ens_mean)
      out2 <- dcrps_dtheta[,2] %*% cbind(1/(2*sqrt(par[3]+par[4]*ens_var)), 
                                         ens_var/(2*sqrt(par[3]+par[4]*ens_var)))
      return(cbind(out1, out2))
    }
    
    # check dimensions of input ensfc and obs objects
    if(any(dim(obs_init) != dim(ensfc_init)[c(1,3)])){
      stop("dimensions of ensfc_init and obs_init do not match")
    }
    # if(train == "moving" & any(dim(obs_init) != dim(ensfc_init)[c(1,3)])){
    #   stop("dimensions of ensfc and obs do not match")
    # }
    if(any(dim(ensfc_init)[c(2,3)] != dim(ensfc)[c(2,3)])){
      stop("dimensions of ensfc_init and ensfc do not match")
    }
    
    # check training parameters
    if(!is.element(train, c("init", "moving"))){
      stop("'train' needs to be 'init' or 'moving'.")
    }
    if(train == "moving"){
      stop("currently only implemented for train = 'init'")
    }
    if(train == "init" & is.null(trainlength)){
      trainlength <- dim(ensfc_init)[1]
    }
    if(trainlength > dim(ensfc_init)[1]){
      trainlength <- dim(ensfc_init)[1]
      message("specified 'trainlength' is too large, using ninit instead")
    }
    
    
    # generate arrays to save parameter values to
    # dimensions: forecast instance in evaluation period; dimension; parameters 
    n <- dim(ensfc)[1]
    d <- dim(ensfc)[3]
    
    
    emos_param <- array(NA, dim = c(n, d, 2))
    
    # iterating over dimensions, estimate EMOS coefficients using minimum CRPS estimation
    emos_coefs <- matrix(NA, nrow = d, ncol = 4)
    for(dd in 1:d){
      # extract cases in training period
      ensfc_tr <- ensfc_init[,,dd]
      ensfc_tr_mean <- apply(ensfc_tr, 1, mean)
      ensfc_tr_var <- apply(ensfc_tr, 1, var)
      obs_tr <- obs_init[,dd]
      olsCoefs <-  as.vector(lm(as.vector(obs_tr) ~ (ensfc_tr_mean))$coef)
      # estimate EMOS coefficients
      # try BFGS
      coefs <- NULL
      if(!emos_plus){
        try(coefs <- optim(par = c(olsCoefs,1,1),
                           fn = objective_fun,
                           ens_mean = ensfc_tr_mean,
                           ens_var = ensfc_tr_var,
                           obs_train = obs_tr,
                           lower_bound = lowerB[dd],
                           upper_bound = upperB[dd],
                           method = "BFGS",
                           gr = gradfun_wrapper)$par,
            silent = F)
      }
      # if BFGS didn't work or emos_plus specified as TRUE, use L-BFGS-B instead
      if(emos_plus | class(coefs) == "try-error"){
        coefs <- optim(par = c(olsCoefs,1,1),
                       fn = objective_fun,
                       ens_mean = ensfc_tr_mean,
                       ens_var = ensfc_tr_var,
                       obs_train = obs_tr,
                       lower_bound = lowerB[dd],
                       upper_bound = upperB[dd],
                       method = "L-BFGS-B",
                       lower = c(-20, -20, 1e-5, 1e-5),
                       upper = rep(20, 4),
                       gr = gradfun_wrapper)$par
      }
      # save coefficients
      coefs[3:4] <- coefs[3:4]^2
      emos_coefs[dd,] <- coefs
    }
    
    # iterating over days in the evaluation period and dimensions:
    # compute EMOS parameters = parameters of (truncated normal) forecast distributions 
    for(nn in 1:n){
      for(dd in 1:d){
        ensfc_nndd <- ensfc[nn,,dd]
        ensfc_nndd_mean <- mean(ensfc_nndd)
        ensfc_nndd_var <- var(ensfc_nndd)
        loc <- c(1, ensfc_nndd_mean) %*% emos_coefs[dd, 1:2]
        sc <- sqrt(c(1, ensfc_nndd_var) %*% emos_coefs[dd, 3:4])
        emos_param[nn,dd,] <- c(loc, sc)
      }
    }
    
    return(list("parameters" = emos_param, "model" = "tnormal", "lower" = lowerB, "upper" = upperB))
  }

  # Setting 4
  if(fcmodel == 4){
    
    if(train != "moving"){
      stop("Should choose train = 'moving' for this setting") 
    }
    
    # extract internal CRPS function for normal distributions
    crps_norm <- scoringRules::crps_norm
    
    # objective function for minimum CRPS estimation of EMOS coefficients
    objective_fun <- function(par, ens_mean, ens_var, obs_train){
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      return(sum(crps_norm(y = obs_train, location = m, scale = s)))
    }
    
    # CRPS gradient
    crpsgrad <- function(y, location, scale){
      z <- y
      if(!identical(location, 0) | !identical(scale, 1)){
        z <- (y - location) / scale
      }
      dmu <- 1-2*pnorm(z)
      dsig <- 2*dnorm(z)-1/sqrt(pi)
      return(cbind(dmu, dsig))
    }
    
    # wrapper for gradient function to use in optim
    gradfun_wrapper <- function(par, obs_train, ens_mean, ens_var){
      m <- cbind(1, ens_mean) %*% par[1:2]
      s <- sqrt(cbind(1, ens_var) %*% par[3:4])
      dcrps_dtheta <- crpsgrad(y = obs_train, location = m, scale = s)
      out1 <- dcrps_dtheta[,1] %*% cbind(1, ens_mean)
      out2 <- dcrps_dtheta[,2] %*% cbind(1/(2*sqrt(par[3]+par[4]*ens_var)), 
                                         ens_var/(2*sqrt(par[3]+par[4]*ens_var)))
      return(cbind(out1, out2))
    }
    
    # generate arrays to save parameter values to
    # dimensions: forecast instance in evaluation period; dimension; parameters 
    n <- dim(ensfc)[1]
    d <- dim(ensfc)[3]
    emos_param <- array(NA, dim = c(n, d, 2))
    
    # iterating over training cases fron 1:n, estimate EMOS coefficients using minimum CRPS estimation iterated over dimensions
    
    emos_coefs <- matrix(NA, nrow = d, ncol = 4)
    train_length <- trainlength
    
    ensfc_all <- abind::abind(ensfc_init, ensfc, along = 1)
    obs_all <- rbind(obs_init, obs)
    
    ninit <- dim(ensfc_init)[1]
    
    for(nn in 1:n){
      
      for(dd in 1:d){
        
        ensfc_tr <- ensfc_all[(ninit + nn - train_length + 1):(ninit + nn),,dd]
        
        ensfc_tr_mean <- apply(ensfc_tr, 1, mean)
        ensfc_tr_var <- apply(ensfc_tr, 1, var)
        
        obs_tr <- obs_all[(ninit + nn - train_length + 1):(ninit + nn),dd]
        
        coefs <- NULL
        if(!emos_plus){
          try(coefs <- optim(par = c(1,1,1,1),
                             fn = objective_fun,
                             ens_mean = ensfc_tr_mean,
                             ens_var = ensfc_tr_var,
                             obs_train = obs_tr,
                             method = "BFGS",
                             gr = gradfun_wrapper)$par,
              silent = TRUE)
        }
        # if BFGS didn't work or emos_plus specified as TRUE, use L-BFGS-B instead
        if(emos_plus | class(coefs) == "try-error"){
          coefs <- optim(par = c(1,1,1,1),
                         fn = objective_fun,
                         ens_mean = ensfc_tr_mean,
                         ens_var = ensfc_tr_var,
                         obs_train = obs_tr,
                         method = "L-BFGS-B",
                         lower = c(-20, -20, 1e-5, 1e-5),
                         upper = rep(20, 4),
                         gr = gradfun_wrapper)$par
        }
        
        emos_coefs[dd,] <- coefs
        
        ensfc_nndd <- ensfc[nn,,dd]
        ensfc_nndd_mean <- mean(ensfc_nndd)
        ensfc_nndd_var <- var(ensfc_nndd)
        loc <- c(1, ensfc_nndd_mean) %*% emos_coefs[dd, 1:2]
        sc <- sqrt(c(1, ensfc_nndd_var) %*% emos_coefs[dd, 3:4])
        emos_param[nn,dd,] <- c(loc, sc)
      }
    }
    
    return(emos_param)
  }
  
}