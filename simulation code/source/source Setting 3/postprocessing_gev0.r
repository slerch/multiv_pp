# -------------------------------------------------------------------------- #

# (univariate) post-processing of ensemble forecast object "ensfc"

# apply post-processing separately in each dimension, and return the parameters
# of the estimated distribution in out of sample evaluation period

# requirements:
#   installation of the 'scoringRules' package with version >= 0.9.2

# input:
#   fcmodel: character string; indicating which model wass used to generate ensfc
#            within this implementation "gev0" needs to be specified
#   ensfc, ensfc_init: ensemble forecasts to be used in training and computation of parameters,
#       should be output of sim_gev0_ens.R
#   obs, obs_init: observations to be used in training,
#       should be output of sim_gev0_obs.R
#   train: string indicating type of training:
#       "init" = use initial training period,
#       "moving" = use last trainlength instances
#   trainlength: number of forecast instances to be used:
#       if train = "init", use first trainlength instances
#       if train = "moving", this parameter is required, but for gev0 this
#                            is not yet implemented
#       if no value is assigned, trainlength is set to the maximum possible number
#   emos_plus: logical indicating whether variance coefficients should be restricted to positive values
#       use this if postprocessing produces NA or NaN parameter values
#       (may happen for example with short training periods or extreme parameter choices)

# output:
#   array of forecast distribution parameters,
#       first dimension: forecast instance in evaulation period
#       second dimension: dimension in forecasting problem
#       third dimension: vector of parameter values




postproc <- function(ens, verobs, train, trainlength = NULL, emos_plus = FALSE)
{
  require(scoringRules)
  require(evd)

  # check whether the same model is used for input and output
  if(ens$model != verobs$model){
    stop("models for forecasts and observations do not match")
  }


  ensfc <- ens$ensfc
  ensfc_init <- ens$ensfc_init

  obs <- verobs$obs
  obs_init <- verobs$obs_init


# Function to compute ensemble mean difference, used as predictor in EMOS for precipitation,
# See Scheuerer 2014
    gini.md <- function(x,na.rm=FALSE)
{
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)]
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n*(n-1))-2*mean(x)*(n+1)/(n-1))
}


# CRPS of GEV censored at 0, for shape parameter unequal to 0, see Scheuerer 2014
# and ensembleMOS implementation
crps.GEVneq0 <- function(param, obs_train, ensMean, ensProb0, ensMeanDiff)
{
   MEAN <- param[1] + param[2]*ensMean + param[3]*ensProb0
   SCALE <- param[4] + param[5]^2*ensMeanDiff
   SHAPE <- param[6]
   LOC <- MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE

   SCdSH <- SCALE/SHAPE
   Gam1mSH <- gamma(1-SHAPE)
   prob0 <- pgev(0, loc=LOC, scale=SCALE, shape=SHAPE)
   probY <- pgev(obs_train, loc=LOC, scale=SCALE, shape=SHAPE)

   T1 <- (obs_train-LOC)*(2*probY-1) + LOC*prob0^2
   T2 <- SCdSH * ( 1-prob0^2 - 2^SHAPE*Gam1mSH*pgamma(-2*log(prob0),1-SHAPE) )
   T3 <- -2*SCdSH * ( 1-probY - Gam1mSH*pgamma(-log(probY),1-SHAPE) )
   return( mean(T1+T2+T3) )
}


# In Scheuerer 2014 and also in the corresponding ensembleMOS implementation
# a closed form CRPS expression is available in case shape parameter unequal to 0
# If it is 0, a weighted average of above CRPS for shape unequal to 0 but within
# some small epsilon interval is computed to approximate the shape equal 0 case.
# For this, eps is set so small value, see below objective function

eps <- 1e-05


 # objective function for minimum CRPS estimation of EMOS coefficients
 
objective_fun <- function(param, obs_train, ensMean, ensProb0, ensMeanDiff)
{
        if (emos_plus) {
            param[4] <- param[4]^2
        }
        if (abs(param[6]) < eps) {
            crps.eps.m <- crps.GEVneq0(c(param[-6], -eps), obs_train,
                ensMean, ensProb0, ensMeanDiff)
            crps.eps.p <- crps.GEVneq0(c(param[-6], eps), obs_train,
                ensMean, ensProb0, ensMeanDiff)
            w.m <- (eps - param[6])/(2 * eps)
            w.p <- (eps + param[6])/(2 * eps)
            res <- w.m * crps.eps.m + w.p * crps.eps.p
        }
        else {
            res <- crps.GEVneq0(param, obs_train, ensMean, ensProb0, ensMeanDiff)
        }
        return(res)
}



    # check dimensions of input ensfc and obs objects
      if(any(dim(obs_init) != dim(ensfc_init)[c(1,3)])){
        stop("dimensions of ensfc_init and obs_init do not match")
      }

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



    # generate array to save EMOS parameter values to
    # dimensions: forecast instance in evaluation period; dimension; parameters
    # Separately define matrix to save location parameter
    # (EMOS estimates mean parameter of GEV, scale and shape parameter,
    # but not directly location parameter, this will derived
    # from the other estimates)
      n <- dim(ensfc)[1]
      d <- dim(ensfc)[3]
      emos_param <- array(NA, dim = c(n, d, 2))
      mu_param <- matrix(NA, nrow=n, ncol=d)

    # iterating over dimensions, estimate EMOS coefficients using minimum CRPS estimation
    # In GEV0 EMOS model, the paramters alpha0, alpha1, alpha2, beta0, beta1, and shape
    # are estimated by optim, so emos_coefs has 6 columns
      emos_coefs <- matrix(NA, nrow = d, ncol = 6)
      
      for(dd in 1:d){
      # extract cases in training period

      # For the Scheuerer/GEV0 EMOS model, ensemble mean, ensemble mean difference and fraction of members predicting
      # zero precipitation are required, thus computed below
        ensfc_tr <- ensfc_init[,,dd]
        ensfc_tr_mean <- apply(ensfc_tr, 1, mean)
        ensfc_Prob0 <- apply(ensfc_tr==0, 1, mean)
        ensfc_meanDiff  <- apply(ensfc_tr, 1, gini.md)
        obs_tr <- obs_init[,dd]

      # In ensembleMOS starting values for the coefficients alpha0 and alpha1 are estimated
      # by a linear model, and only in case "enough" non-zero observations are available
      # In the ensembleMOS implementation "enough" means number of obs > twice the number of parameters
      # in our simple case the number of parameters is 2, so 2*2 non-zero observations are necessary

        pozobs <- (obs_tr > 0)
        if (sum(pozobs) > 2 * (1 + 1)) {
            olsCoefs <- as.vector(lm(as.vector(obs_tr[pozobs]) ~ (ensfc_tr_mean[pozobs]))$coef)
        } else {olsCoefs <- c(1,1)}
        # if not enough non-zero obs set starting values to standard values (taken from ensembleMOS implementation)


      # estimate EMOS coefficients
      # (Unconstrained) BFGS for GEV left-censored at 0 currently not used,
      # as it seems to yield unreliable parameter estimates or even error/warning
      # messages due to non-admissable parameter estimates
      # in some iteration steps

        coefs <- NULL

      # So for now only use L-BFGS-B
        if(emos_plus) {b0 <- sqrt(0.1)} else {b0 <- 0.1}

          coefs <- optim(par = c(olsCoefs,1,b0, sqrt(1),0.1),  # starting values taken from ensembleMOS implementation
                         fn = objective_fun,
                         ensMean = ensfc_tr_mean,
                         ensProb0 = ensfc_Prob0,
                         ensMeanDiff = ensfc_meanDiff,
                         obs_train = obs_tr,
                         method = "L-BFGS-B",
                         # box constraints taken from ensembleMOS implementation
                         lower = c(1e-04, 1e-04, -10, 1e-04, 1e-04, -0.278),
                         upper = c(10, 10, 10, 10, 10, 0.99999),
                         gr = NULL)$par



      # save coefficients
        coefs[5] <- coefs[5]^2
        if(emos_plus){coefs[4] <- coefs[4]^2}
        emos_coefs[dd,] <- coefs
      }


    # iterating over days in the evaluation period and dimensions:
    # compute EMOS parameters = parameters of GEV0 forecast distributions
      for(nn in 1:n){
        for(dd in 1:d){
          ensfc_nndd <- ensfc[nn,,dd]
          ensfc_nndd_mean <- mean(ensfc_nndd)
          ensfc_nndd_prob0 <- mean(ensfc_nndd==0)
          ensfc_nndd_md <- gini.md(ensfc_nndd)
          shape <- emos_coefs[dd,6]
          # Attention: loc is expected value (not location) of non-censored GEV
          # sc is scale parameter of non-censored GEV
          loc <- c(1, ensfc_nndd_mean, ensfc_nndd_prob0) %*% emos_coefs[dd, 1:3]
          sc <- sqrt(c(1, ensfc_nndd_md) %*% emos_coefs[dd, 4:5])
          # mu is location parameter of non-censored GEV, needed later on to evaulate e.g. estimated CDF
          mu <- if(abs(shape) > 10^(-6)) {loc - sc*(gamma(1-shape)-1)/shape} else {loc - sc * 0.57721566490153}
          emos_param[nn,dd,] <- c(loc, sc)
          mu_param[nn,dd] <- mu
        }
      }

    # Additonally to predictive mean and scale of EMOS, save the estimated location and shape paramter
    # for use in multivariate postprocessing calls
     return(list("parameters" = emos_param, "model" = "gev0", "location" = mu_param, "shape" = emos_coefs[,6]))



  }