# -------------------------------------------------------------------------- #

# multivariate post-processing methods, applied to result of univariate post-processing

# the following methods are applied:
#   - EMOS-X: not accounting for any multivariate dependencies,
#             independet samples from forcast distributions in margins
#       EMOS-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       EMOS-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       EMOS-R: random sample
#       EMOS-S: stratifed sampling approach of Hu et al (2016)
#       EMOS-T: fit parametric forecast distribution (eg normal, or gev0) to ensemble,
#               and use quantiles corresponding to ens forecasts as quantile levels
#               for draws from post-processed distribution
#               --> For gev0, currently set in comments, as for gev0 ensemble forecasts
#                   all ensemble members could take the value 0, in that case neither ensemble variance
#                   nor ensemble mean difference can be computed in a reasonable way
#                   so other alternative ensemble statitistics need to be considered in future research
#   - ECC-X: inheriting the multivariate dependency structure from the ensemble forecast
#             draw samples (as in EMOS-X), and re-order
#       ECC-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       ECC-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       ECC-R: random sample
#       ECC-S: stratifed sampling approach of Hu et al (2016)
#       ECC-T: see above, for gev0 set in comments due to issues with all ensemble forecasts possibly 0 sometimes
#   - dECC: dynamic ECC, ECC variant proposed by Bouallegue et al. (2016, MWR)
#       ...
#       there are also different variants here!
#   - GCA: Gaussian copula approach from Pinson and Girard (2012), Möller et al (2013)
#       GCA: with correlation matrix estimated from observations
#   - Schaake shuffle: reordering based on historical observations
#       SSh-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       SSh-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       SSh-R: random sample
#       SSh-S: stratifed sampling approach of Hu et al (2016)
#       SSh-T: see above, for gev0 currently set in comments

# input:
#   method: indicate reordering method: "none", "ECC", "dECC", "GCA", "SSh",
#   variant: variant of re-ordering method to be used:
#       Q: ...
#       R: ...
#   ensfc, ensfc_init, obs, obs_init: output of generators of observations and ensemble fcsts
#     dimensions need to match
#   postproc_out: output of postproc code in postproc_ensfc.R (array of EMOS parameter values)
#   EMOS_sample: give specific EMOS sample to the function to use in re-ordering method
#     to ensure that for the methods with randomness, the same EMOS samples are used
#     otherwise, obtaining a new EMOS-R sample as basis for ECC-R means that different
#     random samples are used in each margin
#     This input is not required, and only used in the run_all function for the simulation





mvpp <- function(method, variant = NULL, ens, verobs, postproc_out, EMOS_sample = NULL, ECC_out = NULL){


  require(evd)
  require(NORTARA)

  # include some checks for inputs

  # check whether the same model is used for input and output
  if(ens$model != verobs$model){
    stop("models for forecasts and observations do not match")
  }

  if(ens$model != postproc_out$model){
    stop("models for forecasts and post-processed forecasts do not match")
  }


  ensfc <- ens$ensfc
  ensfc_init <- ens$ensfc_init

  obs <- verobs$obs
  obs_init <- verobs$obs_init

  # generate array for ouput
  mvppout <- array(NA, dim = dim(ensfc))
  n <- dim(mvppout)[1]
  m <- dim(mvppout)[2]
  d <- dim(mvppout)[3]


  # Quantile function of GEV left-censored at 0
  qgev0 <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
  {
  require(evd)
  pmax(0, qgev(p, loc, scale, shape, lower.tail))
  }

  # CDF of GEV left-censored at 0
  pgev0 <- function(x, loc=0, scale=1, shape=0, lower.tail=TRUE)
  {
  require(evd)
  y <- numeric(length(x))
  id <- x>=0
  y[id] <- pgev(x[id], loc, scale, shape, lower.tail)
  y[!id] <- 0
  return(y)
  }
  
  

  # Ensemble mean difference (Scheuerer, 2014)
  gini.md <- function(x,na.rm=FALSE)
  {
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)]
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n*(n-1))-2*mean(x)*(n+1)/(n-1))
  }


    # EMOS methods: no accounting for multivariate dependencies, simply draw from marginals
    if(method == "EMOS"){
      if(!is.null(EMOS_sample)){
        stop("for method = 'EMOS', input 'EMOS_sample' must be 'NULL'")
      }
      if(is.null(variant)){
        stop("for method = 'EMOS', a variant parameter must be specified")
      }

      if(variant == "R"){
        for(nn in 1:n){
          for(dd in 1:d){
            
           par_scale <- postproc_out$parameters[nn, dd, 2]
           par_loc <- postproc_out$location[nn, dd]
           par_shape <- postproc_out$shape[dd]
            
            # Random number generation from left-censored GEV by applying quantile function
            # of left-censored GEV to uniform random numbers,
            # as currently no function such as "rgev0" exists, only rgev for the non-censored distribution
            mvppout[nn, , dd] <- qgev0(runif(m), loc=par_loc, scale=par_scale, shape=par_shape)
          }
        }
      } else if(variant == "Q"){
        qlevels <- 1:m/(m+1)
        for(nn in 1:n){
          for(dd in 1:d){
          
           par_scale <- postproc_out$parameters[nn, dd, 2]
           par_loc <- postproc_out$location[nn, dd]
           par_shape <- postproc_out$shape[dd]

            mvppout[nn, , dd] <- qgev0(qlevels, loc=par_loc, scale=par_scale, shape=par_shape)
          }
        }
      } else if(variant == "QO"){
        qlevels <- (1:m-0.5)/m
        for(nn in 1:n){
          for(dd in 1:d){
          
           par_scale <- postproc_out$parameters[nn, dd, 2]
           par_loc <- postproc_out$location[nn, dd]
           par_shape <- postproc_out$shape[dd]
           
            mvppout[nn, , dd] <- qgev0(qlevels, loc=par_loc, scale=par_scale, shape=par_shape)
          }
        }
      } else if(variant == "S"){
        breakpoints <- 0:m/m
        qlevels <- runif(m, min = breakpoints[1:m], max = breakpoints[2:(m+1)])
        for(nn in 1:n){
          for(dd in 1:d){
          
           par_scale <- postproc_out$parameters[nn, dd, 2]
           par_loc <- postproc_out$location[nn, dd]
           par_shape <- postproc_out$shape[dd]

            mvppout[nn, , dd] <- qgev0(qlevels, loc=par_loc, scale=par_scale, shape=par_shape)
          }
        }
      }
#      else if(variant == "T"){
#        for(nn in 1:n){
#          for(dd in 1:d){
#            ensfc_tmp <- ensfc[nn, , dd]
#
#           par_scale <- postproc_out$parameters[nn, dd, 2]
#           par_loc <- postproc_out$location[nn, dd]
#           par_shape <- postproc_out$shape[dd]
#
#            # Computing ensemble variance or mean difference can be an issue when all members are 0
#            # Furthermore, mean and variance of ensemble not ideal estimates for location and scale in GEV
#            ens_par <- c(mean(ensfc_tmp), gini.md(ensfc_tmp))
#            qlevels <- pgev0(ensfc_tmp, loc = ens_par[1], scale = ens_par[2], shape=0)
#
#            mvppout[nn, , dd] <- qgev0(qlevels, loc=par_loc, scale=par_scale, shape=par_shape)
#          }
#        }
#      }
      # end of EMOS methods
    }


    # GCA code
    if(method == "GCA"){
      require(MASS)

      # if no EMOS_sample to base GCA on is given, recursively call 'mvpp' to generate such a sample
      if(any(!is.null(c(EMOS_sample, variant)))){
        message("'EMOS_sample' and 'variant' input have no effect for GCA")
      }


    # concatenate obs_init and obs arrays to determine correlation matrix for Gaussian copulas
    obs_all <- rbind(obs_init, obs)

    for(nn in 1:n){
    
      # select forecast cases to estimatate correlation matrix
      #   ... theoretically, correlation matrix could also be estimated from past ensemble forecasts
      #   ... to speed up computation, it would also be possible to only use obs_init

      obs_train <- obs_all[1:(dim(obs_init)[1]+nn-1), ]
      # transform to Gaussian
      obs_train_transformed <- obs_train
      for(dd in 1:d){
      
      par_scale <- postproc_out$parameters[nn, dd, 2]
      par_loc <- postproc_out$location[nn, dd]
      par_shape <- postproc_out$shape[dd]
      
      # Value of CDF at 0
      F0 <- pgev0(0, loc=par_loc, scale=par_scale, shape=par_shape)
      # Randomization of value CDF takes at 0, by drawing from uniform distribution on (0, F(0))
      u <- if(obs_train[nn,dd]==0) runif(1, 0, F0) else pgev0(obs_train[nn,dd], loc=par_loc, scale=par_scale, shape=par_shape)
      
      # tolerance query to ensure that u in (0,1)
      u <- if(u > 0.9999999999999999) 0.9999999999999999 else u
      u <- if(u < 0.0000000000000001) 0.0000000000000001 else u
      
      # Observations transformed to standard normal space
      obs_train_transformed[nn, dd] <- qnorm(u)
      }
      # estimate correlation matrix
      cor_obs <- cor(obs_train_transformed)
      # draw random sample from multivariate normal distribution with this correlation matrix
      mvsample <- mvrnorm(n = m, mu = rep(0,d), Sigma = cor_obs)
      
      # impose dependence structure on post-processed forecasts
      for(dd in 1:d){
        par_scale <- postproc_out$parameters[nn, dd, 2]
        par_loc <- postproc_out$location[nn, dd]
        par_shape <- postproc_out$shape[dd]
        mvppout[nn, , dd] <- qgev0(pnorm(mvsample[, dd]), loc=par_loc, scale=par_scale, shape=par_shape)
      }
    }
      # end of GCA code
    }


    # in random methods: distinguish cases with and without given EMOS_sample, maybe only handle that with sample at first, rest can be included later on
    # if no sample is given, a new one has to be generated, as done for the EMOS methods themselves


  # ECC code
  if(method == "ECC"){
    # if no EMOS_sample to base ECC on is given, recursively call 'mvpp' to generate such a sample
    if(is.null(EMOS_sample)){
      EMOS_sample <- mvpp(method = "EMOS", variant = variant, postproc_out = postproc_out,
                          ens = ens, verobs = verobs)
      message("no 'EMOS_sample' given for ECC, a new one is generated")
    } else if(!is.null(variant)){
      message("'variant' parameter has no influence if EMOS_sample is supplied,
              make sure the EMOS_sample is produced with the desired variant")
    }

    # application of ECC is independent of 'variant' parameter, dependency is only through EMOS_sample
    for(nn in 1:n){
      for(dd in 1:d){
        ensfc_tmp <- ensfc[nn, , dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(ensfc_tmp, ties.method = "random")]
      }
    }
    # end of ECC code
    }

  # Schaake shuffle code
  if(method == "SSh"){
    # if no EMOS_sample to base SSh on is given, recursively call 'mvpp' to generate such a sample
    if(is.null(EMOS_sample)){
      if(is.null(variant)){
        stop("if no 'EMOS_sample' is given, 'variant' has to be specified to generate a new one")
      }
      EMOS_sample <- mvpp(method = "EMOS", variant = variant, postproc_out = postproc_out,
                          ens = ens, verobs = verobs)
      message("no 'EMOS_sample' given for SSh, a new one is generated")
    } else if(!is.null(variant)){
      message("'variant' parameter has no influence if EMOS_sample is supplied,
                make sure the EMOS_sample is produced with the desired variant")
    }

    # concatenate obs_init and obs arrays to sample from available forecast cases later on
    obs_all <- rbind(obs_init, obs)

    # reorder post-processed forecast sample according to past observations
    
    for(nn in 1:n){
    
      # choose set of past forecast cases to determine dependence template
      #   ... this way, a new set of IDs is drawn for every forecast instance
      #   ... this needs to depend on nn in a more suitable manner if there is temporal change in the simulation setup
      
      obs_IDs <- sample(x = 1:(dim(obs_init)[1]+nn-1), size = m, replace = FALSE)
      for(dd in 1:d){
        obs_tmp <- obs_all[obs_IDs, dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(obs_tmp, ties.method = "random")]
      }
    }
    # end of SSh code
  }

  # dECC
  if(method == "dECC"){
    # both EMOS_sample and ECC_out are required in the algorithm later on
    # if no EMOS_sample to base dECC on is given, recursively call 'mvpp' to generate such a sample
    if(is.null(EMOS_sample)){
      EMOS_sample <- mvpp(method = "EMOS", variant = variant, postproc_out = postproc_out,
                          ens = ens, verobs = verobs)
      message("no 'EMOS_sample' given for dECC, new one is generated, requires 'variant'")
    }

    # additionally, dECC requires the application of ECC beforehand; ECC_out is utilized here
    if(is.null(ECC_out)){
      message("no 'ECC_out' given, generating new one from input EMOS_sample")
      if(!is.null(variant)){message("input 'variant' has no effect here,
                                    make sure variant in EMOS_sample is desired one")}
      ECC_out <- mvpp(method = "ECC", postproc_out = postproc_out,
                      ens = ens, verobs = verobs, EMOS_sample = EMOS_sample)
      }

    # estimate correlation matrix
    #   .... for now, only from the init period,
    #   .... later on, it should be possible to use a shifting window
    #   .... therefore, use something like
    #           require(abind)
    #           ensfc_all <- abind(ensfc_init, ensfc, along = 1)
    #   .... to generate required arrays

    # e in the dECC paper, eq (14)
    fcerror <- obs_init - apply(ensfc_init, c(1,3), mean)
    # compute error correlation matrix and root
    R_e <- cor(fcerror)
    eig <- eigen(R_e)
    Re_root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

    # c in the dECC paper, eq (18)
    error_cor <- ECC_out - ensfc

    # breve c, eq (19)
    breve_c <- array(NA, dim = dim(ensfc))
    for(nn in 1:n){
      for(mm in 1:m){
        breve_c[nn, mm, ] <- Re_root %*%  error_cor[nn, mm, ]
      }
    }

    # adjusted ensemble breve x, eq (22)
    breve_x <- ensfc + breve_c

    # apply ECC with breve_x as dependence template (Step 5 in paper)
    for(nn in 1:n){
      for(dd in 1:d){
        ensfc_tmp <- breve_x[nn, , dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(ensfc_tmp, ties.method = "random")]
      }
    }

    # end of dECC code
  }

  return(mvppout)
}


