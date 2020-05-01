# code to generate multivariate observations 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: numeric, indicating which model is used for the observations
#   nout: number of multivariate observations to be generated as evaluation period 
#   ninit: additional initial training period for model estimation purposes
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "obs_init" and "obs" containing the observations
#   dimensions: ninit, d; and nout; d
#   first dimension: forecast instance; second dimension: dimension in multivariate setting


generate_obs <- function(model, nout, ninit, d, ...){
  require(MASS)
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, d)))){
    stop("Input 'nout', 'ninit' and 'd' need to be numeric of length 1")
  }
  
  # initialize output arrays
  obs_init <- array(NA, dim = c(ninit, d))
  obs <- array(NA, dim = c(nout, d))
  
  # Setting 1 (multivariate Gaussian distribution)
  if(model == 1){
    # check if appropriate additional parameters are given
    input <- list(...)
    ind_check <- match(c("rho0"), names(input), nomatch = 0)
    if(ind_check == 0){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste("rho0", collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign additional parameters from input
    rho0 <- input$rho0
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- rho0^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(0, d)
    
    # observations
    obs_init <- mvrnorm(n = ninit, mu = mu, Sigma = S)
    obs <- mvrnorm(n = nout, mu = mu, Sigma = S)
  }
  
  
  # Setting 2 (S1 in the paper)
  if(model == 2){
    
    require(tmvtnorm)
    input <- list(...)
    
    
    # assign additional parameters from input
    rho0 <- input$rho0
    
    if (is.null(input$mu0)) {
      mu <- rep(0, d)
    } else {
      mu <- rep(input$mu0, d)
    }
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- rho0^(abs(i-j))
      }
    }
    
    # specify lower and upper bounds of truncation
    if (is.null(input$lower)){
      lowerB <- rep(-Inf,d)
    } else {
      lowerB <- input$lower
    }
    
    if (is.null(input$upper)){
      upperB <- rep(Inf,d)
    } else {
      upperB <- input$upper
    }
    
    # observations
    obs_init <- rtmvnorm(n = ninit, mean = mu, sigma = S, lower = lowerB, upper = upperB)
    obs <- rtmvnorm(n = nout, mean = mu, sigma = S, lower = lowerB, upper = upperB)
    return(list("obs_init" = obs_init, "obs" = obs, "model" = 'tnormal', "lower" = lowerB, "upper" = upperB))
    
  }
  
  # Setting 4 (3A in the paper)
  if(model == 4){
    # check if appropriate additional parameters are given
    input <- list(...)
    ind_check <- match(c("rho0"), names(input), nomatch = 0)
    if(ind_check == 0){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste("rho0", collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign additional parameters from input
    rho0 <- input$rho0
    
    n_all <- ninit + nout
    
    # generate forecasts
    for(nn in 1:ninit){
      
      # mean vector
      mu <- rep(sin(2*pi*nn/n_all), d)
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho0^(abs(i-j)) + sin(2*pi*nn/n_all)
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR)
      
      tmp <- mvrnorm(n = 1, mu = mu, Sigma = S)
      obs_init[nn,] <- tmp
    }
    
    for(nn in 1:nout){
      
      # mean vector
      mu <- rep(sin(2*pi*nn/n_all), d)
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho0^(abs(i-j)) + sin(2*pi*(ninit + nn)/n_all)
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = 1, mu = mu, Sigma = S)
      obs[nn,] <- tmp
    }
    
  }
  
  # Setting 5 (3B in the paper)
  if(model == 5){
    # check if appropriate additional parameters are given
    input <- list(...)
    ind_check <- match(c("rho0"), names(input), nomatch = 0)
    if(ind_check == 0){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste("rho0", collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign additional parameters from input
    rho0 <- input$rho0
    a <- input$a
    
    n_all <- ninit + nout
    
    # generate forecasts
    for(nn in 1:ninit){
      
      # mean vector
      mu <- rep(0, d)
      
      # time-varying rho
      rho0_nn <- rho0*(1 - a/2) + rho0*a/2*sin(2*pi*nn/n_all) 
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho0_nn^(abs(i-j)) 
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR)
      
      tmp <- mvrnorm(n = 1, mu = mu, Sigma = S)
      obs_init[nn,] <- tmp
    }
    
    for(nn in 1:nout){
      
      # mean vector
      mu <- rep(0, d)
      
      # time-varying rho
      rho0_nn <- rho0*(1 - a/2) + rho0*a/2*sin(2*pi*nn/n_all) 
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho0_nn^(abs(i-j)) 
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = 1, mu = mu, Sigma = S)
      obs[nn,] <- tmp
    }
    
  }
  
  return(list("obs_init" = obs_init, "obs" = obs))
}
