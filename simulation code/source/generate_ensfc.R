# generate multivariate ensemble forecasts 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: character string; indicating which model is used for generating the forecasts
#   nout: number of multivariate observations to be generated as evaluation period
#   ninit: additional initial training period for model estimation purposes
#   nmembers: number of ensemble members
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "ensfc_init" and "ensfc" containing the ensemble forecasts
#   dimensions: ninit, nmembers, d; and nout, nmembers, d
#   content in array dimensions:
#     first dimension: forecast instance 
#     second dimension: ensemble member
#     third dimension: dimension in multivariate setting


generate_ensfc <- function(model, nout, ninit, nmembers, d, ...){
  require(MASS)
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, nmembers, d)))){
    stop("Input 'nout', 'ninit', 'nmembers' and 'd' need to be numeric of length 1")
  }
  m <- nmembers
  
  # initialize output arrays
  ensfc_init <- array(NA, dim = c(ninit, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  
  # Setting 1
  if(model == 1){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- sigma*rho^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(eps, d)
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
  }
  
  # Setting 2 (S1 in the paper)
  if(model == 2){
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- sigma*rho^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(eps, d)
  
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
    
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- rtmvnorm(n = m, mean = mu, sigma = S, lower = lowerB, upper = upperB)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- rtmvnorm(n = m, mean = mu, sigma = S, lower = lowerB, upper = upperB)
      ensfc[nn,,] <- tmp
    }
    return(list("ensfc_init" = ensfc_init, "ensfc" = ensfc, "model" = 'tnormal', "lower" = lowerB, "upper" = upperB))

  }
  
  # Setting 4 (3A in the paper)
  if(model == 4){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    
    n_all <- ninit + nout
    
    # generate forecasts
    for(nn in 1:ninit){
      
      # mean vector
      mu <- rep(sin(2*pi*nn/n_all) + eps, d)
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho^(abs(i-j)) + sin(2*pi*nn/n_all)
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
      
    }
    
    
    for(nn in 1:nout){
      
      # mean vector
      mu <- rep(sin(2*pi*nn/n_all) + eps, d)
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- rho^(abs(i-j)) + sin(2*pi*(ninit + nn)/n_all)
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
    
  }
  
  # Setting 5 (3B in the paper)
  if(model == 5){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    a <- input$a
    
    n_all <- ninit + nout
    
    # generate forecasts
    for(nn in 1:ninit){
      
      # mean vector
      mu <- rep(eps, d)
      
      # time-varying rho
      rho_nn <- rho*(1 - a/2) + rho*a/2*sin(2*pi*nn/n_all) 
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- sigma*rho_nn^(abs(i-j)) 
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
      
    }
    
    for(nn in 1:nout){
      
      # mean vector
      mu <- rep(eps, d)
      
      # time-varying rho
      rho_nn <- rho*(1 - a/2) + rho*a/2*sin(2*pi*nn/n_all) 
      
      # covariance matrix
      R <- matrix(NA, d, d)
      for(i in 1:d){
        for(j in 1:d){
          R[i,j] <- sigma*rho_nn^(abs(i-j)) 
        }
      }
      RtR <- R %*% t(R)
      S <- cov2cor(RtR) 
      
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
    
  }

  return(list("ensfc_init" = ensfc_init, "ensfc" = ensfc))
}
