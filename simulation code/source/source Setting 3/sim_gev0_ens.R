# -------------------------------------------------------------------------- #

# generate multivariate ensemble forecasts from
# a distribution with univariate GEV0 margins
# Note that parts of the notation may differ from the notation in the paper


# input:
#   model: character string, indicating which model is used for the observations
#          currently for Setting 3 GEV0 there is a separate set of functions
#          where string "gev0" needs to be set
#   nout: number of multivariate observations to be generated as evaluation period
#   ninit: additional initial training period for model estimation purposes
#   nmembers: number of ensemble members
#   Additional parameters, depending on the chosen model:
#   For "gev0" this is location, scale and shape parameter
#   Furthermore, the correlation parameter (for observations rho0)
#   defining the multivariate correlation matrix needs to be specified (see paper)


# output:
#   list of two arrays "ensfc_init" and "ensfc" containing the ensemble forecasts
#   dimensions: ninit, nmembers, d; and nout, nmembers, d
#   content in array dimensions:
#     first dimension: forecast instance
#     second dimension: ensemble member
#     third dimension: dimension in multivariate setting


# -------------------------------------------------------------------------- #
library(MASS)
library(evd)
library(NORTARA)
# -------------------------------------------------------------------------- #
# quantile function of GEV left-censored at 0
qgev0 <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
  require(evd)
  pmax(0, qgev(p, loc, scale, shape, lower.tail))
}
# -------------------------------------------------------------------------- #
generate_ensfc <- function(model = "gev0", nout, ninit, nmembers, d, 
                           rho = 0.5, 
                           gev0loc = 0, gev0scale = 1, gev0shape = 0){
  m <- nmembers
  
  # initialize output arrays
  ensfc_init <- array(NA, dim = c(ninit, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  
  # Correlation matrix,
  cor.mat <- toeplitz(rho^(0:(d-1)))
  
  # input of genNORTARA, which allows to simulate from a multivariate
  # distribution with given margins and given target-correlation structure
  # requirements are lists for the different dimensions of specific types,
  # constructed below
  # Currently: all margins identical, namely gev0 (GEV left-censored at 0)
  margin_distr <- rep("qgev0", d)  
  para_list <- vector(mode = "list", length = d)
  
  # The list places need to be named "m1", "m2", ... , "md" for input in the nortara function
  # (this names "m.." are not to be confused with the variable m for the ensemble size
  names(para_list) <- paste0("m", 1:d) 
  para_list <- lapply(para_list, function(x) {
    x <- list(loc = gev0loc, scale = gev0scale, shape = gev0shape)})



  # generate the multivariate ensemble members with above specifications
  # Instead of loop, generate all m d-dim. members at the same time
  # otherwise genNORTARA very slow as covariance matrix
  # has to be recompted in every loop step, which is not necessary

  tmp <- genNORTARA((m * ninit), cor_matrix = cor.mat,
                    invcdfnames = margin_distr, 
                    paramslists = para_list)
                    
  # now arrange the output into m multivariate members of dimension d
  rep_ind <- rep(1:ninit, each = m)
  for(nn in 1:ninit){ensfc_init[nn,,] <- tmp[rep_ind == nn, ]}
  

  tmp <- genNORTARA((m * nout), cor_matrix = cor.mat,
                    invcdfnames = margin_distr, 
                    paramslists = para_list)
                    
   # now arrange the output into m multivariate members of dimension d
  rep_ind <- rep(1:nout, each = m)
  for(nn in 1:nout){ensfc[nn,,] <- tmp[rep_ind == nn, ]}




  return(list("ensfc_init" = ensfc_init, 
              "ensfc" = ensfc, 
              "model" = 'gev0', 
              "loc" = gev0loc, 
              "scale" = gev0scale, 
              "shape" = gev0shape))
}
