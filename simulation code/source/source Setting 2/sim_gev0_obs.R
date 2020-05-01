# ---------------------------------------------------------------------------- #

# code to generate multivariate observations from
# a distribution with univariate GEV0 margins
# Note that parts of the notation may differ from the notation in the paper

# input:
#   model: character string, indicating which model is used for the observations
#          currently for Setting 3 GEV0 there is a separate set of functions
#          where string "gev0" needs to be set
#   nout: number of multivariate observations to be generated as evaluation period
#   ninit: additional initial training period for model estimation purposes
#   d: dimension of the multivariate vectors
#   Additional parameters, depending on the chosen model:
#   For "gev0" this is location, scale and shape parameter
#   Furthermore, the correlation parameter (for observations rho0)
#   defining the multivariate correlation matrix needs to be specified (see paper)

# output:
#   list of two arrays "obs_init" and "obs" containing the observations
#   dimensions: ninit, d; and nout; d
#   first dimension: forecast instance; second dimension: dimension in multivariate setting


# --------------------------------------------------------------------------- #
library(MASS)
library(evd)
library(NORTARA)
# --------------------------------------------------------------------------- #
# quantile function of GEV left-censored at 0
qgev0 <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
  require(evd)
  pmax(0, qgev(p, loc, scale, shape, lower.tail))
}
# -------------------------------------------------------------------------- #
generate_obs <- function(model = "gev0", nout, ninit, d, rho0, 
                         gev0loc = 0, gev0scale = 1, gev0shape = 0)
{
  # initialize output arrays
  obs_init <- array(NA, dim = c(ninit, d))
  obs <- array(NA, dim = c(nout, d))
  # Correlation matrix
  S <- toeplitz(rho0^(0:(d-1)))
  
  # Define input of genNORTARA, which allows to simulate from a multivariate
  # distribution with given margins and given target-correlation structure
  # requirements are lists for the different dimensions of specific types,
  # constructed below
  # Currently: all margins identical, namely gev0 (GEV left-censored at 0)
  margin_distr <- rep("qgev0", d)
  para_list <- vector(mode = "list", length = d)
  
  # The list places need to be named "m1", "m2", ... , "md" for input in the nortara function
  # (this names "m.." are not to be confused with the variable m for the ensemble size
  names(para_list) <- paste0("m", 1:d)
  
  para_list <- lapply(para_list, function(x){
    x <- list(loc = gev0loc, scale = gev0scale, shape = gev0shape)})
    
  # generate the multivariate observations with above specifications
  
  obs_init <- genNORTARA(ninit, cor_matrix = S,
                      invcdfnames = margin_distr,
                      paramslists = para_list)
  obs <- genNORTARA(nout, cor_matrix = S,
                    invcdfnames = margin_distr,
                    paramslists = para_list)
                    
                    
  return(list("obs_init" = obs_init, "obs" = obs, 
              "model" = 'gev0', "loc" = gev0loc, 
              "scale" = gev0scale, "shape" = gev0shape))
}
