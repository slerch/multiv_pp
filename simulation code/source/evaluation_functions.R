## evaluation functions, operating on output of mvpp (or ensfc)

# document..

### scoring rules
# crps
# es
# vs, p = 0.5
# vs, p = 1
# vs_w, p = 0.5
# vs_w, p = 1
#   for 1 fixed pre-specified weighting matrix

## -------------------- ##

## Part 1: multivariate scoring rules

# wrapper function for scoringRules functions es_sample and vs_sample 
#   to apply to specific array formats in simulation study

# code requires scoringRules package, version >= 0.9.2

# input:
#   mvpp_out: output of mvpp() function
#   obs: observations in evaluation period, should match mvpp_out (is not checked)
#   weight: (only in vs): logical whether to use weighting 
#       (this will use a pre-defined weight matrix, see function)
#   p: (only in vs): parameter p in VS

# output:
#   vector of scores (length = nout)

es_wrapper <- function(mvpp_out, obs){
  n <- dim(mvpp_out)[1]
  out <- vector(length = n)
  for(nn in 1:n){
    out[nn] <- es_sample(y = obs[nn,], dat = t(mvpp_out[nn,,]))
  }
  out
}

vs_wrapper <- function(mvpp_out, obs, weight, p){
  if(!is.logical(weight)){stop("'weight' must be logical")}
  d <- dim(mvpp_out)[3]
  if(weight){
    w <- matrix(0, d, d)
    for(i in 1:d){
      for(j in 1:d){
        w[i,j] <- 1/abs(i-j)
      }
    }
    diag(w) <- 0
    nf <- 0.5*sum(w)
    w <- 1/nf*w
  } else{
    w <- NULL
  }
  n <- dim(mvpp_out)[1]
  out <- vector(length = n)
  for(nn in 1:n){
    out[nn] <- vs_sample(y = obs[nn,], dat = t(mvpp_out[nn,,]), w = w, p = p)
  }
  out
}

## -------------------- ##

## Part 2: univariate scoring rules (CRPS)

# wrapper function for scoringRules function crps_sample
#   to apply to specific array formats in simulation study
# code requires scoringRules package, version >= 0.9.2

# input:
#   mvpp_out: output of mvpp() function
#   obs: observations in evaluation period, should match mvpp_out (is not checked)

# output:
#   array of scores (dimensions = forecast case (n); dimension (d))

crps_wrapper <- function(mvpp_out, obs){
  n <- dim(mvpp_out)[1]
  d <- dim(mvpp_out)[3]
  out <- matrix(NA, nrow = n, ncol = d)
  for(nn in 1:n){
    for(dd in 1:d){
      out[nn, dd] <- crps_sample(y = obs[nn,dd], dat = mvpp_out[nn, , dd])
    }
  }
  out
}

