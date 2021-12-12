
#Sample from a multivariate normal distribution based on cholesky decomposition
# Alternative version of rmvnorm to eliminate dependence of mixtools
# on additional package 'mvtnorm'
#Sample from a multivariate normal distribution based on cholesky decomposition

rmvnorm_chol <- function(n, mu=NULL, sigma=NULL) {
  if (is.null(mu)) {
    if (is.null(sigma)) {
      return(rnorm(n)) # return standard norm
    } else {
      mu = rep(0, nrow(sigma))
    }
  } else if (is.null(sigma)) {
    sigma=diag(length(mu))
  }
  p <- length(mu)
  if (p != nrow(sigma) || p != ncol(sigma)) 
    stop("length of mu must equal nrow and ncol of sigma")
  Z <- matrix(rnorm(n*p),nrow=p,ncol=n)
  try(U <- chol(sigma))
  X <- mu + crossprod(U, Z)
  return(t(X))
}
