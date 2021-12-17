#' @name normMixEM
#' @title EM Algorithm for Mixtures of Independent Multivariate Normals
#' @description  This is a normMixEM class, which initializes and runs EM algorithm output for mixtures of independent multivariate normal distributions
#' @param input_data  A matrix/dataframe of size nxp consisting of the data.
#' @param num_components Number of components.
#' @return A list including
#' \itemize{
#'   \item convergence - Convergece status. 0 if converged,1 if not.
#'   \item mu_mat - A k*p matrix includes the mean parameters of k normals.
#'   \item sigma_mat - A k*p matrix includes the variance parameters of k normals.
#'   \item pi_vec - A vector of length k includes the probability parameters of k components.
#'   \item iter - A integer. Number of iterations
#'   \item loglik_list - A number. The log-likelihood of the output result
#'   \item prob_mat - A n*k matrix includes the probabilities of being assigned to k components for each sample.
#'}
#' @examples
#'  EM <- normMixEM$new(input_dat = data,num_components = num_components) ## initialization with $new
#'  res_EM <- EM$run.EM(loglik_tol=1e-5) ## run the main function for EM algorithm with $run.EM
#' @import stats
#' @import mixtools
#' @import mvtnorm
#' @export



#' Create a class normMixEM
#' method to initialize normMixEM class (called with normMixEM$new)
normMixEM <- setRefClass("normMixEM",
                         fields=list(k = "integer",
                                     n = "integer",
                                     dat_mat = "matrix",
                                     pi_vec = "vector",
                                     mu_mat = "matrix",
                                     sigma_mat = "matrix",
                                     prob_mat = "matrix",
                                     loglik = "numeric",
                                     tol = "numeric"))

normMixEM$methods(initialize = function(input_dat,num_components){
  dat_mat <<- input_dat ## Use <<- to assign fields
  k <<- num_components
  n <<- dim(dat_mat)[1]
  p <<- dim(dat_mat)[2]
})

# dat_mat: n * p 
# mu_mat: k * p
# sigma_mat: k * p
# prob_mat: n * k

#' method to initialize parameters for E-M
normMixEM$methods(init.paras = function(){
  tol <<- 1e-100  ## small number to avoid log(zero) 
  init   <<- normalmix_init(dat_mat, k = k)
  pi_vec <<- init$lambda  # k*1
  mu_mat <<- init$mu # k*p 
  sigma_mat <<- init$sigma # k*p
  prob_mat <<- matrix(NA,nrow=n,ncol=k)
  loglik <<- -Inf
})


#' E-step for E-M algorithm
normMixEM$methods(update.prob = function(){
  ## prob_mat contains log-likelihood of each components
  # version 1
  prob_mat <<- sapply(1:k, function(j)
    log(pi_vec[j]+tol) + logdmvnorm(dat_mat, mu_mat[j,], diag(sigma_mat[j,]))
  )
  ## version 2
  # prob_mat <<- sapply(1:k, function(j)
  #   log(pi_vec[j]+tol) + log(mvtnorm::dmvnorm(dat_mat, mu_mat[j,], diag(sigma_mat[j,]))+tol)
  # ) #  mvtnorm::dmvnorm is more stable!!
  max_log_prob <- apply(prob_mat,1,max) ## local change, do not change the field value
  prob_mat <<- exp(prob_mat - max_log_prob) ## re-scale probability (important)
  sum_prob <- apply(prob_mat,1,sum)   ## sum probabilities for each observation
  prob_mat <<- prob_mat/sum_prob      ## normalize the probability E[z|theta]
  loglik <<- sum(max_log_prob + log(sum_prob))  ## evaluate log-likelihood
})

#' M-step for E-M algorithm to update pi_vec
normMixEM$methods(update.pi = function(){
  pi_vec <<- apply(prob_mat,2,mean)# prob_mat:  n*k => k*1
})

#' M-step for E-M algorithm to update mu_mat
normMixEM$methods(update.mu = function(){
  # mu_mat: k * p
  mu_mat <<- (crossprod(prob_mat,dat_mat)/(apply(prob_mat,2,sum)+tol))#apply(dat_mat*prob_mat,2,mean)/(pi_vec+tol)
})

#' M-step for E-M algorithm to update sigma_mat
normMixEM$methods(update.sigma = function(){
  # sigma_mat: k * p
  sigma_mat <<- t(sapply(1:k, function(i) crossprod(prob_mat[,i], 
                                                    t(t(dat_mat)- mu_mat[i,])^2)
  ))
  ## update sigma as weighted average
  sigma_mat <<- sigma_mat/(n*pi_vec+tol)
})


#' check convergence of mixture normal EM
normMixEM$methods(check.tol = function(fmax,fmin,ftol){
  delta = abs(fmax - fmin)
  accuracy = (abs(fmax) + abs(fmin))*ftol
  return(delta < (accuracy + tol))
})

#' main function for E-M algorithm
normMixEM$methods(run.EM = function(max_iter=1000L,loglik_tol=1e-5){
  convergence = 1L
  init.paras()  ## initialize parameter
  loglik_list = NULL
  for(iter in 1:max_iter){
    loglik0 <- loglik ## log-likelihood of previous steps
    update.prob() # E-step
    update.pi()   # M-step for pi_vec
    update.mu()   # M-step for mu_mat
    update.sigma() # M-step for sigma_mat
    loglik_list = c(loglik_list,loglik) # append log-likelihood
    if(check.tol(loglik0,loglik,loglik_tol)){
      convergence = 0 # converged
      break
    }
  }
  return(list(convergence=convergence,mu_mat=mu_mat,sigma_mat=sigma_mat,
              pi_vec=pi_vec,iter=iter,loglik_list=loglik_list, prob_mat=prob_mat))
})



#'@name normalmix_init
#'@title Initialization for normMixEm Algorithm
#'@description Internal initialization function for normMixEm algorithm in this package. Modified based on the normalmix.init{mixtools}.
#'
#'@param x A matrix of size n*p consisting of the data.
#'@param k Number of components. Default as 5.
#'
#'@return A list of parameters for initialization including:
#' \itemize{
#'   \item lambda - A vector of length k includes the probability parameters of k components.
#'   \item mu - A k*p matrix includes the mean parameters of k Normals.
#'   \item sigma - A k*p matrix includes the variance parameters of k Normals.
#'   \item k - Integer. Number of components.  
#'}
#'
normalmix_init = function (x, k = 5) 
{
  n <- nrow(x)
  p <- ncol(x)
  y <- apply(x, 1, mean)
  x <- x[order(y), ]
  x.bin <- list()
  for (j in 1:k) {
    x.bin[[j]] <- x[max(1, floor((j - 1) * n/k)):ceiling(j *  n/k), ]# get part of the rows 
  }
  # init sigma: based on random exp
  # sigma: k*p
  sigma.hyp = lapply(1:k, function(i) (apply(x.bin[[i]],2, var))^(-1))
  sigma = t(sapply(1:k, function(i) 1/rexp(p, rate = sigma.hyp[[i]]))
  )# i * p*p 
  
  # init mu: based on random normal
  #***
  mu.hyp <- lapply(1:k, function(i) apply(x.bin[[i]], 2,mean))
  # mu: k*p
  mu <- t(sapply(1:k, function(i) as.vector(rmvnorm_chol(1,
                                                         mu = as.vector(mu.hyp[[i]]),
                                                         sigma = as.matrix(diag(sigma[i,]))
  )
  )
  )
  )
  # init lambda
  lambda <- runif(k)
  lambda <- lambda/sum(lambda)
  return(list(lambda = lambda, mu = mu, sigma = sigma, k = k))
}


#'@name rmvnorm_chol
#'@title Simulate from a Multivariate Normal Distribution
#'@description Sample from a multivariate normal distribution based on Cholesky decomposition. Modified based on rmvnorm {mixtools}
#'
#'@param n Number of vectors to simulate.
#'@param mu Mean vector.
#'@param sigma Covariance matrix, assumed symmetric and nonnegative definite.
#'
#'@return A n*p matrix. Each row is a multivariate normal vector generated independently based on the input parameters.
#'
#'

# Alternative version of rmvnorm to eliminate dependence of mixtools
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
