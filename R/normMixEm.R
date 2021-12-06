library(mixtools)
library(stats)
library(microbenchmark)
library(dplyr)
library(devtools)

#### DONNOT MODIFY THIS PART 
#### RUN THIS CHUNK BEFORE USING 'normMixEm ' FOR CLUSTERING

## My EM algorithm: Multivariate(independent) norm()
#' Create a class normMixEM
normMixEm <- setRefClass("normMixEM",
                         fields=list(k = "integer",
                                     n = "integer",
                                     dat_mat = "matrix",
                                     pi_vec = "vector",
                                     mu_mat = "matrix",
                                     sigma_mat = "matrix",
                                     prob_mat = "matrix",
                                     loglik = "numeric",
                                     tol = "numeric"))


#' method to initialize normMixEm class (called with normMixEm$new)
#' @param input_data - data to initalize
#' @param num_components - number of components (k)
normMixEm$methods(initialize = function(input_dat,num_components){
  dat_mat <<- input_dat ## Use <<- to assign fields
  k <<- num_components
  n <<- dim(dat_mat)[1]
  p <<- dim(dat_mat)[2]
})

# dat_mat: n * p 
# mu_mat: k * p
# sigma_mat: k * p
# prob_mat: n * k
# rmvnorm_chol
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

#' method to initialize parameters for E-M
normMixEm$methods(init.paras = function(){
  tol <<- 1e-100  ## small number to avoid log(zero) 
  init   <<- normalmix_init(dat_mat, k = k)
  pi_vec <<- init$lambda  # k*1
  mu_mat <<- init$mu # k*p 
  sigma_mat <<- init$sigma # k*p
  prob_mat <<- matrix(NA,nrow=n,ncol=k)
  loglik <<- -Inf
})


#' E-step for E-M algorithm
normMixEm$methods(update.prob = function(){
  prob_mat <<- sapply(1:k, function(j)
    log(pi_vec[j]+tol) + logdmvnorm(dat_mat, mu_mat[j,], diag(sigma_mat[j,]))
  ) #  logdmvnorm(dat_mat, mu_mat[j,] vector, diag(sigma_mat[j,]) matrix)
  
  max_log_prob <- apply(prob_mat,1,max) ## local change, do not change the field value
  prob_mat <<- exp(prob_mat - max_log_prob) ## re-scale probability (important)
  sum_prob <- apply(prob_mat,1,sum)   ## sum probabilities for each observation
  prob_mat <<- prob_mat/sum_prob      ## normalize the probability E[z|theta]
  loglik <<- sum(max_log_prob + log(sum_prob))  ## evaluate log-likelihood
})

#' M-step for E-M algorithm to update pi_vec
normMixEm$methods(update.pi = function(){
  pi_vec <<- apply(prob_mat,2,mean)# prob_mat:  n*k => k*1
})

#' M-step for E-M algorithm to update mu_mat
normMixEm$methods(update.mu = function(){
  # mu_mat: k * p
  mu_mat <<- (crossprod(prob_mat,dat_mat)/(apply(prob_mat,2,sum)+tol))#apply(dat_mat*prob_mat,2,mean)/(pi_vec+tol)
})

#' M-step for E-M algorithm to update sigma_mat
normMixEm$methods(update.sigma = function(){
  # sigma_mat: k * p
  sigma_mat <<- t(sapply(1:k, function(i) crossprod(prob_mat[,i], 
                                                    t(t(dat_mat)- mu_mat[i,])^2)
  ))
  ## update sigma as weighted average
  sigma_mat <<- sigma_mat/(n*pi_vec+tol)
})


#' check convergence of mixture normal EM
normMixEm$methods(check.tol = function(fmax,fmin,ftol){
  delta = abs(fmax - fmin)
  accuracy = (abs(fmax) + abs(fmin))*ftol
  return(delta < (accuracy + tol))
})
#' main function for E-M algorithm

normMixEm$methods(run.EM = function(max_iter=1000L,loglik_tol=1e-5){
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




