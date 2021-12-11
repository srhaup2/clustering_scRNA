library(devtools)
library(microbenchmark)
library(ggplot2)
library(sparsepca)
devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp)

#Simulate 1000x10 Data
m <- 10000
V1 <- rnorm(m, -100,200)
V2 <- rnorm(m, -100, 300)
V3 <- -0.1*V1 + 0.1*V2 + rnorm(m, 0, 100)

X <- cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X <- X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))

#VS prcomp
prcomp_out = prcomp(X, rank. = 3)
rcpp_out = spcaRcpp(X, k = 3, alpha = 0, beta = 0)   #set alpha & beta to 0

#print results
prcomp_out
rcpp_out

#results should be the same
all.equal(prcomp_out$sdev[1:3], rcpp_out$sdev)
all.equal(prcomp_out$rotation, rcpp_out$loadings,
          check.attributes = FALSE)
all.equal(prcomp_out$center, rcpp_out$center)

#runtime of rcpp vs prcomp
runtime = microbenchmark::microbenchmark(rcpp_out,prcomp_out)
autoplot(runtime)

#VS sparsepca
sparsepca_out = sparsepca::spca(X, k = 3, alpha = 0, beta = 0, verbose = F)
sparsepca_out

all.equal(sparsepca_out$sdev, rcpp_out$sdev)
all.equal(sparsepca_out$loadings, rcpp_out$loadings)
all.equal(sparsepca_out$center, rcpp_out$center)

runtime2 = microbenchmark::microbenchmark(rcpp_out,sparsepca_out)
autoplot(runtime2)




