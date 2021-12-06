library(devtools)
library(microbenchmark)
library(ggplot2)
library(sparsepca)
#devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp)

#Simulate 1000x10 Data
set.seed(322)
m <- 10000
V1 <- rnorm(m, -100,200)
V2 <- rnorm(m, -100, 300)
V3 <- -0.1*V1 + 0.1*V2 + rnorm(m, 0, 100)

X <- cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X <- X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))

center=T
spcaObj = list(loadings = NULL,
               eigenvalues = NULL,
               center = center)
#spcaObj$center <- colMeans(X)
spcaObj$center = t(replicate(nrow(X),colMeans(X)))

#X2 = sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE)
X = X - spcaObj$center

# all.equal(X2,X3)
#
# res = microbenchmark(sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE),
#                      X3 = X - center_m, times = 300)
#
# autoplot(res)

ori <- function(){
  spcaObj$center <- colMeans(X)
  X2 = sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE)
  return (X2)
}

test <- function(){
  spcaObj$center = t(replicate(nrow(X),colMeans(X)))
  X2 = X - spcaObj$center
  return (X2)
}

all.equal(ori(),test())
res = microbenchmark(ori(),test(),times=10)
autoplot(res)

k=3
svd_X <- svd(X)
Dmax  <- svd_X$d[1] # l2 norm

B <- svd_X$v[,1:k]
V <- svd_X$v

ori <- function(){
  VD  = sweep(V, MARGIN = 2, STATS = svd_X$d, FUN = "*")
  return (VD)
}

test <- function(){
  test_d = t(replicate(nrow(V),svd_X$d))
  V2 = V * test_d
  return (V2)
}

all.equal(ori(), test())
res2 = microbenchmark(ori(), test(), times = 300)
autoplot(res2)
