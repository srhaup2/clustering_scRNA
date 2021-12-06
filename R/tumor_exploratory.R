list.of.packages <- c("devtools", "microbenchmark", "ggplot2", "ggfortify",
                      "mclust", "sparsepca", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(devtools)
library(microbenchmark)
library(ggplot2)
library(ggfortify)
library(mclust)
library(dplyr)
library(sparsepca)   #original spca package
#devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp) 

set.seed(322)
tumor = read.csv("data/data_500.csv")
tumor2 = tumor[, 3:ncol(tumor)]
TC <- read.csv('data/labels.csv', row.names = NULL)[,2]
nclust = 5
n = length(TC)

spca_out = spcaRcpp(tumor2, k = 20, alpha = 1e-4, beta = 1e-4)
spca_out2 = spcaRcpp(tumor2, k = 20, alpha = 0, beta = 0)
prcomp_out = prcomp(tumor2, rank.=20)
sparse_out = sparsepca::spca(tumor2, k = 20, alpha = 1e-4, beta = 1e-4, verbose = F)
sparse_out2 = sparsepca::spca(tumor2, k = 20, alpha = 0, beta = 0, verbose = F)

#calculate cumulated explained var ratio
variance = spca_out$sdev ^ 2
explained_var_ratio = variance / spca_out$var
cum_explained_ratio = cumsum(explained_var_ratio)
cum_explained_ratio

#similar to autoplot(prcomp_out)
dtp <- data.frame("tumors" = TC, spca_out$scores[, 1:2])
ggplot(data = dtp) +
  geom_point(aes (x = X1, y = X2, col = tumors))+
  theme_minimal()

normMixEm_test <- function(data, num_components= 5L ){
  EM <- normMixEm$new(input_dat = data,num_components = num_components)
  res_EM <- EM$run.EM(loglik_tol=1e-5)
}

#--------------------------------------------------------------------
#   different choices of PCA on EM
#--------------------------------------------------------------------

#1. spcaRcpp (alpha = beta = 1e-4)
res_EM = normMixEm_test(data = spca_out$scores, num_components= 5L)
class <- apply(res_EM$prob_mat, 1, which.max)
print(data.frame(est=class, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM$iter,res_EM$loglik_list,type="l")

#2. spcaRcpp (alpha = beta = 0)
res_EM2 = normMixEm_test(data = spca_out2$scores, num_components= 5L)
class2 <- apply(res_EM2$prob_mat, 1, which.max)
print(data.frame(est=class2, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM2$iter,res_EM2$loglik_list,type="l")

#3. prcomp
res_EM3 = normMixEm_test(data = prcomp_out$x, num_components= 5L)
class3 <- apply(res_EM3$prob_mat, 1, which.max)
print(data.frame(est=class3, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM3$iter,res_EM3$loglik_list,type="l")

#4. sparsepca::pca (alpha = beta = 1e-4)
res_EM4 = normMixEm_test(data = as.matrix(sparse_out$scores), num_components= 5L)
class4 <- apply(res_EM4$prob_mat, 1, which.max)
print(data.frame(est=class4, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM4$iter,res_EM4$loglik_list,type="l")

#5. sparsepca::pca (alpha = beta = 0)
res_EM5 = normMixEm_test(data = as.matrix(sparse_out2$scores), num_components= 5L)
class5 <- apply(res_EM5$prob_mat, 1, which.max)
print(data.frame(est=class5, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM5$iter,res_EM5$loglik_list,type="l")

#--------------------------------------------------------------------
#   Adjusted Rand Index
#--------------------------------------------------------------------
mclust::adjustedRandIndex(TC, class)    #spcaRcpp (alpha = beta = 1e-4)
mclust::adjustedRandIndex(TC, class2)   #spcaRcpp (alpha = beta = 0)
mclust::adjustedRandIndex(TC, class3)   #prcomp
mclust::adjustedRandIndex(TC, class4)   #sparsepca (alpha = beta = 1e-4)
mclust::adjustedRandIndex(TC, class5)   #sparsepca (alpha = beta = 0)




#--------------------------------------------------------------------
#   Speed Test for spcaRcpp
#--------------------------------------------------------------------

s = Sys.time()
prcomp_out = prcomp(tumor2, rank. = 14)
#prcomp_out
Sys.time() - s

s = Sys.time()
rcpp_out = spcaRcpp(tumor2, k = 14, alpha = 1e-4, beta = 1e-4)
#rcpp_out
Sys.time() - s

s = Sys.time()
sparse_out = sparsepca::spca(tumor2, k = 14, alpha = 1e-4, beta = 1e-4, verbose = F)
#sparse_out
Sys.time() - s

# mb = microbenchmark(#prcomp(tumor2, rank. = 14),
#                     spcaRcpp(tumor2, k = 14, alpha = 0, beta = 0),
#                     sparsepca::spca(tumor2, k = 14, alpha = 0, beta = 0,verbose = F),
#                     times = 100)
# autoplot(mb)





