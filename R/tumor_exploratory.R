library(devtools)
library(microbenchmark)
library(ggplot2)
library(ggfortify)
library(sparsepca)   #original spca package
#devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp) 

set.seed(322)
tumor = read.csv("data/data.csv")
tumor2 = tumor[, 2:200]
TC <- read.csv('data/labels.csv', row.names = NULL)[,2]
nclust = 5
n = length(TC)

spca_out = spcaRcpp(tumor2, k = 14, alpha = 1e-4, beta = 1e-4)
spca_out2 = spcaRcpp(tumor2, k = 14, alpha = 0, beta = 0)
prcomp_out = prcomp(tumor2)

#calculate cumulated explained var ratio
variance = spca_out$sdev ^ 2
explained_var_ratio = variance / spca_out$var
cum_explained_ratio = cumsum(explained_var_ratio)

#similar to autoplot(prcomp_out)
dtp <- data.frame("tumors" = TC, spca_out$scores[, 1:2])
ggplot(data = dtp) +
  geom_point(aes (x = X1, y = X2, col = tumors))+
  theme_minimal()

normMixEm_test <- function(data, num_components= 5L ){
  EM <- normMixEm$new(input_dat = data,num_components = num_components)
  res_EM <- EM$run.EM(loglik_tol=1e-5)
}

#on spca data
res_EM = normMixEm_test(data = spca_out$scores, num_components= 5L)
class <- apply(res_EM$prob_mat, 1, which.max)
print(data.frame(est=class, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
plot(1L:res_EM$iter,res_EM$loglik_list,type="l")

#on pca data
res_EM2 = normMixEm_test(data = spca_out2$scores, num_components= 5L)
class2 <- apply(res_EM2$prob_mat, 1, which.max)
print(data.frame(est=class2, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
plot(1L:res_EM2$iter,res_EM2$loglik_list,type="l")








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





