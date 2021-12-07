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
library(knitr)
library(tidyverse) #data manipulation
library(sparsepca)   #original spca package
#devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp)

#--------------------------------------------------------------------
#   Prepare Dataset & Subsampling
#--------------------------------------------------------------------
# tumor = read.csv("data/data.csv")
# tumor_reduced = tumor[,2:2001]
# write.csv(tumor_reduced,"data/tumor_reduced.csv", row.names = FALSE)

tumor_reduced = read.csv("data/tumor_reduced.csv")
tumor_var = apply(tumor_reduced, 2, var)
tumor_var2 = tumor_reduced[ , order(tumor_var, decreasing = T) ]  
# tumor2 = tumor_var2[, 1:200]
tumor2 = tumor_reduced[1:200]
tumor3 = tumor_reduced[,1:500]

####tumor2 = tumor[, 3:ncol(tumor)]
TC <- read.csv('data/labels.csv', row.names = NULL)[,2]
nclust = 5
n = length(TC)

spca_out = spcaRcpp(tumor2, k = 21, alpha = 1e-4, beta = 1e-4)
spca_out2 = spcaRcpp(tumor2, k = 21, alpha = 0, beta = 0)
spca_out3 = spcaRcpp(tumor3, k = 21, alpha = 1e-4, beta = 1e-4)
prcomp_out = prcomp(tumor2, rank.=21)
sparse_out = sparsepca::spca(tumor2, k = 21, alpha = 1e-4, beta = 1e-4, center = T, scale =T, verbose = F)
sparse_out2 = sparsepca::spca(tumor2, k = 21, alpha = 0, beta = 0, verbose = F)

#calculate culmulated explained var ratio
explained_var <- function(obj){
  variance = obj$sdev ^ 2
  explained_var_ratio = variance / obj$var
  cum_explained_ratio = cumsum(explained_var_ratio)
  return (cum_explained_ratio)
}
explained_var(spca_out)

#similar to autoplot(prcomp_out)
plot_cluster<-function(obj){
  dtp <- data.frame("tumors" = TC, obj$scores[, 1:2])
  ggplot(data = dtp) +
    geom_point(aes (x = X1, y = X2, col = tumors))+
    xlab("PC1") + ylab("PC2") +
    ggtitle ("Cluster of Tumor scRNA-seq") + 
    theme_minimal()
}
plot_cluster(spca_out2)

normMixEm_test <- function(data, num_components= 5L ){
  EM <- normMixEm$new(input_dat = data,num_components = num_components)
  res_EM <- EM$run.EM(loglik_tol=1e-5)
}

#--------------------------------------------------------------------
#   var selection with EM
#--------------------------------------------------------------------
res_EM0 = normMixEm_test(data = as.matrix(tumor2), num_components= 5L)
class0 <- apply(res_EM0$prob_mat, 1, which.max)
kable(data.frame(est=class0, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM0$iter,res_EM0$loglik_list,type="l")

#--------------------------------------------------------------------
#   different choices of PCA on EM
#--------------------------------------------------------------------
#1. spcaRcpp (alpha = beta = 1e-4)
set.seed(202112)
s = Sys.time()
res_EM = normMixEm_test(data = spca_out$scores, num_components= 5L)
Sys.time() - s
class <- apply(res_EM$prob_mat, 1, which.max)
kable(data.frame(est=class, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))

table(class,TC)
#plot(1L:res_EM$iter,res_EM$loglik_list,type="l")

#perm_class = class  #STORE SEED! DO NOT MODIFY
kable(data.frame(est=perm_class, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
adjustedRandIndex(perm_class, TC)

#2. spcaRcpp (alpha = beta = 0)
res_EM2 = normMixEm_test(data = spca_out2$scores, num_components= 5L)
class2 <- apply(res_EM2$prob_mat, 1, which.max)
kable(data.frame(est=class2, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))
#plot(1L:res_EM2$iter,res_EM2$loglik_list,type="l")

#3. prcomp
res_EM3 = normMixEm_test(data = as.matrix(prcomp_out$x), num_components= 5L)
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
mclust::adjustedRandIndex(TC, perm_class)    #spcaRcpp (alpha = beta = 1e-4)
mclust::adjustedRandIndex(TC, class2)   #spcaRcpp (alpha = beta = 0)
mclust::adjustedRandIndex(TC, class3)   #prcomp
mclust::adjustedRandIndex(TC, class4)   #sparsepca (alpha = beta = 1e-4)
mclust::adjustedRandIndex(TC, class5)   #sparsepca (alpha = beta = 0)

spcaEM <- function(data){
  res_EM = normMixEm_test(data = data$scores, num_components= 5L)
  class <- apply(res_EM$prob_mat, 1, which.max)
  return (class)
}

full_class = spcaEM(full_out)
ARI_EM = numeric(10)
for (i in 1:10) {
  ARI_EM[i] = adjustedRandIndex(spcaEM (spca_out),
                              TC)
}

#--------------------------------------------------------------------
#   Speed Test spca + mixtools
#--------------------------------------------------------------------
s = Sys.time()
mix_out = mixtools::mvnormalmixEM (spca_out$scores, arbvar = F,k=5)
Sys.time() - s

class_mix <- apply(mix_out$posterior, 1, which.max)
table(class_mix,TC)
kable(data.frame(est=class_mix, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))



#--------------------------------------------------------------------
#   Speed Test for spcaRcpp
#--------------------------------------------------------------------

s = Sys.time()
prcomp_out = prcomp(tumor2, rank. = 21)
#prcomp_out
Sys.time() - s

s = Sys.time()
rcpp_out = spcaRcpp(tumor2, k = 21, alpha = 1e-4, beta = 1e-4)
#rcpp_out
Sys.time() - s

s = Sys.time()
sparse_out = sparsepca::spca(tumor2, k = 21, alpha = 1e-4, beta = 1e-4, verbose = F)
#sparse_out
Sys.time() - s


#--------------------------------------------------------------------
#   Microbenchmark spcaRcpp vs sparsepca
#--------------------------------------------------------------------
sparsepca <- function(data) return (sparsepca::spca(data, k = 21, alpha = 1e-4, beta = 1e-4,verbose = F))
rcpp <- function(data) return (spcaRcpp(data, k = 21, alpha = 1e-4, beta = 1e-4))

mb = microbenchmark(sparsepca(tumor2),
                    rcpp(tumor2),
                    times = 100)
autoplot(mb)

#--------------------------------------------------------------------
#   spcaRcpp + kmeans
#--------------------------------------------------------------------
set.seed(20211205)
s = Sys.time()
#spca_out = spcaRcpp(tumor2, k = 21, alpha = 1e-4, beta = 1e-4)
res_kmeans = kmeans_clust(spca_out$scores, k =5, init.method = "gkmeans++")
Sys.time() - s
PC <- res_kmeans$clusters[,1]
kable(data.frame(est=PC, true=TC) %>% count(true,est) %>% mutate(freq=n/sum(n)))

table(PC,TC)

spcaKmeans <- function(scores){
  res_kmeans = kmeans_clust(scores, k =5, init.method = "gkmeans++")
  class <- res_kmeans$clusters[,1]
  return (class)
}

ARI_kmeans = numeric(10)
for (i in 1:10) {
  ARI_kmeans[i] = adjustedRandIndex(spcaKmeans(spca_out$scores),
                              TC)
}

res_kmeans2 = numeric(10)
for (i in 1:10) {
  res_kmeans2[i] = adjustedRandIndex(spcaKmeans(tumor2),
                                    TC)
}

s = Sys.time()
res_kmeans = kmeans_clust(spca_out3$scores, k =5, init.method = "gkmeans++")
Sys.time() - s


#--------------------------------------------------------------------
#   microbenchmark 
#--------------------------------------------------------------------
# mb_km_em = microbenchmark(spcaEM(spca_out),
#                           spcaKmeans(spca_out),
#                           times = 50)


#suicide attempt
s = Sys.time()
full_out = sparsepca::spca(tumor_reduced, k = 100, verbose = F)
Sys.time() - s



