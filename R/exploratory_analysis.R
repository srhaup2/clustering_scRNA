library(devtools)
library(microbenchmark)
library(ggplot2)
library(sparsepca)
#devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp)

rna = read.csv("data/nestorowa_corrected_log2_transformed_counts.txt", sep="")
dim(rna)
rna2 = rna[1:1000,1:300]

s = Sys.time()
prcomp_out = prcomp(rna, rank.=10)
#prcomp_out
Sys.time() - s

s = Sys.time()
spca_out = spcaRcpp(rna, k=10,alpha = 0, beta = 0)
#spca_out
Sys.time() - s

s = Sys.time()
ori_spca_out = sparsepca::spca(rna, k=10,alpha = 0, beta = 0)
#ori_spca_out
Sys.time() - s


res = microbenchmark(ori_spca_out, spca_out)
autoplot(res)

res2 = microbenchmark(prcomp_out, spca_out)
autoplot(res2)


s = Sys.time()
spca_out2 = spcaRcpp(rna2, k=10)
#spca_out
Sys.time() - s

s = Sys.time()
ori_spca_out2 = sparsepca::spca(rna2, k=10, verbose = F)
#ori_spca_out
Sys.time() - s

res3 = microbenchmark(spca_out2 = spcaRcpp(rna2, k=10),
                      ori_spca_out2 = sparsepca::spca(rna2, k=10,verbose=F))
autoplot(res3)

res4 = microbenchmark(spca_out3 = spcaRcpp(rna2, k=10),
                      prcomp_out3 = prcomp(rna2, rank. = 10))
autoplot(res4)





