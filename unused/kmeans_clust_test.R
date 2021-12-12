library(tidyverse) #data manipulation
library(mclust) #adjusted Rand Index
library(microbenchmark) #speed comparison
#library(clusteringscRNA) #my k means function

### Brain data

load("data/tumor_reduced.rda")
load("data/TC.rda")

tumor_reduced = tumor_reduced[,1:500]

# Speed

microbenchmark(
  kmeans(tumor_reduced, centers = 5, nstart = 1, algorithm = "Lloyd"),
  kmeans_clust(tumor_reduced, k = 5, nstart = 1, init.method = "random"),
  kmeans_clust(tumor_reduced, k = 5, nstart = 1, init.method = "kmeans++"),
  kmeans_clust(tumor_reduced, k = 5, nstart = 1, init.method = "gkmeans++"),
  times = 10
)

# Accuracy (takes a long time to run!)
res1 = numeric(10)
res2 = numeric(10)
res3 = numeric(10)
res4 = numeric(10)

for (i in 1:10) {
  res1[i] = adjustedRandIndex(kmeans(tumor_reduced, centers = 5, nstart = 10, algorithm = "Lloyd")$cluster,
                    TC)
}

min(res1)
max(res1)
median(res1)
mean(res1)
sd(res1)



for (i in 1:10) {
  res2[i] = adjustedRandIndex(kmeans_clust(tumor_reduced, k = 5, nstart = 10, init.method = "random")$clusters[,1],
                              TC)
}

min(res2)
max(res2)
median(res2)
mean(res2)
sd(res2)


for (i in 1:10) {
  res3[i] = adjustedRandIndex(kmeans_clust(tumor_reduced, k = 5, nstart = 10, init.method = "kmeans++")$clusters[,1],
                              TC)
}

min(res3)
max(res3)
median(res3)
mean(res3)
sd(res3)


for (i in 1:10) {
  res4[i] = adjustedRandIndex(kmeans_clust(tumor_reduced, k = 5, nstart = 10, init.method = "gkmeans++")$clusters[,1],
                              TC)
}

min(res4)
max(res4)
median(res4)
mean(res4)
sd(res4)




