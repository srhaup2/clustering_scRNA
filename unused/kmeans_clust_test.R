library(tidyverse) #data manipulation
library(mclust) #adjusted Rand Index
library(microbenchmark) #speed comparison
library(clusteringscRNA) # our package

### Brain data

# if you want to make it smaller
tumor_reduced2 = tumor_reduced[,1:500]

# Speed
microbenchmark(
  kmeans(tumor_reduced2, centers = 5, nstart = 1, algorithm = "Lloyd"),
  kmeans_clust(tumor_reduced2, k = 5, nstart = 1, init.method = "random"),
  kmeans_clust(tumor_reduced2, k = 5, nstart = 1, init.method = "kmeans++"),
  kmeans_clust(tumor_reduced2, k = 5, nstart = 1, init.method = "gkmeans++"),
  times = 10
)

# Accuracy (takes a long time to run!)
res1 = matrix(0,2,10)
rownames(res1) = c("ARI", "WCSSE")
res2 = res1
res3 = res1
res4 = res1

for (i in 1:10) {
  m = kmeans(tumor_reduced2, centers = 5, nstart = 5, iter.max = 25, algorithm = "Lloyd")
  res1[1,i] = adjustedRandIndex(m$cluster,TC)
  res1[2,i] = m$tot.withinss
}

print("ARI")
min(res1[1,])
median(res1[1,])
mean(res1[1,])
max(res1[1,])

print("WCSSE")
min(res1[2,])
median(res1[2,])
mean(res1[2,])
max(res1[2,])




for (i in 1:10) {
  m = kmeans_clust(tumor_reduced2, k = 5, nstart = 5, iter.max = 25, init.method = "random")
  res2[1,i] = adjustedRandIndex(m$clusters[,1],TC)
  res2[2,i] = m$wcsse
}

print("ARI")
min(res2[1,])
median(res2[1,])
mean(res2[1,])
max(res2[1,])

print("WCSSE")
min(res2[2,])
median(res2[2,])
mean(res2[2,])
max(res2[2,])



for (i in 1:10) {
  m = kmeans_clust(tumor_reduced2, k = 5, nstart = 5, iter.max = 25, init.method = "kmeans++")
  res3[1,i] = adjustedRandIndex(m$clusters[,1],TC)
  res3[2,i] = m$wcsse
}

print("ARI")
min(res3[1,])
median(res3[1,])
mean(res3[1,])
max(res3[1,])

print("WCSSE")
min(res3[2,])
median(res3[2,])
mean(res3[2,])
max(res3[2,])


for (i in 1:10) {
  m = kmeans_clust(tumor_reduced2, k = 5, nstart = 5, iter.max = 25, init.method = "gkmeans++")
  res4[1,i] = adjustedRandIndex(m$clusters[,1],TC)
  res4[2,i] = m$wcsse
}

print("ARI")
min(res4[1,])
median(res4[1,])
mean(res4[1,])
max(res4[1,])

print("WCSSE")
min(res4[2,])
median(res4[2,])
mean(res4[2,])
max(res4[2,])





