library(tidyverse) #data manipulation
library(mclust) #adjusted Rand Index
library(microbenchmark) #speed comparison
#library(clusteringscRNA) #my k means function

### Brain data

brain = read_csv("data/data_500.csv") %>%
  select(-1) %>%
  rename(sample = X)
brain_labs = read_csv("data/labels.csv") %>%
  rename(sample = "...1")

# Speed

microbenchmark(
  kmeans(brain[,-1], centers = 5, nstart = 1, algorithm = "Lloyd"),
  kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "random"),
  kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "kmeans++"),
  kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "gkmeans++"),
  times = 10
)

# Accuracy (takes a long time to run!)
res1 = numeric(10)
res2 = numeric(10)
res3 = numeric(10)
res4 = numeric(10)

for (i in 1:10) {
  res1[i] = adjustedRandIndex(kmeans(brain[,-1], centers = 5, nstart = 1, algorithm = "Lloyd")$cluster,
                    brain_labs$Class)
}

# > min(res1)
# [1] 0.5809174
# > max(res1)
# [1] 0.8136526
# > median(res1)
# [1] 0.7699028
# > mean(res1)
# [1] 0.7328307
# > sd(res1)
# [1] 0.08393852


for (i in 1:10) {
  res2[i] = adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "random")$clusters[,1],
                              brain_labs$Class)
}

# > min(res2)
# [1] 0.4614844
# > max(res2)
# [1] 0.8963029
# > median(res2)
# [1] 0.7707366
# > mean(res2)
# [1] 0.7439101
# > sd(res2)
# [1] 0.1328754


for (i in 1:10) {
  res3[i] = adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "kmeans++")$clusters[,1],
                              brain_labs$Class)
}

# > min(res3)
# [1] 0.4339616
# > max(res3)
# [1] 0.9045912
# > median(res3)
# [1] 0.7619758
# > mean(res3)
# [1] 0.6882781
# > sd(res3)
# [1] 0.1733534


for (i in 1:10) {
  res4[i] = adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 1, init.method = "gkmeans++")$clusters[,1],
                              brain_labs$Class)
}

# > min(res4)
# [1] 0.5606479
# > max(res4)
# [1] 0.8963029
# > median(res4)
# [1] 0.8405325
# > mean(res4)
# [1] 0.823571
# > sd(res4)
# [1] 0.09931881




