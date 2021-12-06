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
  kmeans_clust(brain[,-1], k = 5, nstart = 1)
)

# Accuracy (takes a long time to run!)
res1 = numeric(10)
res2 = numeric(10)
res3 = numeric(10)

for (i in 1:10) {
  res1[i] = adjustedRandIndex(kmeans(brain[,-1], centers = 5, nstart = 5, algorithm = "Lloyd")$cluster,
                    brain_labs$Class)
}

for (i in 1:10) {
  res2[i] = adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 5, init.method = "random")$clusters[,1],
                              brain_labs$Class)
}

for (i in 1:10) {
  res3[i] = adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 5, init.method = "kmeans++")$clusters[,1],
                              brain_labs$Class)
}

# results




