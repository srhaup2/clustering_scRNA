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

# Accuracy

adjustedRandIndex(kmeans(brain[,-1], centers = 5, nstart = 10, algorithm = "Lloyd")$cluster,
                  brain_labs$Class)


adjustedRandIndex(kmeans_clust(brain[,-1], k = 5, nstart = 10)$clusters[,1],
                  brain_labs$Class)


### Mouse data
