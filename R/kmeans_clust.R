#'kmeans_clust
#'
#' Clusters data according to k-means algorithm
#'
#'@param X n x p data matrix
#'
#'@param k number of clusters
#'
#'@param nstart number of times to perform k means on the data
#'
#'@param iter.max max iterations per run
#'
#'@param init.method method for centroid initialization - choose from random, kmeans++, greedy kmeans++ (gkmeans++)
#'
#'@return A list
#'\itemize{
#'   \item clusters - n x p+1 matrix of cluster assignments where the first column is cluster assignments
#'   \item iter - number of iterations
#'   \item centroids - k x p+1 matrix of centroids where the first column is cluster assignments
#'   \item wcsse - min within-cluster SSE over all nstart iterations
#'}
#'
#'@examples
#'kmeans_clust(tumor_reduced, k = 3, nstart = 5, init.method = "kmeans++")
#'
#'@import stats
#'@import Rcpp
#'@export
#'

kmeans_clust = function(X, k, nstart = 1L, iter.max = 10L, init.method = "random"){
  ### Setup ###
  
  # get matrix dims
  n = nrow(X)
  p = ncol(X)
  
  # normalize data (helps with clustering)
  X = scale(X, center = TRUE, scale = apply(X,2,sd))
  X[is.nan(X)] = 0
  
  
  ### Repeat k-means method nstart times ###
  for (i in 1:nstart){
    
    ### Initial clusters ###
    if (init.method == "random"){
      # Method 1 - random assignment (easier but worse clusters)
      
      # random centroids (by row number)
      rand = sample(1:n, k, replace = FALSE)
      centroids_i = cbind(seq_along(rand),X[rand, ])
      # cluster assignments for centroids
      cluster_vec = rep(0,n)
      cluster_vec[rand] = seq_along(rand)
      clusters_i = cbind(cluster_vec, X) #column 1 is cluster ID
    } else if (init.method == "kmeans++") {
      # Method 2 - initialize using kmeans++ (better clusters)
      
      #random centroid
      rand = numeric(k)
      centroids_i = matrix(0,k,p)
      open_pos = 1:n
      rand[1] = sample(open_pos, 1)
      open_pos = open_pos[-rand[1]]
      centroids_i[1,] = X[rand[1], ]
      for (j in 2:k) {
        #compute distances, get min for each point
        if (j == 2){ # only one centroid
          dists = dist_squared(X, t(as.matrix(centroids_i[1:j-1,])))
          min_dists = dists[open_pos]
          prob_vec = min_dists/sum(min_dists)
        } else { #many centroids
          dists = dist_squared(X, centroids_i[1:j-1,])
          min_dists = apply(dists,1,which.min)[open_pos]
          prob_vec = min_dists/sum(min_dists)
        }
        # select next centroid with prob proportional to distance
        rand[j] = sample(open_pos, 1, prob = prob_vec)
        open_pos = open_pos[-rand[j]]
        centroids_i[j,] = X[rand[j], ]
      }
      # cluster assignments for centroids
      centroids_i = cbind(seq_along(rand), centroids_i)
      cluster_vec = rep(0,n)
      cluster_vec[rand] = seq_along(rand)
      clusters_i = cbind(cluster_vec, X) #column 1 is cluster ID
      
    } else if (init.method == "gkmeans++"){
      # Method 3 - greedy kmeans++ (even better clusters)
      
      #random centroid
      rand = numeric(k)
      centroids_i = matrix(0,k,p)
      open_pos = 1:n
      rand[1] = sample(open_pos, 1)
      open_pos = open_pos[-rand[1]]
      centroids_i[1,] = X[rand[1], ]
      for (j in 2:k) {
        #compute distances, get min for each point
        if (j == 2){ # only one centroid
          dists = dist_squared(X, t(as.matrix(centroids_i[1:j-1,])))
          min_dists = dists[open_pos]
          prob_vec = min_dists/sum(min_dists)
        } else { # many centroids
          dists = dist_squared(X, centroids_i[1:j-1,])
          min_dists = apply(dists,1,which.min)[open_pos]
          prob_vec = min_dists/sum(min_dists)
        }
        # sample l = 2 + log(k) new centroids (this is default in scikitlearn implementation)
        samp = sample(open_pos, ifelse(floor(2 + log(k)) > length(open_pos), length(open_pos),floor(2 + log(k))), prob = prob_vec)
        tmp_wcsse = numeric(ifelse(floor(2 + log(k)) > length(open_pos), length(open_pos),floor(2 + log(k))))
        for (l in (1:length(samp))) {
          dists = dist_squared(X,rbind(centroids_i[1:j-1,],X[samp[l],]))
          min_dists = apply(dists,1,which.min) # closest centroid to each point
          clusters_i = cbind(min_dists, X)
          pt_to_cent = diag(dists[,clusters_i[,1]])
          tmp_wcsse[l] = sum(pt_to_cent)
        }
        
        # select next centroid based on lowest wcsse
        rand[j] = samp[which.min(tmp_wcsse)]
        open_pos = open_pos[-rand[j]]
        centroids_i[j,] = X[rand[j], ]
        
      }
      # cluster assignments for centroids
      centroids_i = cbind(seq_along(rand), centroids_i)
      cluster_vec = rep(0,n)
      cluster_vec[rand] = seq_along(rand)
      clusters_i = cbind(cluster_vec, X) #column 1 is cluster ID
      
    }
    
    ### K-means algorithm ###
    
    # LLoyd's method (faster and simpler)
    iter_i = 0
    conv = 0
    
    while(conv == 0 && iter_i <= iter.max) {
      # increment
      iter_i = iter_i + 1
      
      # calculate distance
      dists = dist_squared(clusters_i[,-1], centroids_i[,-1]) # n x k squared distance matrix
      min_dists = apply(dists,1,which.min) # closest centroid to each point
      
      # re-assign clusters
      clusters_before = clusters_i[,1]
      clusters_i[,1] = min_dists
      
      # calculate new centroids
      clust_sums = rowsum(clusters_i[,-1],clusters_i[,1]) # colsums by group
      clust_size = as.matrix(aggregate(as.data.frame(clusters_i[,2]), as.data.frame(clusters_i[,1]), length))[,2]
      centroids_i = cbind(centroids_i[,1],clust_sums/clust_size)
      
      
      # if assignments didnt change, we reached convergence
      if (all(clusters_i[,1] == clusters_before)) {
        conv = 1
      }
    }
    
    # calc wcsse for this run of k-means
    pt_to_cent = diag(dists[,clusters_i[,1]])
    wcsse_i = sum(pt_to_cent)
    
    # if this is our best clustering so far, reassign return values
    if (i == 1) {
      wcsse = wcsse_i
      clusters = clusters_i
      iter = iter_i
      centroids = centroids_i
    } else if (wcsse_i < wcsse){
      wcsse = wcsse_i
      clusters = clusters_i
      iter = iter_i
      centroids = centroids_i
    }
    
  }
  
  return(list(clusters = clusters, iter = iter, centroids = centroids, wcsse = wcsse))
}
