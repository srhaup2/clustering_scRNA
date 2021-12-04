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
#'@return A list
#'\itemize{
#'   \item clusters - n x p+1 matrix of cluster assignments where the first column is cluster assignments
#'   \item iter - number of iterations
#'   \item centroids - k+1-length vector of centroids where the first column is cluster assignments
#'   \item wcsse - min within-cluster SSE over all nstart iterations
#'}
#'
#'@examples
#'
#'
#'@export
#'

kmeans_clust = function(X, k, nstart){
  ### Setup ###

  # get matrix dims
  n = nrow(X)
  p = ncol(X)

  # normalize data (helps with clustering)
  X = scale(X, center = TRUE, scale = apply(X,2,sd))


  ### Repeat k-means method nstart times ###
  for (i in 1:nstart){
    ### Initial clusters ###

    # Method 1 - random assignment (easier but worse clusters)
    rand = sample(1:n, k, replace = FALSE)
    centroids_i = cbind(seq_along(rand),X[rand, ])
    cluster_vec = rep(0,n)
    cluster_vec[rand] = seq_along(rand)
    clusters_i = cbind(cluster_vec, X) #column 1 is cluster ID

    # Method 2 - initialize using kmeans++ (better clusters)



    ### K-means algorithm ###

    # Method 1 - LLoyd's method (faster and simpler)
    iter_i = 0
    conv = 0

    while(conv == 0) {
      # increment
      iter_i = iter_i + 1

      # calculate distance
      dists = as.matrix(dist(rbind(clusters_i[,-1], centroids_i[,-1])))[1:n,(n+1):(n+k)]^2 # n x k squared distance matrix
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
    dists2 = as.matrix(dist(rbind(clusters_i[,-1], centroids_i[,-1])))[1:n,(n+1):(n+k)]^2 # n x k squared distance matrix
    pt_to_cent = diag(dists2[,clusters_i[,1]])
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

    # Method 2 - Hartigan-Wong method (better clusters)

  }

  return(list(clusters = clusters, iter = iter, centroids = centroids, wcsse = wcsse))
}
