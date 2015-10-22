#' @title Topological Anomaly Detection
#' @description Unsupervised anomaly detection that targets clusters of
#'   anomalies
#' @details Given n observations, consider an n-by-n distance matrix M. Let d
#'   denote the rq-th percentile of distances in M1. Set all distances greater
#'   than d to 0 and all other entries to 1. Call this new matrix M2. Define a
#'   graph g whose adjacency matrix is given by M2. Identify all connected
#'   components of g and count the number of observations in each component. Any
#'   component containing fewer than p*n observations is flagged as anomalous,
#'   and all observations in that component are determined to be anomalies.
#'
#'   For each observation flagged as an anomaly, determine the shortest distance
#'   from that observation to a non-anomalous observation. This distance is the
#'   anomaly score for that observation.
#'
#' @param dat Input data, as a matrix or dataframe. Each row corresponds to a
#'   unique observation, columns correspond to variables.
#'
#' @param rq If not using the KNN method, this is the "graph resolution" as a
#'   quantile of the observed inter-observation distances. Defaults to .1
#'   (10\%).
#'
#' @param p The "background" threshold as a percentile of the number of
#'   observations in the input data. Defaults to .1, i.e. observations in
#'   connected components containing 10\% or fewer of the dataset will be
#'   flagged as outliers.
#'
#' @param knn Boolean, whether or not to represent the data as a nearest-
#'   neighbor graph (as opposed to a thresholded distance graph). Defaults to
#'   FALSE, i.e. not using KNN method.
#'
#' @param k If using KNN method, determines the number of neighbors considered
#'   in constructing the nearest neighbors graph. Defaults to 3.
#'
#' @param return_scores Boolean value that determines whether or not outleir
#'   scores will be calculated and returned. If scores are not needed, setting
#'   this parameter to FALSE will impart speed gains on the classification.
#'   Defaults to TRUE.
#'
#' @param return_graph Boolean value that determines whether or not to return
#'   the graph representation of the data (as an igraph::graph object). Defaults
#'   to FALSE.
#'
#' @param return_dist Boolean value that determines whether or not to return the
#'   inter-observation distance matrix (as a dist object). If knn=FALSE and
#'   return_scores=FALSE this will require an additional separate calculation.
#'   Defaults to FALSE.
#'
#' @param verbose Boolean that determiens whether or not to run in "verbose"
#'   mode. Defaults to FALSE.
#'
#' @param rqSample This is an optional performance tweak for large datasets. If non-null
#'   and knn=FALSE, determines the max size of the sample of inter-observation
#'   distances to be used to determine graph resolution from the rq quantile
#'   parameter. If NULL, all distances are used. Defaults to NULL. If rqSample < n(n-1)/2,
#'   setting this parameter does not modify the behavior of the function.
#'
#' @param method Determines the distance metric (to be passed to proxy::dist) to be used to calculate inter-observation
#'   distances. For a list of all available techniques, run: names(proxy::pr_DB$get_entries())
#'
#' @return Numeric vector giving number of characters in each element of the
#'   character vector.  Missing strings have missing length.
#'
#' @import igraph
#' @import proxy
#' @import cccd
#' @export
#' @examples
#' library(igraph)
#' data(iris)
#' pca = prcomp(iris[,-5] ,center=TRUE, scale=TRUE)
#' colpal=colorRampPalette(c("yellow","red"))(10)
#'
#' #'  ## Normal algorithm ##
#' tad = TAD(iris[,-5], return_graph=TRUE, verbose=TRUE, return_scores=TRUE)
#' g = tad$graph
#'
#' # TAD-PCA plot
#' outl_discrete = ceiling(10*tad$scores/max(tad$scores))
#' V(g)$color=4
#' V(g)$color[tad$outliers] = colpal[outl_discrete]
#' V(g)$name[-tad$outliers] = ""
#' #V(g)$name[tad$outliers] = round(tad$scores,2)
#' par(mfrow=c(1,1))
#' plot(g, layout=pca$x[,1:2])
#' pairs(pca$x, col=V(g)$color, pch=16)
#'
#' ## KNN variant ##
#' tad = TAD(iris[,-5], knn=TRUE, k=5, p=.05, return_graph=TRUE, verbose=TRUE, rqSample=2250, return_dist=TRUE, return_scores=TRUE)
#' g = tad$graph
#'
#' # TAD-PCA plot
#' outl_discrete = ceiling(10*tad$scores/max(tad$scores))
#' V(g)$color=3
#' V(g)$color[tad$outliers] = colpal[outl_discrete]
#' par(mfrow=c(1,1))
#' plot(g, layout=pca$x[,1:2])
#' pairs(pca$x, col=V(g)$color, pch=16)

#library(igraph)
#library(cccd) # for nearest-neighbor variant
TAD = function(dat,
               rq = .1,
               p = .1,
               knn = FALSE,
               k=3,
               return_scores = TRUE,
               return_graph  = FALSE,
               return_dist   = FALSE,
               verbose=FALSE,
               rqSample=NULL,
               method='euclidean'
){
  start = Sys.time()
  n=nrow(dat)

  if(!knn | return_scores | return_dist){
    if(verbose) {
      print(Sys.time()-start)
      start = Sys.time()
      print (c("n:",n))
      print ( "Calculating distance matrix..." )
    }
    d = dist(dat[,-5], method=method)
  }

  if(knn){
    g = tryCatch(
      nng(dat, k=k, use.fnn=TRUE),
      error=function(e){
        nng(dat, k=k)
      })
  } else {
    if(verbose) {
      print(Sys.time()-start)
      start = Sys.time()
      print ("Calculating graph resolution (from distance quantile)...")
    }
    if( is.null(rqSample) ){
      r = quantile(d, rq)
    } else {
      if(length(d) <= rqSample ){
        r = quantile(d, rq)
      } else {
        r = quantile(sample(d, rqSample, replace=FALSE), rq)
      }
    }

    if(verbose){
      print(Sys.time()-start)
      start = Sys.time()
      print (c("Graph resolution: r=",r))
      print ( "Extracting edges..." )
    }

    edges = t(sapply(which(d<=r), function(x) row_col_from_condensed_index(n,x)))

    if(verbose) {
      print(Sys.time()-start)
      start = Sys.time()
      print ( "Building graph..." )
    }

    # This method ensures that node index corresponds with node name
    g = igraph::graph.empty(n, directed=FALSE)
    g[from=edges[,1], to=edges[,2]]=1
  }
  igraph::V(g)$name = 1:n

  if(verbose) {
    print(Sys.time()-start)
    start = Sys.time()
    print ( "Finding connected components..." )
  }

  #components = clusters(g, mode="weak")
  components = clusters(g, mode="strong")
  csize = components$csize[components$membership]
  if(verbose) {
    print(Sys.time()-start)
    #start = Sys.time()
    print ( "Flagging outliers..." )
  }

  outliers = V(g)$name[csize<=p*n]
  retval = list(outliers=outliers)

  if(return_graph){
    V(g)$outlier = csize<=p*n
    retval$graph = g}
  if(return_dist) retval$dist = d
  if(return_scores & !is.null(outliers)){
    inliers  = c(1:n)[-outliers]
    m=length(outliers)
    scores = rep(NA, m)
    for(i in 1:(m)){
      ix = condensed_index_from_row_col(outliers[i], inliers, n) # i,j,n
      scores[i] = min(d[ix])
    }
    retval$scores = scores
    V(g)$score = NA
    V(g)$score[outliers] = scores
  }
  retval
}

