#' Perform TCA based clustering
#'
#' @param esetRUV expression set RUV
#' @param genes gene names
#' @param clusterK number of clusters K
#'
#' @return list with all you need
#'
#' @importFrom EDASeq normCounts
#' @importFrom BoolNet binarizeTimeSeries
#'
#' @rdname additional-TCA
#'
#' @export
#'
clusterWithTca <- function(esetRUV, genes, clusterK=9) {
  normalizedCountsMatrix <- normCounts(esetRUV$setEmpirical)
  tca <- clusterTimes(normMatrix = normalizedCountsMatrix, days = as.character(esetRUV$setEmpirical$days),
                      genes = genes, clusterK = 9, takeLog = T)

  binaryCenter <- binarizeTimeSeries(tca@centers)
  clusters <- data.frame(gene=names(tca@cluster), cluster=tca@cluster)

  return(list(binaryCenter=binaryCenter,
              clusters=clusters,
              tca=tca))
}


#' Perform TCA based clustering
#'
#' @inheritParams clusterWithTca
#' @param normMatrix normalized matrix
#' @param days days
#' @param algo algorithm
#' @param takeLog logical; need log?
#'
#' @return TCA object
#'
#' @importFrom EDASeq normCounts
#' @importFrom BoolNet binarizeTimeSeries
#' @importFrom TCseq timeclust
#'
#' @rdname additional-TCA
#'
#' @export
#'
clusterTimes <- function(normMatrix, days, genes, algo="km", clusterK=9, takeLog=FALSE){
  colnames(normMatrix) <- days
  tc <- data.frame(normMatrix, stringsAsFactors = FALSE, check.names = FALSE)
  tc <- sapply(unique(names(tc)),function(col) rowMeans(tc[names(tc) == col]))
  tc <- data.frame(tc, stringsAsFactors = FALSE)
  sheredGenes <- intersect(row.names(tc), genes)
  tc <- tc[sheredGenes, , drop=FALSE]

  tc <- data.matrix(tc)
  if (takeLog)
    tc <- log2(tc+1)
  TCseq::timeclust(tc, algo = algo, k = clusterK, dist = "euclidean", centers = NULL,
            standardize = TRUE)
}

#' Visualize the cluster produced
#'
#' @param tca the results of clusterTimes function
#' @param cols subdivide the column in _col_ cols
#'
#' @return NULL
#' @importFrom TCseq timeclustplot
#'
#' @rdname additional-TCA
#'
#' @export
#'
watchCluster <- function(tca, cols=3) {
  p <- TCseq::timeclustplot(tca, value = "expression", cols = cols)
}
