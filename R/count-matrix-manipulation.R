#' Average Replicates Gene Wise
#'
#' @param matrix with columns names from which infer the factors
#'
#' @return averaged matrix
#'
#' @export
#'
averageReplicate <- function(matrix) {
  replicates <- detectReplicates(colnames(matrix))
  l <- tapply(seq_along(replicates), replicates, function(idx) {
    rowMeans(matrix[, idx, drop=F])
  })
  am <- do.call(cbind, l)
  colnames(am) <- paste0("d",colnames(am))
  am
}

detectReplicates <- function(names) {
  sn <- strsplit(names, "_")
  sn <- sapply(sn, function(x) x[2])
  factor(sn, levels=as.character(sort(as.numeric(unique(sn)))))
}

