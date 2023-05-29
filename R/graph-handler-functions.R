#' From edge matrix to edge List
#'
#' @param edgeMat edge matrix
#'
#' @export
edgeList <- function(edgeMat){
  if (!(is.matrix(edgeMat))){
    return(edgeMat)
  }
  if (NCOL(edgeMat) < 2)
    stop("edge matrix do not have 2 columns.")

  if(NCOL(edgeMat)>2) {
    warning("More than 2 columns found. Only col 1 and 2 will be used.")
    edgeMat <- edgeMat[,1:2]
  }

  edgeVector <- NULL
  for (i in 1:NROW(edgeMat)){
    edgeVector<-c(edgeVector,edgeMat[i,])
  }
  return(edgeVector)
}
