combineActiveAndConstitutive <- function(x) {
  if (is.character(x))
    return(x)
  c(x$activeAtDay, x$constitutive, x$constitutive2)
}

#' @importFrom network as.edgelist
createSifOFLigandReceptor <- function(ligandsReceptorCumulative, receptorGraph,
                                      cell, condition) {
  # library(network)
  removeTimes <- NULL
  for (x in names(ligandsReceptorCumulative)) {
    l <- unlist(lapply(ligandsReceptorCumulative[[x]], length))
    if (any(l==0))
      removeTimes <- c(removeTimes, x)
  }
  if (!is.null(removeTimes))
    ligandsReceptorCumulative[[removeTimes]] <- NULL

  networklist <- lapply(ligandsReceptorCumulative, extractDayNetwork, graphNEL=receptorGraph)
  networkEdgeList <- lapply(networklist, as.edgelist)
  lapply(networkEdgeList, function(edgesMat) {
    nodes <- attr(edgesMat, "vnames")
    data.frame(src=nodes[edgesMat[,1]], dest=nodes[edgesMat[,2]], stringsAsFactors = F)
  })
}

loadCellRData <- function(RData, objname) {
  load(RData)
  get(objname)
}
