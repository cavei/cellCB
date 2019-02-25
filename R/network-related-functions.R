#' Given a list of ligand receptor it extracts the network
#'
#' @param lrday a list with ligand and receptors
#' @param graphNEL a graphNEL object representing the network
#'
#' @return a list
#'    \item{g}{the graph}
#'    \item{activeEdges}{the active edges}
#'    \item{autocrine}{connected nodes}
#'
#' @importFrom reshape2 melt
#' @importFrom graph subGraph degree
#' @rdname network-handling
#' @export
#'
extractDayLigandReceptorNetwork <- function(lrday, graphNEL) {
  nodes <- unique(do.call(c, lrday))
  ig <- graph::subGraph(nodes, graphNEL)
  nodesDegree <- graph::degree(ig)
  if (identical(names(nodesDegree$inDegree), names(nodesDegree$outDegree)))
    nd <- nodesDegree$inDegree + nodesDegree$outDegree
  connectedNodes <- names(which(nd != 0))
  g = graph::subGraph(connectedNodes, ig)
  activeEdges <- reshape2::melt(graph::edges(g))
  activeEdges <- cbind(src=as.character(activeEdges[,2]), dest=as.character(activeEdges[,1]))
  list(g = g, activeEdges=activeEdges, autocrine=connectedNodes)
}

#' Given a list of ligand receptor it extracts the network
#'
#' @inheritParams extractDayLigandReceptorNetwork
#'
#' @return a network
#' @importFrom graph subGraph degree
#' @importFrom network as.network set.network.attribute set.vertex.attribute
#' @importFrom methods as
#'
#' @rdname network-handling
#' @export
#'
extractDayNetwork <- function(lrday, graphNEL) {
  nodes <- unique(do.call(c, lrday))
  ig <- graph::subGraph(nodes, graphNEL)
  nodesDegree <- graph::degree(ig)
  if (identical(names(nodesDegree$inDegree), names(nodesDegree$outDegree)))
    nd <- nodesDegree$inDegree + nodesDegree$outDegree
  connectedNodes <- names(which(nd != 0))
  g = graph::subGraph(connectedNodes, ig)
  adjMat <- as(g, "matrix")
  net <- as.network(adjMat)
  network::set.network.attribute(net, "vertex.pid", row.names(adjMat))
  network::set.vertex.attribute(net, "vertex.pid", row.names(adjMat))
  net
}

getDailyNodeSize <- function(mat) {
  ma <- apply(mat, 2, computeNodeSize)
  row.names(ma) <- row.names(mat)
  ma
}

computeNodeSize <- function(x) {
  # size <- cut(x, breaks = quantile(unique(x), probs = seq(0, 1, 0.1)),
  #   include.lowest = TRUE, labels = 1:10)
  size <- cut(x, breaks = c(0,5,55,Inf),include.lowest = TRUE, labels = c(1,3,10))
  as.numeric(as.character(size))
}

