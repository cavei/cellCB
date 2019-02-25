#' Convert to Gene
#'
#' @param list a list
#' @return term2gene
#'
#' @rdname gskb-handling
#'
#' @export
#'
convertGskb2Term2Gene <- function(list) {
  lext <- lapply(list, function(x) {
    ext <- cbind(x[1], x[3:length(x)])
    data.frame(ont=ext[,1], gene=ext[,2], stringsAsFactors = F)
  })
  term2gene <- do.call(rbind, lext)
  row.names(term2gene) <- NULL
  term2gene$ont <- factor(term2gene$ont)
  term2gene
}

#' Convert to Name
#'
#' @inheritParams convertGskb2Term2Gene
#' @return a data.frame with term2name
#'
#' @rdname gskb-handling
#'
#' @export
#'
convertGskb2Term2Name <- function(list) {
  lext <- lapply(list, function(x) {
    c(x[1], x[2])
  })
  term2name <- do.call(rbind, lext)
  data.frame(ont=factor(term2name[,1]), name=term2name[,2], stringsAsFactors = F)
}

#' create summarizing barplots from enrichment lists
#'
#' @param enrichmentList a list
#' @param single whether to plot single page plots or cumulative plots
#' @param title suffix names
#' @param ncol the column
#' @param font.size the font size
#' @param showCategory max number of category to show
#'
#' @return a data.frame with term2name
#'
#' @importFrom gridExtra grid.arrange
#' @rdname gskb-handling
#' @export
#'
summaryBarplot <- function(enrichmentList, single=T, title="cluster", ncol=3, font.size = 5, showCategory=100) {
  # library(gridExtra)
  # library(clusterProfiler)

  # if (!is.null(names(enrichmentList))) {
  #   nms <- names(enrichmentList)
  # } else {
  #   nms <- seq_along(enrichmentList)
  # }
  #
  # a <- lapply(nms, function(cluster) {
  #   if (nrow(as.data.frame(enrichmentList[[cluster]]))==0)
  #     return(NULL)
  #   clusterProfiler::barplot(enrichmentList[[cluster]], showCategory=showCategory, font.size = font.size, title = paste(title, cluster))
  # })
  # if (all(sapply(a, is.null)))
  #   return(NULL)
  # a <- a[(!sapply(a, is.null))]
  # if (single)
  #   return(a)
  # gridExtra::grid.arrange(grobs = a, ncol = ncol)
}
