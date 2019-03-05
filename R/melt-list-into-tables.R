#' Melt list of tables
#'
#' Tables are the output of tables
#'
#' @param list a list of tables
#'
#' @return a data frame with
#'    \item{gene}{the gene name}
#'    \item{cell}{the cell}
#'    \item{cpmWin}{the number of wins it collected}
#'
#' @rdname melt-list-into-tables
#' @export
#'
meltTablesWithListNames <- function(list) {
  out <- NULL
  for (n in names(list)){
    elem <- list[[n]]
    if (!is.null(elem) & !nrow(elem)==0) {
      out <- rbind(out, data.frame(gene=n, cell=names(elem), cmpWin=as.numeric(elem), stringsAsFactors = F))
    }
  }
  out
}

#' Melt list of summarized loosers
#'
#' @inheritParams meltTablesWithListNames
#'
#' @return a data frame with
#'    \item{gene}{the gene name}
#'    \item{cell}{the cell who win}
#'    \item{loosers}{the semi-colon separated list of loosers}
#'
#' @rdname melt-list-into-tables
#' @export
#'
meltTablesWithLoosers <- function(list) {
  out <- NULL
  for (n in names(list)){
    elem <- list[[n]]
    if (!is.null(elem) & !NROW(elem)==0) {
      out <- rbind(out, data.frame(gene=n, cell=names(elem), loosers=elem, stringsAsFactors = F))
    }
  }
  out
}

#' Melt list of vectors
#'
#' Tables are the output of tables
#'
#' @inheritParams meltTablesWithListNames
#'
#' @return a data frame with list content and the name of the element in the column
#'
#' @rdname melt-list-into-tables
#' @export
#'
meltListOfVectors <- function(list) {
  out <- NULL
  for (n in names(list)){
    elem <- list[[n]]
    if (!is.null(elem) & length(elem)!=0 ) {
      out <- rbind(out, data.frame(nm = n, content=elem, stringsAsFactors = F))
    }
  }
  out
}



