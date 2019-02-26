#' Merge columns
#'
#' @param li a list objects dataframe like
#' @param col the name of the column
#' @param cutOff set to NA if lower than cutOff
#'
#' @return a matrix like
#'
#' @export
mergeColumns <- function(li, col="p.adjust", cutOff=0.05) {
  li <- resolveAndOrder(li)
  mat <- do.call(cbind, lapply(li, function(o) o[,col]))
  mat <- apply(mat, 2, as.numeric)
  row.names(mat) <- row.names(li[[1]])
  sel <- apply(na2false(mat<=cutOff), 1, any)
  mat[sel, , drop=F]
}

resolveAndOrder <- function(li) {
  allnames <- unname(unique(unlist(lapply(li, row.names))))
  lapply(li, function(o) {
    absents <- setdiff(allnames, row.names(o))
    ab <- as.data.frame(matrix(NA, nrow=length(absents), ncol=ncol(o)))
    row.names(ab) <- absents
    colnames(ab) <- colnames(o)
    out <- rbind(o, ab)
    out[allnames, , drop=F]
  })
}

na2false <- function(x) {
  x[is.na(x)] <- FALSE
  x
}

na2true <- function(x) {
  x[is.na(x)] <- TRUE
  x
}

cliPvalues <- function(x, value=0.05) {
  x[na2true(x > value)] <- 1
  x
}
