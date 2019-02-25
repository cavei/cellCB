#' Extract the contrasts with a given log fold change (lfc) from anova-like (edgeR) lfc table
#'
#' @param x a vector of logFoldChanges
#' @param cmp the vector of contrasts, typically colnames of lfc table
#' @param lfcThr the lfc threshold (in absolute value)
#'
#' @return a data.frame with
#'    \item{gene}{the gene}
#'    \item{num}{the numerator of the contrast}
#'    \item{den}{the denominator of the contrast}
#'    \item{winner}{the winner of the contrast}
#'    \item{looser}{the loosers of the contrast}
#'
#' @rdname log-fold-change-handling
#'
#' @importFrom stats dist
#'
#' @export
anyCmpSigByLogFC <- function(x, cmp, lfcThr=2) {
  v <- as.numeric(x)
  idx <- which(abs(v) >= lfcThr)
  if (length(idx)!=0) {
    first <- data.frame(lfc= v[idx],
                        num=sub(pattern = "logFC.cell", "", cmp[idx]),
                        den="ec",
                        stringsAsFactors = F)
    winners <- first$num; winners[first$lfc < 0] <- "ec"
    loosers <- first$num; loosers[first$lfc > 0] <- "ec"

    first$winner <- winners
    first$looser <- loosers
  } else {
    first<-NULL
  }

  d <- as.matrix(dist(cbind(v,1)))
  d[lower.tri(d)]<-NA
  idxCalc <- which(d >= lfcThr, arr.ind = T)
  if (nrow(idxCalc)!=0) {
    add <- unname(apply(idxCalc, 1, function(i) {
      num <- i[1]; den=i[2]
      data.frame(lfc=v[num]-v[den],
                 num=sub(pattern = "logFC.cell", "", cmp[num]),
                 den=sub(pattern = "logFC.cell", "", cmp[den]),
                 winner = ifelse(v[num]-v[den] <0,
                                 sub(pattern = "logFC.cell", "", cmp[den]),
                                 sub(pattern = "logFC.cell", "", cmp[num])),
                 looser = ifelse(v[num]-v[den] <0,
                                 sub(pattern = "logFC.cell", "", cmp[num]),
                                 sub(pattern = "logFC.cell", "", cmp[den])),
                 stringsAsFactors = F)
    }))
    add <- do.call(rbind, add)
  } else {
    add <- NULL
  }
  rbind(first,add)
}


#' Extract All Genes with given log fold change (lfc) from anova-like lfc table
#'
#' @inheritParams anyCmpSigByLogFC
#'
#' @return TRUE / FALSE
#' @rdname log-fold-change-handling
#'
#' @importFrom stats dist
#'
#' @export
#'
anyLogFC <- function(x, lfcThr=2) {
  v <- as.numeric(x)
  any(abs(v) >= lfcThr) | any (dist(cbind(v,1)) >= lfcThr)
}
