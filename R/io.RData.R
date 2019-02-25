#' Load and extract normalized expression from RData.
#'
#' @param RData filnename with RData
#' @param objname the name of the loaded object
#'
#' @return normalized gene expression matrix
#'
#' @rdname io-RData
#' @importFrom EDASeq normCounts
#' @export
#'
getNormExpressionFromRData <- function(RData, objname) {
  # library(findAWay)
  load(RData)
  loadObj <- get(objname)
  EDASeq::normCounts(loadObj$esetRUV$setEmpirical)
}

#' Load and collect RPKM from RData.
#'
#' @inheritParams getNormExpressionFromRData
#'
#' @return rpkm expression matrix
#'
#' @rdname io-RData
#' @export
#'
getRpkmExpressionFromRData <- function(RData, objname) {
  # library(findAWay)
  load(RData)
  loadObj <- get(objname)
  loadObj$cell.data$rpkm
}

#' Load and collect DEGs Pvalues with LogFold Change from RData.
#'
#' @inheritParams getNormExpressionFromRData
#' @param pADJ.thr the threshold for adjusted pvalue, default 0.05
#' @param logfcthr the absolute threshold value for lfc, default 1.5
#'
#' @return a three column matix with gene, pvalue and lfc
#'
#' @rdname io-RData
#' @export
#'
extractDegsPvaluesWithLFCFromRData <- function(RData, objname, pADJ.thr = 0.05, logfcthr = 1.5) {
  load(RData)
  loadObj <- get(objname)
  extractDegsPvaluesWithLFC(loadObj$esetRUV, pADJ.thr, logfcthr)
}

## extractDegsPvaluesWithLFC ??????????

#' Load and collect Active at Days from RData.
#'
#' @inheritParams getNormExpressionFromRData
#'
#' @return a list of genes Expressed at a given Day.
#'
#' @rdname io-RData
#' @export
#'
getDEGsPvaluesFromRData <- function(RData, objname) {
  # library(findAWay)
  load(RData)
  loadObj <- get(objname)
  getDEGsPvaluesFromObj(loadObj)
}

#' Load and collect Active Genes at Days from RData.
#'
#' @inheritParams getNormExpressionFromRData
#'
#' @return a list of genes active at a given Day.
#'
#' @rdname io-RData
#' @export
#'
getActiveAtDaysFromRData <- function(RData, objname) {
  # library(findAWay)
  load(RData)
  loadObj <- get(objname)
  getActiveAtDaysFromObj(loadObj)
}
