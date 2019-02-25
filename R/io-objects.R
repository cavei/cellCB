#' Collect DEGs Pvalues with LogFold Change from RData.
#'
#' @param esetRUV an expression set
#' @param pADJ.thr the threshold for adjusted pvalue, default 0.05
#' @param logfcthr the absolute threshold value for lfc, default 1
#'
#' @return a three column matix with gene, pvalue and lfc
#'
#' @rdname io-object
#' @export
#'
extractDegsPvaluesWithLFC <- function(esetRUV, pADJ.thr=0.05, logfcthr=1) {
  logfcsIdx <- grep("logFC", colnames(esetRUV$DEGs$table))

  logfcFilter <- apply(esetRUV$DEGs$table[, logfcsIdx,drop=F], 1, anyLogFC, lfcThr=logfcthr)
  padjFilter <- esetRUV$DEGs$table$FDR <= pADJ.thr
  selectGenes <- padjFilter & logfcFilter

  esetRUV$DEGs$table[selectGenes, c("PValue","FDR"), drop=F]
}

#' Collect DEGs LogFold Change from RData.
#' @inheritParams extractDegsPvaluesWithLFC
#'
#' @return a three column matix with gene, pvalue and lfc
#'
#' @rdname io-object
#' @export
#'
extractDegsWithLFC <- function(esetRUV, pADJ.thr=0.05, logfcthr=1) {
  logfcsIdx <- grep("logFC", colnames(esetRUV$DEGs$table))

  logfcFilter <- apply(esetRUV$DEGs$table[, logfcsIdx,drop=F], 1, anyLogFC, lfcThr=logfcthr)
  padjFilter <- esetRUV$DEGs$table$FDR <= pADJ.thr
  selectGenes <- padjFilter & logfcFilter

  genes <- row.names(esetRUV$DEGs$table[selectGenes, , drop=F])
  genes
}

#' Load and collect Active at Days from Obj.
#'
#' @param obj expression set
#' @param adjusted logical, if true retrieve adjusted pvalues
#'
#' @return pvalue vector
#'
#' @rdname io-object
#' @export
#'
getDEGsPvaluesFromObj <- function(obj, adjusted=T) {
  # (extractDegsPvaluesWithLFC(fap_ko$esetRUV, pADJ.thr = 0.01, logfcthr = 2))
  table <- obj$esetRUV$DEGsTable$table
  if (adjusted) {
    pvalue <- as.numeric(table$FDR)
  } else {
    pvalue <- as.numeric(table$PValue)
  }
  names(pvalue) <- row.names(table)
  pvalue
}

#' Load and collect Active at Days from Obj.
#'
#' @inheritParams getDEGsPvaluesFromObj
#'
#' @return a list of genes active at Days.
#'
#' @rdname io-object
#' @export
#'
getActiveAtDaysFromObj <- function(obj) {
  # library(findAWay)
  clusterCell <- obj$DEGSclusters$tca@cluster
  binCell <- obj$DEGSclusters$binaryCenter$binarizedMeasurements
  cellExpressionAtDays <- lapply(colnames(binCell), defineCellExpression, binCell = binCell,
                                 clusterCell = clusterCell)
  names(cellExpressionAtDays) <- colnames(binCell)
  lapply(cellExpressionAtDays, function(x) x$activeAtDay)
}


#' Define Cell Expression
#'
#' Formely in findAWay R package
#'
#' @param day the day
#' @param binCell binary cluster representation
#' @param clusterCell the gene organized in clusters
#' @param measurableGenes the total measurable genes
#'
#' @return list with activeAtDay, non_avtiveAtDay and constitutive
#'
#' @export
defineCellExpression <- function(day, binCell, clusterCell, measurableGenes=NULL) {
  constitutive = NULL
  dayi <- which(colnames(binCell)==day) # binCell
  if (length(dayi)==0)
    return(list(activeAtDay=NULL, non_avtiveAtDay=NULL, constitutive=NULL))
  singleDaySel <- binCell[,dayi] > 0
  names(singleDaySel) <- row.names(binCell)
  not_singleDaySel <- !singleDaySel

  singleDaySel <- as.numeric(names(which(singleDaySel)))
  not_singleDaySel <- as.numeric(names(which(not_singleDaySel)))

  genes_of_the_day <- names(clusterCell[clusterCell %in% singleDaySel])
  genes_nonActive <- names(clusterCell[clusterCell %in% not_singleDaySel])

  if (!is.null(measurableGenes)) {
    genes_of_the_day <- intersect(genes_of_the_day, measurableGenes)
    genes_nonActive <- intersect(genes_nonActive, measurableGenes)
    constitutive <- setdiff(measurableGenes, c(genes_of_the_day, genes_nonActive))
  }

  list(activeAtDay=genes_of_the_day, non_avtiveAtDay=genes_nonActive, constitutive=constitutive)
}

