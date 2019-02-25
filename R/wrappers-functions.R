getDEGsAcrossPopulation <- function(esetRUV, pADJ.thr=0.05, logfcthr=1, varName="days") {
  logfcsIdx <- grep("logFC", colnames(esetRUV$DEGsTable$table))
  logfcFilter <- apply(esetRUV$DEGsTable$table[, logfcsIdx,drop=F], 1, anyLogFC, lfcThr=logfcthr)
  padjFilter <- esetRUV$DEGsTable$table$FDR <= pADJ.thr

  selectGenes <- padjFilter & logfcFilter
  genes <- esetRUV$DEGsTable$table[selectGenes, "FDR", drop=F]
  genes
}

#' @importFrom graphics hist par
#' @importFrom stats as.formula
#' @importFrom utils write.table
#'
fullAnalysis <- function(cell, condition, allCounts, allColData, clusterK=9, rcFilter=100, unvariantNumber=5000, k=1,
                         formula = "~days + W_1", pADJ.thr=0.05, logfcthr=1, writeRes=F, varName="days") {

  cellCond <- filterExpression(cell=cell, condition=condition,
                               expr = allCounts, annotations = allColData,
                               unvariantNumber = unvariantNumber, rcFilter = rcFilter, varName=varName)

  esetRUV <- edgeR.anovaOnTimeWithRUV(cellCond$eset, cellCond$unvariant, k = k, formula = formula, varName=varName)

  normalizedCountsMatrix <- normCounts(esetRUV$setEmpirical)
  logfcsIdx <- grep("logFC", colnames(esetRUV$DEGs$table))

  logfcFilter <- apply(esetRUV$DEGs$table[, logfcsIdx,drop=F], 1, anyLogFC, lfcThr=logfcthr)
  padjFilter <- esetRUV$DEGs$table$FDR <= pADJ.thr

  selectGenes <- padjFilter & logfcFilter
  genes <- row.names(esetRUV$DEGs$table[selectGenes, , drop=F])

  tca <- clusterTimes(normMatrix = normalizedCountsMatrix, days = as.character(esetRUV$setEmpirical$days),
                      genes = genes, clusterK = 9, takeLog = T)

  binaryCenter <- binarizeTimeSeries(tca@centers)
  clusters <- data.frame(gene=names(tca@cluster), cluster=tca@cluster)

  ## Write results
  if(writeRes) {
    write.table(normalizedCountsMatrix, file=paste0(cell,"-",condition,"-filter",rcFilter,"-normalized-RUV.txt"), sep="\t", quote=F)
    write.table(clusters, file=paste0(cell,"-",condition,"-clusters.txt"), col.names=F, row.names=F, sep="\t", quote=F)
    write.table(tca@centers, file=paste0(cell,"-",condition,"-centers.txt"), sep="\t", quote=F)
    write.table(tca@data, file=paste0(cell,"-",condition,"-data.txt"), sep="\t", quote=F)
    write.table(binaryCenter$binarizedMeasurements, file=paste0(cell,"-",condition,"-centers-binary.txt"), sep="\t", quote=F)
  }

  return(list(tca=tca, binaryCenter=binaryCenter, esetRUV=esetRUV, normalizedCounts = normalizedCountsMatrix))
}

#' @importFrom graphics hist par
#' @importFrom stats as.formula
#' @importFrom utils write.table
#' @importFrom EDASeq normCounts
#'
getCrossPopulationDEGS <- function(cell, condition, allCounts, allColData, rcFilter=100, unvariantNumber=5000, k=1,
                                   formula = "~days + W_1", pADJ.thr=0.05, logfcthr=1, writeRes=F, varName="days") {

  cellCond <- filterExpression(cell=cell, condition=condition,
                               expr = allCounts, annotations = allColData,
                               unvariantNumber = unvariantNumber, rcFilter = rcFilter, varName=varName)

  esetRUV <- edgeR.anovaOnTimeWithRUV(cellCond$eset, cellCond$unvariant, k = k, formula = formula, varName=varName)

  logfcsIdx <- grep("logFC", colnames(esetRUV$DEGs$table))

  logfcFilter <- apply(esetRUV$DEGs$table[, logfcsIdx,drop=F], 1, anyLogFC, lfcThr=logfcthr)
  padjFilter <- esetRUV$DEGs$table$FDR <= pADJ.thr

  selectGenes <- padjFilter & logfcFilter
  genesFdr <- esetRUV$DEGs$table[selectGenes, "FDR", drop=F]
  normalizedCountsMatrix <- EDASeq::normCounts(esetRUV$setEmpirical)
  ## Write results
  if(writeRes) {
    write.table(normalizedCountsMatrix, file=paste0(cell,"-",condition,"-filter",rcFilter,"-normalized-RUV-notca.txt"), sep="\t", quote=F)
    write.table(genesFdr, file=paste0(cell,"-",condition,"-filter",rcFilter,"-normalized-RUV-degs.txt"), sep="\t", quote=F)
  }

  return(list(esetRUV=esetRUV, degs=genesFdr, normalizedCounts = normalizedCountsMatrix))
}

#' @importFrom graphics hist par
#' @importFrom stats as.formula
#' @importFrom EDASeq plotRLE plotPCA
#'
plotDiagnostic <- function(cellCond, k=1, formula="~days + W_1",
                           type=c("all", "plotRLE", "plotPCA", "histPvalue"), varName="days") {

  esetRUV <- edgeR.anovaOnTimeWithRUV(cellCond$eset, cellCond$unvariant, k = k,
                                      formula = formula, varName=varName)

  cls <- as.factor(pData(esetRUV$setEmpirical)[, varName])
  colors <- brewer.pal(length(unique(cls)), "Spectral")
  if (type[1]=="all") {
    par(mfrow=c(1,1))
    EDASeq::plotRLE(esetRUV$setEmpirical, outline=FALSE, col=colors[cls], las=3, cex.axis=0.6)
    EDASeq::plotPCA(esetRUV$setEmpirical, col=colors[cls])
    hist(esetRUV$DEGs$table$PValue)
    # meanVarPlot(esetRUV$setEmpirical, log=TRUE)
  } else if (type[1]=="plotRLE") {
    par(mfrow=c(1,1))
    EDASeq::plotRLE(esetRUV$setEmpirical, outline=FALSE, col=colors[cls], las=3, cex.axis=0.6)
  } else if (type[1]=="plotPCA") {
    par(mfrow=c(1,1))
    EDASeq::plotPCA(esetRUV$setEmpirical, col=colors[cls])
  } else if (type[1]=="histPvalue") {
    par(mfrow=c(1,1))
    hist(esetRUV$DEGs$table$PValue)
  } else {
    stop("Unkonow plot.")
  }
}

#' @importFrom graphics hist par
#' @importFrom EDASeq plotRLE plotPCA
#' @importFrom RColorBrewer brewer.pal
runWithDiagnosticPlots <- function(cellCond, k=1, formula="~days + W_1",
                                   type=c("all", "plotRLE", "plotPCA", "histPvalue"), varName="days") {

  if (k==0) {
    esetRUV <- edgeR.anovaOnTime(cellCond$eset, formula = formula, varName=varName)
  } else {
    esetRUV <- edgeR.anovaOnTimeWithRUV(cellCond$eset, cellCond$unvariant, k = k,
                                        formula = formula, varName=varName)
  }
  cls <- as.factor(pData(esetRUV$setEmpirical)[, varName])
  colors <- RColorBrewer::brewer.pal(length(unique(cls)), "Spectral")
  if (type[1]=="all") {
    par(mfrow=c(1,1))
    EDASeq::plotRLE(esetRUV$setEmpirical, outline=FALSE, col=colors[cls], las=3, cex.axis=0.6)
    EDASeq::plotPCA(esetRUV$setEmpirical, col=colors[cls])
    hist(esetRUV$DEGs$table$PValue)
    # meanVarPlot(esetRUV$setEmpirical, log=TRUE)
  } else if (type[1]=="plotRLE") {
    par(mfrow=c(1,1))
    EDASeq::plotRLE(esetRUV$setEmpirical, outline=FALSE, col=colors[cls], las=3, cex.axis=0.6)
  } else if (type[1]=="plotPCA") {
    par(mfrow=c(1,1))
    EDASeq::plotPCA(esetRUV$setEmpirical, col=colors[cls])
  } else if (type[1]=="histPvalue") {
    par(mfrow=c(1,1))
    hist(esetRUV$DEGs$table$PValue)
  } else {
    stop("Unkonow plot.")
  }
  invisible(esetRUV)
}
