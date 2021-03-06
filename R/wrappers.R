#' @importFrom utils read.table
runWithParameters <- function(dataFile, annotationFile, cell, condition, runWithK=1, runWithFormula="~days + W_1",
                              removeSamples=NULL, varName="days") {

  allCounts <- read.table(dataFile,header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table(annotationFile, header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  if (!is.null(removeSamples)) {
    if (!identical(row.names(allColData), colnames(allCounts)))
      stop("Error counts and colData")
    notFound <- setdiff(removeSamples, row.names(allColData))
    if (length(notFound)>0)
      warning(paste0("The following samples were not found: ", paste(notFound, collapse = ", ")),
              call. = FALSE)
    found <- intersect(removeSamples, row.names(allColData))
    if (length(found)>0) {
      idx <- match(found, row.names(allColData))
      allCounts <- allCounts[,-idx, drop=F]
      allColData <- allColData[-idx, , drop=F]
    }
  }

  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts, annotations = allColData, rcFilter = 100, unvariantNumber=5000)
  dim(cell.data$eset)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)

  DEGS <- extractDegsWithLFC(esetRUV, pADJ.thr = 0.05, logfcthr = 1.5)
  DEGSclusters <- clusterWithTca(esetRUV, DEGS, clusterK = 9)
  list(DEGSclusters=DEGSclusters, DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

#' @importFrom utils read.table
runDignosticWithParameters <- function(cell, condition, runWithK=1, runWithFormula="~days + W_1",
                                       removeSamples=NULL, varName="days") {
  allCounts <- read.table("data/countData.txt",header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table("data/colData.txt", header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  if (!is.null(removeSamples)) {
    if (!identical(row.names(allColData), colnames(allCounts)))
      stop("Error counts and colData")
    notFound <- setdiff(removeSamples, row.names(allColData))
    if (length(notFound)>0)
      warning(paste0("The following samples were not found: ", paste(notFound, collapse = ", ")),
              call. = FALSE)
    found <- intersect(removeSamples, row.names(allColData))
    if (length(found)>0) {
      idx <- match(found, row.names(allColData))
      allCounts <- allCounts[,-idx, drop=F]
      allColData <- allColData[-idx, , drop=F]
    }
  }

  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts, annotations = allColData, rcFilter = 100, unvariantNumber=5000)
  set.seed(1234)
  runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)
}

#' @importFrom utils read.table
runWithParametersTotal <- function(cell, condition, runWithK=1, runWithFormula="~days + W_1",varName="days") {
  allCounts <- read.table("data/total-countData.txt",header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table("data/t-colData.txt", header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts, annotations = allColData, rcFilter = 100, unvariantNumber=5000)
  dim(cell.data$eset)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all",varName=varName)

  DEGS <- extractDegsWithLFC(esetRUV, pADJ.thr = 0.05, logfcthr = 1.5)
  DEGSclusters <- clusterWithTca(esetRUV, DEGS, clusterK = 9)
  list(DEGSclusters=DEGSclusters, DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

# save(fap, file=paste0("Rdatas/",cell,"-",condition,".RData"))

#' @importFrom utils read.table
runDignosticWithParametersTotal <- function(cell, condition, runWithK=1, runWithFormula="~days + W_1",varName="days") {
  allCounts <- read.table("data/total-countData.txt",header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table("data/t-colData.txt", header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts, annotations = allColData, rcFilter = 100, unvariantNumber=5000)
  set.seed(1234)
  runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)
}

runWithParametersCell <- function(cell, condition, runWithK=1, runWithFormula="~days + W_1", removeSamples=NULL, varName="cell") {
  allCounts <- read.table("data/countData.txt",header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table("data/colData.txt", header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  if (!is.null(removeSamples)) {
    if (!identical(row.names(allColData), colnames(allCounts)))
      stop("Error counts and colData")
    notFound <- setdiff(removeSamples, row.names(allColData))
    if (length(notFound)>0)
      warning(paste0("The following samples were not found: ", paste(notFound, collapse = ", ")),
              call. = FALSE)
    found <- intersect(removeSamples, row.names(allColData))
    if (length(found)>0) {
      idx <- match(found, row.names(allColData))
      allCounts <- allCounts[,-idx, drop=F]
      allColData <- allColData[-idx, , drop=F]
    }
  }

  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts,
                                          annotations = allColData, rcFilter = 100,
                                          unvariantNumber=5000, varName=varName)
  dim(cell.data$eset)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)

  DEGS <- getDEGsAcrossPopulation(esetRUV, pADJ.thr = 0.05, logfcthr = 1.5)
  list(DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

#' @importFrom utils read.table
#' @importFrom EDASeq counts
.computeRPKM <- function(cell, condition, gene_length="../m38.p5-ensembl-91-renge-length/gene_length.txt") {
  warning("Deprecated version")
  allCounts <- read.table("data/countData.txt",header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table("data/colData.txt", header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  cell.data <- filterExpression(cell = cell, condition = condition, expr = allCounts, annotations = allColData,
                                          rcFilter = NULL, unvariantNumber=NULL)
  lengths <- read.table(gene_length, header=T, row.names=1,
                        sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  usableGenes <- intersect(row.names(cell.data$eset), lengths$gene)
  lengths <- lengths[lengths$gene %in% usableGenes, ]
  maxMeanGeneLength <- tapply(seq_len(NROW(lengths)), lengths$gene, function(idx) max(lengths$mean[idx]))

  expr <- EDASeq::counts(cell.data$eset)[usableGenes, ]
  maxMeanGeneLength <- maxMeanGeneLength[usableGenes]
  y <- DGEList(counts=expr)
  y$genes <- data.frame(Length=as.numeric(maxMeanGeneLength))
  edgeR::rpkm(y, gene_length="Length")
}


#' Run a RUV with parameters
#'
#' @param dataFile the name of the file that contains the expression dataset
#' @param annotationFile the name of the file that contains annotations
#' @param cell the cell sample name to pick
#' @param condition the condition
#' @param runWithK run RUV using k covariates
#' @param runWithFormula a string coercible with as.formula
#' @param removeSamples remove the listed samples
#' @param varName the variable to used to evaluate unwanted variance
#' @param pADJ.thr the thresholf for pvalue adjusted
#' @param logfcthr absolute log fold change
#' @param rcFilter filter expression for genes minimal rowCounts in atLeastIn samples
#' @param atLeastIn filter expression for genes minimal rowCounts in atLeastIn samples
#' @param rpkmMean filter expression for genes with an average rpkm of 3
#' @param clusterK try clusterize according to k clusters
#' @param gene_length filename of file containg gene lengths
#'
#' @return DEGSclusters, DEGS, esetRUV and cell.data
#'
#' @importFrom utils read.table
#' @export
#'
runWithParametersRpkm <- function(dataFile, annotationFile, cell, condition, runWithK=1, runWithFormula="~days + W_1",
                                  removeSamples=NULL, varName="days", pADJ.thr = 0.05, logfcthr = 1.5,
                                  rcFilter = 100, atLeastIn = 2,rpkmMean = 3, clusterK = 9, gene_length="gene_length.txt") {

  allCounts <- read.table(dataFile, header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  allColData <- read.table(annotationFile, header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  if (!is.null(removeSamples)) {
    if (!identical(row.names(allColData), colnames(allCounts)))
      stop("Error counts and colData")
    notFound <- setdiff(removeSamples, row.names(allColData))
    if (length(notFound)>0)
      warning(paste0("The following samples were not found: ", paste(notFound, collapse = ", ")),
              call. = FALSE)
    found <- intersect(removeSamples, row.names(allColData))
    if (length(found)>0) {
      idx <- match(found, row.names(allColData))
      allCounts <- allCounts[,-idx, drop=F]
      allColData <- allColData[-idx, , drop=F]
    }
  }

  cell.data <- filterExpressionByRpkmMean(cell = cell, condition = condition, expr = allCounts, annotations = allColData,
                                          rcFilter = rcFilter, atLeastIn = atLeastIn, rpkmMean = rpkmMean,
                                          unvariantNumber=3000,
                                          varName=varName,
                                          gene_length=gene_length)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)

  DEGS <- extractDegsWithLFC(esetRUV, pADJ.thr = pADJ.thr, logfcthr = logfcthr)
  DEGSclusters <- clusterWithTca(esetRUV, DEGS, clusterK = clusterK)
  list(DEGSclusters=DEGSclusters, DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

#' Run in total muscle [DEPRECATED]
#'
#' @inheritParams runWithParametersRpkm
#'
#' @importFrom utils read.table
runWithParametersTotalRpkm <- function(dataFile, annotationFile, cell, condition, runWithK=1, runWithFormula="~days + W_1",
                                       varName="days", pADJ.thr = 0.05, logfcthr = 1.5,
                                       rcFilter = 100, atLeastIn = 2,rpkmMean = 3,
                                       clusterK = 9, gene_length="gene_length.txt") {
  warning("Deprecated")
  allCounts <- read.table(dataFile ,header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table(annotationFile, header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  cell.data <- filterExpressionByRpkmMean(cell = cell, condition = condition, expr = allCounts, annotations = allColData,
                                          rcFilter = rcFilter, atLeastIn = atLeastIn, rpkmMean = rpkmMean,
                                          unvariantNumber=3000,
                                          varName=varName, gene_length=gene_length)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all",varName=varName)

  DEGS <- extractDegsWithLFC(esetRUV, pADJ.thr = pADJ.thr, logfcthr = logfcthr)
  DEGSclusters <- clusterWithTca(esetRUV, DEGS, clusterK = clusterK)
  list(DEGSclusters=DEGSclusters, DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

#' @importFrom utils read.table
runWithParametersCellRpkm <- function(dataFile, annotationFile, cell, condition, runWithK=1, runWithFormula="~days + W_1", removeSamples=NULL,
                                      varName="cell", pADJ.thr = 0.05, logfcthr = 1.5, rcFilter = 100, atLeastIn = 2,rpkmMean = 3,
                                      clusterK = 9, gene_length="gene_length.txt") {

  warning("Deprecated")
  allCounts <- read.table(dataFile,header=T, row.names=1,
                          sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  allColData <- read.table(annotationFile, header=T, row.names=1,
                           sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

  if (!is.null(removeSamples)) {
    if (!identical(row.names(allColData), colnames(allCounts)))
      stop("Error counts and colData")
    notFound <- setdiff(removeSamples, row.names(allColData))
    if (length(notFound)>0)
      warning(paste0("The following samples were not found: ", paste(notFound, collapse = ", ")),
              call. = FALSE)
    found <- intersect(removeSamples, row.names(allColData))
    if (length(found)>0) {
      idx <- match(found, row.names(allColData))
      allCounts <- allCounts[,-idx, drop=F]
      allColData <- allColData[-idx, , drop=F]
    }
  }

  cell.data <- filterExpressionByRpkmMean(cell = cell, condition = condition, expr = allCounts,
                                          annotations = allColData, rcFilter = rcFilter,
                                          atLeastIn = atLeastIn, rpkmMean = rpkmMean,
                                          unvariantNumber=3000, varName=varName, gene_length=gene_length)

  set.seed(1234)
  esetRUV <- runWithDiagnosticPlots(cell.data, k = runWithK, formula=runWithFormula, type="all", varName=varName)

  DEGS <- getDEGsAcrossPopulation(esetRUV, pADJ.thr = pADJ.thr, logfcthr = logfcthr)
  list(DEGS=DEGS, esetRUV=esetRUV, cell.data=cell.data)
}

# runMarkdown <- function(file, params=NULL) {
#   # params=list(region = "west", start = as.Date("2015-02-01"))
#   rmarkdown::render(file, params = params)
# }
