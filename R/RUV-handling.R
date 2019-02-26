#' Create an expression set using EDASeq function
#'
#' @param obj a list with $e expression and $ann annotation
#'
#' @return expression set (eset)
#'
#' @importFrom EDASeq newSeqExpressionSet
#' @rdname RUV-handling
#'
createSet <- function(obj) {
  obj$ann$days <- as.factor(paste0("d",obj$ann$time))
  pData <- obj$ann
  EDASeq::newSeqExpressionSet(obj$e, phenoData=pData)
}

#' Create an expression set using EDASeq function
#'
#' @param expression expression data
#' @param annotation annotation data
#'
#' @return expression set (eset)
#'
#' @importFrom EDASeq newSeqExpressionSet
#'
#' @rdname RUV-handling
#' @export
#'
createSeqExpressionSet <- function(expression, annotation) {
  if (is.null(annotation$time))
    stop("Annotation must have a 'time' column.")
  annotation$days <- as.factor(paste0("d",annotation$time))
  newSeqExpressionSet(expression, phenoData=annotation)
}

#' Compute empirically invariant genes
#'
#' @param set an expression set
#' @param number the number of desired invariant genes
#' @inheritParams filterExpression
#'
#' @return a vector of invariant gene
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp glmFit glmLRT
#' @importFrom stats model.matrix
#' @importFrom EDASeq counts
#' @importFrom utils tail
#'
#' @rdname RUV-handling
#'
#' @export
#'
computeEmpiricalUnvariants <- function(set, number=5000, varName="days") {
  design <- model.matrix(as.formula(paste0("~", varName)), data=pData(set))
  y <- edgeR::DGEList(counts=counts(set), group=pData(set)[,varName])
  y <- edgeR::calcNormFactors(y, method="upperquartile")
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  fit <- edgeR::glmFit(y, design)
  lrt <- edgeR::glmLRT(fit, coef=2:NCOL(design))
  top <- edgeR::topTags(lrt, n=nrow(set))$table
  number=min(number, NROW(top))
  row.names(tail(top, number))
}

#' Compute residuals
#'
#' @inheritParams computeEmpiricalUnvariants
#'
#' @return a vector residuals
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp glmFit
#' @importFrom stats model.matrix residuals
#' @importFrom EDASeq counts
#' @rdname RUV-handling
#'
#' @export
#'
computeResiduals <- function(set) {
  design <- model.matrix(~days, data=pData(set))
  y <- DGEList(counts=EDASeq::counts(set), group=pData(set)$days)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  residuals(fit, type="deviance")
}

#' Filter expression
#'
#' @inheritParams filterExpression
#'
#' @return a list with filtered expression, filtered annotation
#'
#' @importFrom EDASeq betweenLaneNormalization
#'
#' @rdname RUV-handling
#'
#' @export
#'
filterExps <- function(cell, condition, expr, annotations, rcFilter=NULL) {
  warning("Function deprecated")
  sel <- annotations$cell == cell & annotations$condition == condition
  a <- annotations[sel, , drop=F]
  e <- as.numeric(expr[, row.names(a), drop=F])

  if (!is.null(rcFilter)) {
    keep <- apply(e>=rcFilter, 1, any)
    e <- e[keep, , drop=F]
  }
  eNorm <- EDASeq::betweenLaneNormalization(e, which="upper")
  return(list(
    e = e,
    ann = a,
    eNorm = eNorm,
    methods = paste0("normalization=upper\n",
                     "cellType=",cell,
                     "\ncondition=",condition,
                     "\ngenes with less then ", rcFilter, " in one sample were filtered out.")
  ))
}

#' Filter expression
#'
#' @param cell string for cell name (can be ALL or a vector)
#' @param condition string for cell condition (can be a vector with more than one)
#' @param expr data matrix
#' @param annotations annotation data referred to the matrix
#' @param rcFilter minimal read count accepted in atLeastIn samples
#' @param atLeastIn default 1 sample
#' @param unvariantNumber number of invariant genes; default 5000
#' @param varName the variable that set up the contrast to compute invariants; default "days"
#'
#' @return a filtered eset and other details
#'
#' @importFrom EDASeq betweenLaneNormalization
#' @rdname RUV-handling
#'
#' @export
#'
filterExpression <- function(cell, condition, expr, annotations, rcFilter=NULL,
                             atLeastIn=1, unvariantNumber=5000, varName="days") {

  if (cell=="ALL") {
    sel <- annotations$condition %in% condition
  } else {
    sel <- annotations$cell %in% cell & annotations$condition %in% condition
  }

  a <- annotations[sel, , drop=F]
  e <- data.matrix(expr[, row.names(a), drop=F])
  if (!is.null(rcFilter)) {
    keep <- apply(e>=rcFilter, 1, sum) >= atLeastIn
    e <- e[keep, , drop=F]
  }

  eset <- createSeqExpressionSet(e,a)
  eset <- EDASeq::betweenLaneNormalization(eset, which="upper")
  if (!is.null(unvariantNumber)) {
    unvariant <- computeEmpiricalUnvariants(eset, unvariantNumber, varName)
  } else {
    unvariant <- NULL
  }

  return(list(
    eset = eset,
    unvariant=unvariant,
    methods = paste0("normalization=upper\n",
                     "cellType=",cell,
                     "\ncondition=",condition,
                     "\ngenes with less then ", rcFilter, " in one sample were filtered out.")
  ))
}

#' Filter expression using RPKM
#'
#' @inheritParams filterExpression
#' @inheritParams computeRPKM
#' @param rpkmMean minimal RPKM row mean accepted to keep the row; NULL
#'
#' @return a filtered eset and other details
#'
#' @importFrom EDASeq betweenLaneNormalization
#' @rdname RUV-handling
#'
#' @export
#'
filterExpressionByRpkmMean <- function(cell, condition, expr, annotations, rcFilter=NULL,
                                       atLeastIn=1, rpkmMean=NULL, unvariantNumber=5000,
                                       varName="days", gene_length="gene_length.txt") {
  if (cell=="ALL") {
    sel <- annotations$condition %in% condition
  } else {
    sel <- annotations$cell %in% cell & annotations$condition %in% condition
  }

  a <- annotations[sel, , drop=F]
  e <- data.matrix(expr[, row.names(a), drop=F])

  eset <- createSeqExpressionSet(e,a)
  rpkm <- computeRPKM(eset, gene_length = gene_length)

  keepRpmk <- rep(TRUE, nrow(rpkm))
  if (!is.null(rpkmMean)) {
    keepRpmk <- rowMeans(rpkm) >= rpkmMean
  }
  usableGenes <- row.names(rpkm[keepRpmk, , drop=F])

  keep <- rep(TRUE, nrow(e))
  if (!is.null(rcFilter)) {
    keep <- apply(e>=rcFilter, 1, sum) >= atLeastIn
  }

  e <- e[keep, , drop=F]
  usableGenes <- intersect(row.names(e), usableGenes)
  e <- e[usableGenes, , drop=F]

  eset <- createSeqExpressionSet(e,a)
  eset <- betweenLaneNormalization(eset, which="upper")
  if (!is.null(unvariantNumber)) {
    unvariant <- computeEmpiricalUnvariants(eset, unvariantNumber, varName)
  } else {
    unvariant <- NULL
  }

  return(list(
    eset = eset,
    unvariant=unvariant,
    rpkm = rpkm,
    methods = paste0("normalization=upper\n",
                     "cellType=",cell,
                     "\ncondition=",condition,
                     "\ngenes with less then ", rcFilter, " in one sample were filtered out.")
  ))
}

#' Compute RPKM from eset
#'
#' @param eset expression set
#' @param gene_length filename of file containg gene lengths
#'
#' @importFrom EDASeq counts
#' @importFrom edgeR DGEList rpkm
#' @importFrom utils read.table
#'
#' @rdname RUV-handling
#'
#' @export
#'
computeRPKM <- function(eset, gene_length="../m38.p5-ensembl-91-renge-length/gene_length.txt") {
  lengths <- read.table(gene_length, header=T, row.names=1,
                        sep = "\t", quote = "\"", check.names = FALSE,
                        comment.char = "", stringsAsFactors = FALSE)
  usableGenes <- intersect(row.names(eset), lengths$gene)
  lengths <- lengths[lengths$gene %in% usableGenes, ]
  maxMeanGeneLength <- tapply(seq_len(NROW(lengths)), lengths$gene, function(idx) max(lengths$mean[idx]))

  expr <- EDASeq::counts(eset)[usableGenes, ]
  maxMeanGeneLength <- maxMeanGeneLength[usableGenes]
  y <- edgeR::DGEList(counts=expr)
  y$genes <- data.frame(Length=as.numeric(maxMeanGeneLength))
  edgeR::rpkm(y, gene_length="Length")
}

