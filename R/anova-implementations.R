#' ANOVA like analysis (edgeR) with RUV on times
#'
#' @param eset expression set
#' @param unvariant set of unvariant genes
#' @param k number of covariate
#' @param formula a string that can be converted in formula
#' @param varName the name of the Pdata variable to use for the anova
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp glmFit estimateDisp topTags
#' @importFrom RUVSeq RUVg
#' @importFrom stats model.matrix
#' @importFrom Biobase pData
#' @importFrom EDASeq counts
#'
#' @return a list with top Tags table and setEmpirical
#' @export
#'
edgeR.anovaOnTimeWithRUV <- function(eset, unvariant, k=1, formula="~days + W_1", varName="days") {

  setEmpirical <- RUVSeq::RUVg(eset, unvariant, k=k)

  design <- model.matrix(as.formula(formula), data=pData(setEmpirical))
  y <- edgeR::DGEList(counts=EDASeq::counts(setEmpirical), group=pData(setEmpirical)[, varName])
  y <- edgeR::calcNormFactors(y, method="upperquartile")
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  if (is.na(y$common.dispersion) & all(is.na(y$tagwise.dispersion))) {
    y <- edgeR::DGEList(counts=EDASeq::counts(setEmpirical), group=pData(setEmpirical)[, varName])
    y <- edgeR::calcNormFactors(y, method="upperquartile")
    y1 <- y
    y1$samples$group <- 1
    y0 <- edgeR::estimateDisp(y1[unvariant,], trend="none", tagwise=FALSE)
    y$common.dispersion <- y0$common.dispersion
  }
  fit <- edgeR::glmFit(y, design)
  coefs <- grep(c("W_"), colnames(design), invert = T)[-1]
  lrt <- edgeR::glmLRT(fit, coef=coefs)
  list(DEGsTable=edgeR::topTags(lrt, n=NROW(setEmpirical)), setEmpirical=setEmpirical)
}

#' ANOVA like analysis (edgeR)
#'
#' @inheritParams edgeR.anovaOnTimeWithRUV
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp glmFit estimateDisp topTags
#' @importFrom RUVSeq RUVg
#' @importFrom stats model.matrix
#' @importFrom EDASeq counts
#'
#' @export
#'
edgeR.anovaOnTime <- function(eset, k=1, formula="~days", varName="days") {
  # setEmpirical <- RUVg(eset, unvariant, k=k)
  setEmpirical <- eset
  design <- model.matrix(as.formula(formula), data=pData(setEmpirical))
  y <- edgeR::DGEList(counts=EDASeq::counts(setEmpirical), group=pData(setEmpirical)[, varName])
  y <- edgeR::calcNormFactors(y, method="upperquartile")
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  if (is.na(y$common.dispersion) & all(is.na(y$tagwise.dispersion))) {
    y <- edgeR::DGEList(counts=EDASeq::counts(setEmpirical), group=pData(setEmpirical)[, varName])
    y <- edgeR::calcNormFactors(y, method="upperquartile")
    y1 <- y
    y1$samples$group <- 1
    # y0 <- edgeR::estimateDisp(y1[unvariant,], trend="none", tagwise=FALSE)
    warning("Not sure about this thing. Check carefully line 65 anova-implementation.R")
    y0 <- edgeR::estimateDisp(y1, trend="none", tagwise=FALSE)
    y$common.dispersion <- y0$common.dispersion
  }
  fit <- edgeR::glmFit(y, design)
  coefs <- grep(c("W_"), colnames(design), invert = T)[-1]
  lrt <- edgeR::glmLRT(fit, coef=coefs)
  list(DEGsTable=edgeR::topTags(lrt, n=NROW(setEmpirical)), setEmpirical=setEmpirical)
}

#' ANOVA like analysis (edgeR)
#'
#' @inheritParams edgeR.anovaOnTimeWithRUV
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp glmFit estimateDisp topTags
#' @importFrom RUVSeq RUVg
#' @importFrom stats model.matrix
#' @importFrom Biobase pData
#' @importFrom EDASeq counts
#'
#' @export
#'
edgeR.allCmpOnTimeWithRUV <- function(eset, unvariant, k=1, formula="~days + W_1", varName="days") {
  setEmpirical <- RUVg(eset, unvariant, k=k)

  design <- model.matrix(as.formula(formula), data=pData(setEmpirical))
  y <- DGEList(counts=EDASeq::counts(setEmpirical), group=pData(setEmpirical)[, varName])
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  coefs <- grep(c("W_"), colnames(design), invert = T)[-1]
  DEGS <- lapply(coefs, function(coef) {
    edgeR::topTags(glmLRT(fit, coef=coef), n=NROW(setEmpirical))
  })
  list(DEGsTable=DEGS, setEmpirical=setEmpirical)
}
