
create.GO_DB <- function(ont="BP", keyType = "SYMBOL", OrgDb="org.Mm.eg.db"){
  requireNamespace(OrgDb)
  OrgDb=get("org.Mm.eg.db")
  clusterProfiler:::get_GO_data(OrgDb, ont, keyType)
}

#' Perform enricher with a slight change
#'
#' @param gene gene
#' @param go_bp_full_data go biological processed
#' @param pvalueCutoff pvalue cut off
#' @param pAdjustMethod BH
#' @param universe the full list of genes
#' @param qvalueCutoff 0.2
#' @param minGSSize 10
#' @param maxGSSize 500
#' @param readable FALSE
#' @param pool FALSE
#' @param ont BP
#' @param keyType SYMBOL
#' @param OrgDb org.Mm.eg.db
#'
#' @importFrom DOSE setReadable
#'
#' @export
#'
my.enricher <- function (gene, go_bp_full_data, pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", universe, qvalueCutoff = 0.2, minGSSize = 10,
                         maxGSSize = 500, readable = FALSE, pool = FALSE, ont="BP",
                         keyType = "SYMBOL", OrgDb="org.Mm.eg.db") {
  if (missing(universe))
    universe <- NULL

  requireNamespace("org.Mm.eg.db")
  OrgDb=get("org.Mm.eg.db")

  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = go_bp_full_data)
  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(OrgDb)
  if (readable) {
    res <- DOSE::setReadable(res, OrgDb)
  }
  res@ontology <- ont
  return(res)
}

#' Perform compare cluster
#'
#' @param geneClusters gene
#' @param go_bp_full_data go biological processed
#' @param data NULL
#' @param \dots other argument
#'
#' @importFrom magrittr %>%
#' @importFrom stats formula
#' @importFrom plyr dlply llply ldply rename .
#'
#' @export
#'
my.compareCluster <- function (geneClusters, go_bp_full_data, data = "", ...) {
  # library(magrittr)
  # library(plyr)
  fun_name <- "my.enricher"
  fun <- my.enricher
  if (typeof(geneClusters) == "language") {
    if (!is.data.frame(data)) {
      stop("no data provided with formula for compareCluster")
    }
    else {
      genes.var = all.vars(geneClusters)[1]
      grouping.formula = gsub("^.*~", "~", as.character(as.expression(geneClusters)))
      geneClusters = dlply(.data = data, formula(grouping.formula),
                           .fun = function(x) {
                             as.character(x[[genes.var]])
                           })
    }
  }
  clProf <- llply(geneClusters, .fun = function(i) {
    x = suppressMessages(fun(i, go_bp_full_data=go_bp_full_data, ...))
    if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
      as.data.frame(x)
    }
  })
  clusters.levels = names(geneClusters)
  clProf.df <- ldply(clProf, rbind)
  if (nrow(clProf.df) == 0) {
    stop("No enrichment found in any of gene cluster, please check your input...")
  }
  clProf.df <- rename(clProf.df, c(.id = "Cluster"))
  clProf.df$Cluster = factor(clProf.df$Cluster, levels = clusters.levels)
  if (is.data.frame(data) && grepl("+", grouping.formula)) {
    groupVarName <- strsplit(grouping.formula, split = "\\+") %>% unlist %>% gsub("~", "", .) %>% gsub("^\\s*", "",.) %>% gsub("\\s*$", "", .)
    groupVars <- sapply(as.character(clProf.df$Cluster), strsplit, split = "\\.") %>% do.call(rbind, .)
    for (i in seq_along(groupVarName)) {
      clProf.df[, groupVarName[i]] <- groupVars[, i]
    }
    i <- which(colnames(clProf.df) %in% groupVarName)
    j <- (1:ncol(clProf.df))[-c(1, i)]
    clProf.df <- clProf.df[, c(1, i, j)]
  }
  return(list(df = clProf.df, geneClusters = geneClusters,
              fun = fun_name, .call = match.call(expand.dots = TRUE)))
}
