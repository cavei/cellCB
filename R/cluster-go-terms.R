#' Filter significance table pickeing 1 term for cluster
#'
#' @param sigTable the table with significant results
#' @param clusters the GO organized in clusters
#' @param lowerCut excludes the GO with fewer than lowerCut significant times
#' @param upperCut excludes the GO with more than upperCut significant times
#'
#' @export
#'
sigTable_filter_by_clusters <- function(sigTable, clusters, lowerCut = 2, upperCut=Inf) {
  terms <- exclude_go_terms_by_counts(sigTable, lowerCut=lowerCut, upperCut=upperCut)
  if (any(is.na(clusters[terms])))
    stop("Some terms are not present in clusters")

  cls <- clusters[terms]
  keep_go <- unname(guess_best_sets(cls, sigTable))
  list(keep_go=keep_go, grp=groupByCluster(cls))
}


remove_words <- function(x) {
  patterns = c("positive regulation of ", "negative regulation of ",
               " signaling pathway", " signal transduction", " system development",
               "activation of ", "regulation of ",
               "cellular response to ", "response to ",
               "system", "development ", " process", " stimulus", " protein",
               "branching", "involved ")

  sapply(x, function(string) {
    out <- string
    for (ptt in patterns){
      out <- gsub(ptt, "", out)
    }
    out
  }
  )
}

#' @importFrom stats as.dist cutree hclust
#' @importFrom graphics abline plot
#' @importFrom utils adist
cluster.goTerms <- function(terms, h=NULL, k=NULL) {
  # library("DescTools")
  semanticDist <- adist(terms)

  labels <- remove_words(colnames(semanticDist))
  normalizer <- matrix(1, ncol=length(labels), nrow=length(labels),
                       dimnames = list(labels, labels))
  coordinates <- which(upper.tri(semanticDist), arr.ind = T)
  for (i in seq_len(NROW(coordinates))) {
    x <- coordinates[i,1]
    y <- coordinates[i,2]
    # nf <- abs(nchar(labels[x])-nchar(labels[y]))
    nf <- max(nchar(labels[x]),nchar(labels[y]))
    # normalizer[x, y] <- ifelse(nf!=0, nf, 1)
    normalizer[x, y] <- nf
    normalizer[y, x] <- normalizer[x, y]
  }

  hc <- hclust(as.dist(semanticDist/normalizer), method = "complete")
  plot(hc)

  if (!is.null(h)){
    abline(h=h, col="red")
    return(cutree(hc, h=h))
  }

  if (!is.null(k))
    return(cutree(hc, k=k))
}

levSim <- function(pure_terms){
  warning("
library(DescTools)
DescTools::StrDist(pure_terms, method=\"normlevenshtein\")
          ")
  # library("DescTools")
  # DescTools::StrDist(pure_terms, method="normlevenshtein")
}

groupByCluster <- function(cls) {tapply(names(cls), cls,function(x) x)}

#' @importFrom utils head
guess_best_sets <- function(clusters, df) {
  clsGrp <- groupByCluster(clusters)
  sapply(clsGrp, function(descr) {
    clsdf <- df[df$Description %in% descr, , drop=F]
    bestIdx <- head(order(clsdf$"p.adjust"), 1)
    clsdf$"Description"[bestIdx]
  })
}

exclude_go_terms_by_counts <- function(df, lowerCut=1, upperCut=Inf) {
  names(which(table(df$Description) >= lowerCut & table(df$Description) <= upperCut))
}

# exclude_class_go_term_by_counts <- function(df, class, lowerCut=1, upperCut=Inf) {
#   names(which(table(df$Description) >= lowerCut & table(df$Description) <= upperCut))
# }
