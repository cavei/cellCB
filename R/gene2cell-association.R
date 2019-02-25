associateGeneToCells <- function(cellExpressionTable, winCmpToBeExclusive=4) {
  selEx <- cellExpressionTable$cmpWin >= winCmpToBeExclusive
  exclusive <- cellExpressionTable[selEx, ,drop=F]
  pleiotropic <- cellExpressionTable[!selEx, ,drop=F]
  selPl <- pleiotropic$gene %in% exclusive$gene
  pleioExcluded <- pleiotropic[selPl, , drop=F]
  pleiotropic <- pleiotropic[!selPl, , drop=F]
  list(cellExclusive=exclusive, cellPleiotropic=pleiotropic, excluded=pleioExcluded)
}

tables2maps <- function(geneCellTable) {
  full <- rbind(geneCellTable$cellExclusive, geneCellTable$cellPleiotropic)[, c("gene", "cell")]
  tapply(seq_along(full$gene), as.factor(full$cell), function(idx) full$gene[idx])
}

filter_dailyList <- function(l, genes) {
  lapply(l , function(x) {
    intersect(x, genes)
  })
}
