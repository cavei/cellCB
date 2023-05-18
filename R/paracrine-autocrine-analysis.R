#' get Cell Expression Signature
#'
#' @param RData an data RData
#' @param objname the loaded object name
#' @param cols columns to use
#'
#' @rdname paracrine-autocrine-functions
#'
#' @export
#'
getCellExpressionSignature <- function(RData="rpkm-Rdatas/all-wt.RData", objname="all_wt", cols=1:4) {
  load(RData)
  loadObj <- get(objname)

  lfcs <- loadObj$esetRUV$DEGsTable$table[,cols]
  sig <- apply(lfcs, 1, anyCmpSigByLogFC, cmp = colnames(lfcs), lfcThr = 1.5)

  tableWithLoosers <- lapply(sig, function(s) {
    if (!is.null(s))
      tapply(s$looser, s$winner, function(x) paste(x, collapse=";"))
  })

  tables <- lapply(sig, function(s) {table(s$winner)})

  cellExpressionTable <- meltTablesWithListNames(tables)
  cellExpressionTableWitLoosers <- meltTablesWithLoosers(tableWithLoosers)

  if (!identical(cellExpressionTable$gene, cellExpressionTableWitLoosers$gene))
    stop("Table names does not match")

  cellExpressionTable$looser <- cellExpressionTableWitLoosers$loosers
  # save(cellExpressionTable, file="rpkm-Rdatas/cellExpressionTable.RData")
  cellExpressionTable
}

#' Define gene behaviour in the cell
#'
#' @param activeAtDays list of gene active at days
#' @param activeAtDaysTot list of gene active in total tissue
#' @param constitutiveWide constitutive active genes
#'
#' @rdname paracrine-autocrine-functions
#'
#' @export
#'
extractGenesBehavious <- function(activeAtDays, activeAtDaysTot, constitutiveWide) {
  # activeAtDays <- filter_dailyList(activeAtDays, constitutive) ### Check This Carefully
  allGeneTimeDependentInBatch <- unique(do.call(c, activeAtDays))

  constitutive <- setdiff(constitutiveWide, allGeneTimeDependentInBatch) # genes that are time dipendent are not constitutive
  activeAtDaysTotCell <- filter_dailyList(activeAtDaysTot, constitutive) # from total muscle list I extracted those that are expressed by the cell

  allGeneTimeDependentInTotCell <- unique(do.call(c, activeAtDaysTotCell))
  constitutive <- setdiff(constitutive, allGeneTimeDependentInTotCell) # cell specific genes constitutively expressed genes
  list(activeAtDays=activeAtDays,
       activeAtDaysTotCell=activeAtDaysTotCell,
       constitutive=constitutive)
}

#' Define ligand receptor behaviour in the cell
#'
#' @param genesBehaviours list of gene active at days from extractGenesBehavious
#' @param db the ligand receptor database
#'
#' @rdname paracrine-autocrine-functions
#'
#' @export
#'
extractLigandReceptorBehaviours <- function(genesBehaviours, db) {
  activeAtDays <- genesBehaviours$activeAtDays
  activeAtDaysTotCell <- genesBehaviours$activeAtDaysTotCell
  constitutive <- genesBehaviours$constitutive

  ligandReceptorDayBatch <- lapply(activeAtDays, function(dgenes) {
    list(ligandOfTheDay = intersect(dgenes, db$ligands),
         receptorOfTheDay = intersect(dgenes, db$receptors))
  })

  ligandReceptorDayTot <- lapply(activeAtDaysTotCell, function(dgenes) {
    list(ligandOfTheDay = intersect(dgenes, db$ligands),
         receptorOfTheDay = intersect(dgenes, db$receptors))
  })

  ligandReceptorConstitutive <- list(ligandOfTheDay = intersect(constitutive, db$ligands),
                                     receptorOfTheDay = intersect(constitutive, db$receptors))
  list(ligandReceptorDayBatch=ligandReceptorDayBatch,
       ligandReceptorDayTot=ligandReceptorDayTot,
       ligandReceptorConstitutive=ligandReceptorConstitutive)
}

#' Assign color to the nodes
#'
#' @inheritParams wiseMergeTimePoints
#' @param paletteName the name of the palette to use.
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
assignColorCode <- function(ligandReceptorBehaviours, paletteName="Greens") {
  ligandReceptorDayBatch=ligandReceptorBehaviours$ligandReceptorDayBatch
  ligandReceptorDayTot=ligandReceptorBehaviours$ligandReceptorDayTot
  ligandReceptorConstitutive=ligandReceptorBehaviours$ligandReceptorConstitutive

  # library(RColorBrewer)
  ligandCols <- rev(RColorBrewer::brewer.pal(8, paletteName))[c(1,5,8)]
  receptorCols <- rev(RColorBrewer::brewer.pal(8, "Blues"))[c(1,5,8)]

  topLigand <- unique(unname(do.call(c, lapply(ligandReceptorDayBatch, function(x) x$ligandOfTheDay))))
  totLigand <- unique(unname(do.call(c, lapply(ligandReceptorDayTot, function(x) x$ligandOfTheDay))))
  ligandsCols <- c(rep(ligandCols[1], length(topLigand)),
                   rep(ligandCols[2], length(totLigand)))
  names(ligandsCols) <- c(topLigand, totLigand)

  topReceptor <- unique(unname(do.call(c, lapply(ligandReceptorDayBatch, function(x) x$receptorOfTheDay))))
  totReceptor <- unique(unname(do.call(c, lapply(ligandReceptorDayTot, function(x) x$receptorOfTheDay))))
  constitutiveReceptor <- ligandReceptorConstitutive$receptorOfTheDay
  receptorsCols <- c(rep(receptorCols[1], length(topReceptor)),
                     rep(receptorCols[2], length(totReceptor)),
                     rep(receptorCols[3], length(constitutiveReceptor)))
  names(receptorsCols) <- c(topReceptor, totReceptor, constitutiveReceptor)

  bothLandR <- intersect(names(ligandsCols), names(receptorsCols))
  if (length(bothLandR) > 0)
    warning(paste0("Some genes are both ligand and receptor. See ...\n\t",paste0(bothLandR, collapse = ", "), "\n"))
  list(ligandsCols=ligandsCols, receptorsCols=receptorsCols, bothLandR=bothLandR)
}

#' Assign color to the nodes
#'
#' @inheritParams wiseMergeTimePoints
#' @inheritParams renderHTMLvideo
#'
#' @rdname paracrine-autocrine-functions
#' @export
assignCellCode <- function(ligandReceptorBehaviours, cell) {
  ligandReceptorDayBatch=ligandReceptorBehaviours$ligandReceptorDayBatch
  ligandReceptorDayTot=ligandReceptorBehaviours$ligandReceptorDayTot
  ligandReceptorConstitutive=ligandReceptorBehaviours$ligandReceptorConstitutive

  topLigand <- unique(unname(do.call(c, lapply(ligandReceptorDayBatch, function(x) x$ligandOfTheDay))))
  totLigand <- unique(unname(do.call(c, lapply(ligandReceptorDayTot, function(x) x$ligandOfTheDay))))
  constitutiveLigand <- ligandReceptorConstitutive$ligandOfTheDay
  ligandsLabel <- c(rep(paste0("b_", cell), length(topLigand)),
                    rep(paste0("t_", cell), length(totLigand)),
                    rep(paste0("c_", cell), length(constitutiveLigand)))
  names(ligandsLabel) <- c(topLigand, totLigand, constitutiveLigand)

  topReceptor <- unique(unname(do.call(c, lapply(ligandReceptorDayBatch, function(x) x$receptorOfTheDay))))
  totReceptor <- unique(unname(do.call(c, lapply(ligandReceptorDayTot, function(x) x$receptorOfTheDay))))
  constitutiveReceptor <- ligandReceptorConstitutive$receptorOfTheDay

  receptorsLabel <- c(rep(paste0("b_",cell), length(topReceptor)),
                      rep(paste0("t_", cell), length(totReceptor)),
                      rep(paste0("c_", cell), length(constitutiveReceptor)))
  names(receptorsLabel) <- c(topReceptor, totReceptor, constitutiveReceptor)

  bothLandR <- intersect(names(ligandsLabel), names(receptorsLabel))
  if (length(bothLandR) > 0)
    warning(paste0("Some genes are both ligand and receptor. See ...\n\t",paste0(bothLandR, collapse = ", "), "\n"))
  list(ligandsLabels=ligandsLabel, receptorsLabels=receptorsLabel, bothLandR=bothLandR)
}

#' Assign color to the nodes
#'
#' @inheritParams renderHTMLvideo
#' @param ligands vector of ligands
#' @param receptors vector of receptors
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
resolveWarnings <- function(lrCols, ligands, receptors){
  receptorsCols = removeElement(lrCols$receptorsCols, ligands)
  ligandsCols = removeElement(lrCols$ligandsCols, receptors)
  bothLandR <- intersect(names(ligandsCols), names(receptorsCols))
  if (length(bothLandR) > 0)
    warning(paste0("Some genes are both ligand and receptor. See ...\n\t",paste0(bothLandR, collapse = ", "), "\n"))
  list(ligandsCols=ligandsCols, receptorsCols=receptorsCols, bothLandR=bothLandR)
}

#' Assign color to the nodes
#'
#' @param namedVector named vector
#' @param nms names
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
removeElement <- function(namedVector, nms) {
  namedVector[!names(namedVector) %in% nms]
}

#' Merge time points wisely
#'
#' @param ligandReceptorBehaviours list of ligand and receptor ativation
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
wiseMergeTimePoints <- function(ligandReceptorBehaviours) {
  ligandReceptorDayBatch=ligandReceptorBehaviours$ligandReceptorDayBatch
  ligandReceptorDayTot=ligandReceptorBehaviours$ligandReceptorDayTot
  ligandReceptorConstitutive=ligandReceptorBehaviours$ligandReceptorConstitutive

  ligandsReceptorCumulative <- lapply(names(ligandReceptorDayBatch),function(day) {
    list(ligandOfTheDay = c(ligandReceptorDayBatch[[day]]$ligandOfTheDay,
                            ligandReceptorDayTot[[day]]$ligandOfTheDay,
                            ligandReceptorConstitutive$ligandOfTheDay),
         receptorOfTheDay = c(ligandReceptorDayBatch[[day]]$receptorOfTheDay,
                              ligandReceptorDayTot[[day]]$receptorOfTheDay,
                              ligandReceptorConstitutive$receptorOfTheDay))
  })
  names(ligandsReceptorCumulative) <- names(ligandReceptorDayBatch)
  ligandsReceptorCumulative
}

#' Merge time points wisely
#'
#' @inheritParams paraWiseMultiMergeTimePoints
#' @param ligCell cell that produce the ligands
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
paraWiseMergeTimePoints <- function(ligandReceptorBehavioursCell, ligCell="fap", recCell="ec") {
  ligandProvidingCell.batch <- ligandReceptorBehavioursCell[[ligCell]]$ligandReceptorDayBatch
  ligandProvidingCell.tot <- ligandReceptorBehavioursCell[[ligCell]]$ligandReceptorDayTot
  ligandProvidingCell.const <- ligandReceptorBehavioursCell[[ligCell]]$ligandReceptorConstitutive

  receptorProvidingCell.batch <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorDayBatch
  receptorProvidingCell.tot <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorDayTot
  receptorProvidingCell.const <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorConstitutive

  days <- intersect(names(ligandProvidingCell.batch), names(receptorProvidingCell.batch))

  ligandsReceptorCumulative <- lapply(days,function(day) {
    list(ligandOfTheDay = c(ligandProvidingCell.batch[[day]]$ligandOfTheDay,
                            ligandProvidingCell.tot[[day]]$ligandOfTheDay),
         receptorOfTheDay = c(receptorProvidingCell.batch[[day]]$receptorOfTheDay,
                              receptorProvidingCell.tot[[day]]$receptorOfTheDay,
                              receptorProvidingCell.const$receptorOfTheDay))
  })
  names(ligandsReceptorCumulative) <- days
  ligandsReceptorCumulative
}

#' Merge time points wisely
#'
#' @param ligandReceptorBehavioursCell list of ligand enad receptor ativation
#' @param ligCells cells that produce the ligands
#' @param recCell cell that produce the receptors
#' @param removeAutocrine if TRUE remove aoutocrine loops
#' @param useConstitutiveLigands if TRUE use constitutive ligands
#'
#' @rdname paracrine-autocrine-functions
#' @export
#'
paraWiseMultiMergeTimePoints <- function(ligandReceptorBehavioursCell,
                                         ligCells=c("ec", "fap", "mp", "inf", "per"), recCell="ec",
                                         removeAutocrine=TRUE, useConstitutiveLigands=FALSE) {

  cells <- ligCells
  if (removeAutocrine)
    cells <- ligCells[!ligCells %in% recCell]

  ligandProvidingCell.batch <- lapply(cells, function(cell) ligandReceptorBehavioursCell[[cell]]$ligandReceptorDayBatch)
  names(ligandProvidingCell.batch) <- cells
  ligandProvidingCell.tot <-   lapply(cells, function(cell) ligandReceptorBehavioursCell[[cell]]$ligandReceptorDayTot)
  names(ligandProvidingCell.tot) <- cells
  ligandProvidingCell.const <-  lapply(cells, function(cell) ligandReceptorBehavioursCell[[cell]]$ligandReceptorConstitutive)
  names(ligandProvidingCell.const) <- cells

  receptorProvidingCell.batch <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorDayBatch
  receptorProvidingCell.tot <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorDayTot
  receptorProvidingCell.const <- ligandReceptorBehavioursCell[[recCell]]$ligandReceptorConstitutive

  ligandDays <- unique(unlist(lapply(ligandProvidingCell.batch, names)))
  days <- intersect(ligandDays, names(receptorProvidingCell.batch))

  ligandsReceptorCumulative <- lapply(days,function(day) {
    ligandOfTheDay.batch <- unique(unlist(lapply(ligandProvidingCell.batch, function(cellLig) cellLig[[day]]$ligandOfTheDay)))
    ligandOfTheDay.tot   <- unique(unlist(lapply(ligandProvidingCell.tot, function(cellLig) cellLig[[day]]$ligandOfTheDay)))

    ligandOfTheDay = c(ligandOfTheDay.batch, ligandOfTheDay.tot)

    if (useConstitutiveLigands){
      ligandFromCell.constitutive <- unique(unlist(lapply(ligandProvidingCell.const, function(x) x$ligandOfTheDay)))
      ligandOfTheDay = c(ligandOfTheDay, ligandFromCell.constitutive)
    }

    receptorOfTheDay = c(receptorProvidingCell.batch[[day]]$receptorOfTheDay,
                         receptorProvidingCell.tot[[day]]$receptorOfTheDay,
                         receptorProvidingCell.const$receptorOfTheDay)

    list(ligandOfTheDay = ligandOfTheDay, receptorOfTheDay = receptorOfTheDay)
  })
  names(ligandsReceptorCumulative) <- days
  ligandsReceptorCumulative
}

#' Render Html Video
#'
#' @param ligandsReceptorCumulative list of ligand receptors per day
#' @param receptorGraph the graph receptor
#' @param lrCols dictionary of ligand receptor colors
#' @param lrCell dictionary of ligand receptor cells
#' @param cell cell name
#' @param condition cell condition
#'
#' @rdname paracrine-autocrine-functions
#' @importFrom network get.vertex.attribute network.vertex.names %v% %v%<-
#' @export
#'
renderHTMLvideo <- function(ligandsReceptorCumulative, receptorGraph, lrCols, lrCell=NULL, cell, condition) {
  # @importFrom ndtv render.d3movie
  # @importFrom networkDynamic networkDynamic

  removeTimes <- NULL
  for (x in names(ligandsReceptorCumulative)) {
    l <- unlist(lapply(ligandsReceptorCumulative[[x]], length))
    if (any(l==0))
      removeTimes <- c(removeTimes, x)
  }
  if (!is.null(removeTimes))
    ligandsReceptorCumulative[[removeTimes]] <- NULL

  networklist <- lapply(ligandsReceptorCumulative, extractDayNetwork, graphNEL=receptorGraph)
  tnet<-networkDynamic::networkDynamic(network.list=networklist, vertex.pid="vertex.pid")

  # timeAutocrine <- lapply(ligandsReceptorCumulative, extractDayLigandReceptorNetwork,
  #                         graphNEL=receptorGraph)
  #
  # active <- lapply(names(timeAutocrine), function(day) {
  #   mat <- cbind(timeAutocrine[[day]]$activeEdges, day)
  #   data.frame(src=mat[,1], dest=mat[,2], day=mat[,3], time=as.numeric(sub("d", "", mat[,3])), stringsAsFactors = F)
  # })
  # df <- do.call(rbind, active)
  # network::get.vertex.attribute(tnet, "vertex.pid")

  nodeNames <- network::get.vertex.attribute(tnet, "vertex.pid")
  network::network.vertex.names(tnet) <- nodeNames

  # ligs <- unique(df$src)
  # colorsLig <- ligandsCols[ligs]
  #
  # recs <- unique(df$dest)
  # colorsRec <- receptorsCols[recs]
  colors <- c(lrCols$ligandsCols, lrCols$receptorsCols)
  tnet%v%"color" <- colors[nodeNames]

  if (!is.null(lrCell)){
    attrib <- c(lrCell$ligandsNames, lrCell$receptorsNames)
    tnet%v%"cells" <- attrib[nodeNames]
  }



  vertex.tooltip = function(slice){paste("<b>Gene:</b>", (slice%v%"vertex.names"), "<br>",
                                         "<b>Cells:</b>", (slice%v%"cells"))}

  ndtv::render.d3movie(tnet, displaylabels=TRUE, vertex.col="color",  vertex.tooltip= vertex.tooltip, output.mode = 'HTML', main=paste(cell, condition))
}

#' Render the video
#' @inheritParams renderHTMLvideo
#'
#' @rdname paracrine-autocrine-functions
#' @importFrom network get.vertex.attribute network.vertex.names %v% %v%<-

#' @export
#'
renderHTMLvideo2.0 <- function(ligandsReceptorCumulative, receptorGraph, cell, condition) {
# import from ndtv render.d3movie
# @importFrom networkDynamic networkDynamic
  removeTimes <- NULL
  for (x in names(ligandsReceptorCumulative)) {
    l <- unlist(lapply(ligandsReceptorCumulative[[x]], length))
    if (any(l==0))
      removeTimes <- c(removeTimes, x)
  }
  if (!is.null(removeTimes))
    ligandsReceptorCumulative[[removeTimes]] <- NULL

  networklist <- lapply(ligandsReceptorCumulative, extractDayNetwork, graphNEL=receptorGraph)
  tnet<-networkDynamic::networkDynamic(network.list=networklist, vertex.pid="vertex.pid")

  nodeNames <- network::get.vertex.attribute(tnet, "vertex.pid")
  network::network.vertex.names(tnet) <- nodeNames

  # vertex.tooltip = function(slice){paste("<b>Gene:</b>", (slice%v%"vertex.names"), "<br>",
  #                                        "<b>Cells:</b>", (slice%v%"cells"))}
  vertex.tooltip = function(slice){paste("<b>Gene:</b>", (slice%v%"vertex.names"))}

  ndtv::render.d3movie(tnet, displaylabels=TRUE,
                       # vertex.col="color",
                       vertex.tooltip= vertex.tooltip,
                       output.mode = 'HTML', main=paste(cell, condition))
}
