#' @importFrom ggplot2 fortify ggplot aes_string geom_point scale_color_continuous
#'   guide_colorbar ylab ggtitle scale_size
#' @importFrom DOSE theme_dose
dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=10, keep="all", split = NULL,
                             font.size=12, title = "", orderBy="GeneRatio",
                             decreasing=TRUE) {

  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    size <- "Count"
  } else {
    message("invalid x, setting to 'GeneRatio' by default")
    x <- "GeneRatio"
    size <- "Count"
  }

  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)

  if (!orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "GeneRatio"`')
    orderBy <- "GeneRatio"
  }

  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) + scale_size(range=c(3, 8))

}


#####################
personal.dotplot.compareClusterResult <- function(object, x=~Cluster, colorBy="p.adjust", showCategory=5,
                                                  by="geneRatio", split=NULL, includeAll=TRUE,
                                                  font.size=12, title="",
                                                  filterSets=NULL, relabel_description=NULL,
                                                  relabel_names = NULL,
                                                  name_orders = NULL,
                                                  output_df=FALSE) {

  df <- personal.fortify.compareClusterResult(object, showCategory=showCategory, by=by,
                                              includeAll=includeAll, split=split, filterSets=filterSets)

  if (!is.null(relabel_names) & !(is.null(relabel_description)))
    stop("relabel conflict. Specify one of relabel_names or relabel description")

  if (!is.null(relabel_description)){
    df <- rename_description(df, relabel_description)
  }

  if (!is.null(relabel_names)){
    df <- rename_description_map(df, relabel_names)
  }

  #### THIS code need to be carefully reviewed becase of the sorting
  if (!is.null(name_orders)) {
    if (!all(unique(as.character(df$Description)) %in% name_orders)) {
      missing <- setdiff(unique(as.character(df$Description)), name_orders)
      stop(paste0("the following category are missin from imput: ",
                  paste(missing, collapse = ", ")))
    }
    sorted_descr <- name_orders
  } else {
    sorted_descr <- as.character(df$Description)[order(df$Cluster,-df$GeneRatio)]
  }
  df$Description <- factor(df$Description, levels = rev(unique(sorted_descr)))

  if (!is.null(name_orders)) {
    df$Description <- factor(df$Description, levels = rev(unique(sorted_descr)))
  }

  if (output_df)
    return(df)

  # enrichplot:::plotting.clusterProfile(df, x=x, type="dot", colorBy=colorBy, by=by, title=title, font.size=font.size)
  p <- ggplot(df, aes_string(x = "Cluster", y = "Description", size = "GeneRatio")) +
    geom_point(aes_string(color = colorBy))

  p + scale_color_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
    ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
    scale_size_continuous(range = c(3, 8)) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

}

rename_description <- function(df, relabel_description) {
  dict <- meltListOfVectors(relabel_description)
  new_names <- dict$nm
  names(new_names) <- dict$content
  df$Description <- factor(sapply(as.character(df$Description), function(d) new_names[d]))
  df
}


rename_description_map <- function(df, map) {
  df$Description <- factor(sapply(as.character(df$Description), function(d) map[d]))
  df
}


#' @importFrom ggplot2 fortify
#' @importFrom plyr ddply
#' @importFrom plyr mdply
#' @importFrom plyr .
personal.fortify.compareClusterResult <- function(model, data, showCategory=5, by="geneRatio",
                                                  split=NULL, includeAll=TRUE, filterSets=NULL) {
  # library("plyr")
  # library("ggplot2")
  clProf.df <- as.data.frame(model)
  if (!(is.null(filterSets))){
    clProf.df <- clProf.df[clProf.df$Description %in% filterSets, , drop=F]
  }

  .split <- split

  ## get top 5 (default) categories of each gene cluster.
  if (is.null(showCategory)) {
    result <- clProf.df
  } else {
    Cluster <- NULL # to satisfy codetools

    topN <- function(res, showCategory) {
      plyr::ddply(.data = res,
            .variables = .(Cluster),
            .fun = function(df, N) {
              if (length(df$Count) > N) {
                if (any(colnames(df) == "pvalue")) {
                  idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                } else {
                  ## for groupGO
                  idx <- order(df$Count, decreasing=T)[1:N]
                }
                return(df[idx,])
              } else {
                return(df)
              }
            },
            N=showCategory
      )

    }

    if (!is.null(.split) && .split %in% colnames(clProf.df)) {
      lres <- split(clProf.df, as.character(clProf.df[, .split]))
      lres <- lapply(lres, topN, showCategory = showCategory)
      result <- do.call('rbind', lres)
    } else {
      result <- topN(clProf.df, showCategory)
    }

  }

  ID <- NULL
  if (includeAll == TRUE) {
    result = subset(clProf.df, ID %in% result$ID)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$Count != 0, ]
  result$Description <- factor(result$Description,
                               levels=rev(GOlevel[,2]))


  if (by=="rowPercentage") {
    Description <- Count <- NULL # to satisfy codetools
    result <- ddply(result,
                    .(Description),
                    transform,
                    Percentage = Count/sum(Count),
                    Total = sum(Count))

    ## label GO Description with gene counts.
    x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
    y <- sapply(x[,3], paste, ")", sep="")
    result$Description <- y

    ## restore the original order of GO Description
    xx <- result[,c(2,3)]
    xx <- unique(xx)
    rownames(xx) <- xx[,1]
    Termlevel <- xx[as.character(GOlevel[,1]),2]

    ##drop the *Total* column
    result <- result[, colnames(result) != "Total"]

    result$Description <- factor(result$Description,
                                 levels=rev(Termlevel))

  } else if (by == "count") {
    ## nothing
  } else if (by == "geneRatio") {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
    result$GeneRatio = gsize/gcsize
    cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
    lv <- unique(cluster)[order(as.numeric(unique(result$Cluster)))]
    result$Cluster <- factor(cluster, levels = lv)
  } else {
    ## nothing
  }
  return(result)
}

