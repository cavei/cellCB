% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additionalTCA.R
\name{clusterWithTca}
\alias{clusterWithTca}
\alias{clusterTimes}
\alias{watchCluster}
\title{Perform TCA based clustering}
\usage{
clusterWithTca(esetRUV, genes, clusterK = 9)

clusterTimes(
  normMatrix,
  days,
  genes,
  algo = "km",
  clusterK = 9,
  takeLog = FALSE
)

watchCluster(tca, cols = 3)
}
\arguments{
\item{esetRUV}{expression set RUV}

\item{genes}{gene names}

\item{clusterK}{number of clusters K}

\item{normMatrix}{normalized matrix}

\item{days}{days}

\item{algo}{algorithm}

\item{takeLog}{logical; need log?}

\item{tca}{the results of clusterTimes function}

\item{cols}{subdivide the column in _col_ cols}
}
\value{
list with all you need

TCA object
}
\description{
Perform TCA based clustering

Perform TCA based clustering

Visualize the cluster produced
}
