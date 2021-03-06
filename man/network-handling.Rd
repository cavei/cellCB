% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network-related-functions.R
\name{extractDayLigandReceptorNetwork}
\alias{extractDayLigandReceptorNetwork}
\alias{extractDayNetwork}
\title{Given a list of ligand receptor it extracts the network}
\usage{
extractDayLigandReceptorNetwork(lrday, graphNEL)

extractDayNetwork(lrday, graphNEL)
}
\arguments{
\item{lrday}{a list with ligand and receptors}

\item{graphNEL}{a graphNEL object representing the network}
}
\value{
a list
   \item{g}{the graph}
   \item{activeEdges}{the active edges}
   \item{autocrine}{connected nodes}

a network
}
\description{
Given a list of ligand receptor it extracts the network

Given a list of ligand receptor it extracts the network
}
