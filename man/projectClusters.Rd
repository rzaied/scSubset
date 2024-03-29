% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoClustering.R
\name{projectClusters}
\alias{projectClusters}
\title{Function to project downsampled subsets against the parent dataset and score projection quality}
\usage{
projectClusters(seuratObjectsList)
}
\arguments{
\item{seuratObjectsList}{Reactive value containg list of downsampled seurat objects
with reduced dimensions (PCA data), scaled counts, and cluster labels that corresponds across subsets}
}
\value{
Returns a bar plot showing projection quality scores across subsets
}
\description{
Function to project downsampled subsets against the parent dataset and score projection quality
}
