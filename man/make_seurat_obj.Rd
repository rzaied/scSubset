% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadData.R
\name{make_seurat_obj}
\alias{make_seurat_obj}
\title{Function to create seaurat object}
\usage{
make_seurat_obj(seurat.data, dataset_name)
}
\arguments{
\item{seurat.data}{Reactive value containing data from single cell dataset}

\item{dataset_name}{Reactive value containing name of project for seurat}
}
\value{
Returns a Reactive value containing list of downsampled seurat objects
 with reassigned cluster labels that correspond across downsampled subsets
}
\description{
Function to create seaurat object
}
