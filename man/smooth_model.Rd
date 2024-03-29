% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth_model.R
\name{smooth_model}
\alias{smooth_model}
\title{Main smoothLRC algorithm}
\usage{
smooth_model(x, u_init, v_init, w, index, lambda, epsilon, maxiter)
}
\arguments{
\item{x}{sparse matrix; assay matrix}

\item{u_init}{matrix; u initialization matrix}

\item{v_init}{matrix; v initialization matrix}

\item{w}{matrix; distance matrix. If null, computed to specification with k nearest neighbors.}

\item{index}{w indicies.}

\item{lambda}{positive numeric; penalization parameter.}

\item{epsilon}{positive numeric; convergence criterion.}

\item{maxiter}{positive integer; maximum desired iterations}
}
\value{
SummarizedExperiment object with u, v and cluster labels.
}
\description{
Main smoothLRC algorithm
}
\examples{

sce <- example_sce()
x <- SummarizedExperiment::assay(sce)
coords <- SummarizedExperiment::colData(sce)[, c("col", "row")]
k <- 10
init <- smooth_init(x, coords, k)
lambda <- 1
epsilon <- 1e-3
maxiter <- 5
smooth_model(x, init$u0, init$v0, init$w, init$index, lambda, epsilon, maxiter)

}
