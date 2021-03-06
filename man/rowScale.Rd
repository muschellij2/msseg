% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rowScale.R
\name{rowScale}
\alias{rowScale}
\title{Row Centering and Scaling}
\usage{
rowScale(x, center = TRUE, scale = TRUE, add_attr = TRUE, rows = NULL,
  cols = NULL, na.rm = TRUE, ...)
}
\arguments{
\item{x}{numeric matrix}

\item{center}{Should the rows be centered}

\item{scale}{Should the rows be scaled}

\item{add_attr}{Add the center/scale attributes}

\item{rows}{A vector indicating subset of rows to operate over.
If \code{NULL}, no subsetting is done.}

\item{cols}{A vector indicating subset of columns to operate over.
If \code{NULL}, no subsetting is done.}

\item{na.rm}{If \code{TRUE}, missing values are removed first, otherwise not.}

\item{...}{Arguments to pass to \code{\link{rowSds}}}
}
\value{
Matrix of centered/scaled values
}
\description{
This function is a fast way to scale rows of a matrix
}
\examples{
x = matrix(rnorm(10*1000, mean = 4, sd = 5), ncol = 10)
cx = rowScale(x)
all(abs(rowMeans(cx)) < 1e-8)
all(abs(matrixStats::rowSds(cx) - 1) < 1e-8)
}
