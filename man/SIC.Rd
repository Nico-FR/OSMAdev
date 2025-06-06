% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIC.R, R/SIC2.R
\name{SIC}
\alias{SIC}
\title{Structural impact score between 2 matrices}
\usage{
SIC(
  mutated.mat,
  wildtype.mat,
  bin.width,
  vp.start,
  vp.stop = NULL,
  start = NULL,
  stop = NULL,
  verbose = TRUE
)

SIC(
  mutated.mat,
  wildtype.mat,
  bin.width,
  vp.start,
  vp.stop = NULL,
  start = NULL,
  stop = NULL,
  verbose = TRUE
)
}
\arguments{
\item{mutated.mat, wildtype.mat}{mutated (query) or wildtype (control, subject) matrices as \code{dgCMatrix} or \code{matrix} object for only one chromosome.}

\item{bin.width}{Bin width of the matrix in base pair.}

\item{vp.start}{Start of the view point in base pair.}

\item{vp.stop}{Stop/end of the view point in base pair. Default is NULL to only use the bin where vp.start is located (ie: vp.stop = vp.start + bin.width).}

\item{start, stop}{Area in bp of the chromosome to compute impact score. Default is NULL to use the entire chromosome (i.e. entire matrix).}

\item{verbose}{if TRUE, show information when view point does not match a bin annotations.}
}
\value{
integer

integer
}
\description{
Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between 2 matrices.
Fold change can be compute:
\itemize{
\item from one bin to all bins in the matrix (using view point parameters),
\item to one bin to all bins in the matrix (using start/stop parameters).
}

Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between 2 matrices.
Fold change is computed for each bin of the matrices provided
}
\examples{
# to do

# to do

}
