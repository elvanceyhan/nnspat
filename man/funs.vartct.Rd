% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{funs.vartct}
\alias{funs.vartct}
\alias{var.tctI}
\alias{var.tctIII}
\alias{var.tctIV}
\title{Functions for Variances of Cell Counts in the Types I, III and IV TCTs}
\usage{
var.tctI(ct, covN)

var.tctIII(ct, covN)

var.tctIV(ct, covN)
}
\arguments{
\item{ct}{A nearest neighbor contingency table}

\item{covN}{The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.}
}
\value{
Each of these functions returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the variances of
the entries in the TCT for the corresponding type of cell-specific test.
The row and column names are inherited from \code{ct}.
}
\description{
Three functions: \code{var.tctI}, \code{var.tctIII} and \code{var.tctIV}.

These functions return the variances of \eqn{T_{ij}} values for \eqn{i,j=1,\ldots,k} in the TCT in matrix form which
is of the same dimension as TCT for types I, III and IV tests. 
The argument \code{covN} must be the covariance between \eqn{N_{ij}} values which are obtained from the NNCT by row-wise
vectorization.
These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.

See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{var.tct}} and \code{\link{var.nnct}}
}
