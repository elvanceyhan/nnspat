% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{Znnsym2cl}
\alias{Znnsym2cl}
\title{NN Symmetry Test with Normal Approximation for Two Classes}
\usage{
Znnsym2cl(
  dat,
  lab,
  type = "dixon",
  alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95
)
}
\arguments{
\item{dat}{The data set in one or higher dimensions, each row corresponds to a data point.}

\item{lab}{The \code{vector} of class labels (numerical or categorical)}

\item{type}{The type of the NN symmetry test with default=\code{"dixon"}.
Takes on values \code{"dixon"} and \code{"pielou"} for Dixon's and Pielou's (first type) NN symmetry test}

\item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.}

\item{conf.level}{Level of the upper and lower confidence limits, default is \code{0.95}, 
for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}}
}
\value{
A \code{list} with the elements
\item{statistic}{The \eqn{Z} test statistic for Pielou's first type of NN symmetry test}
\item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
\item{conf.int}{Confidence interval for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}} in the \eqn{2 \times 2} NNCT
at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
\item{estimate}{Estimate, i.e., the difference of the off-diagonal entries of the \eqn{2 \times 2} NNCT, \eqn{N_{12}-N_{21}}.}
\item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
\eqn{E(N_{12})-E(N_{21})} in the \eqn{2 \times 2} NNCT, which is 0 for this function.}
\item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
\item{method}{Description of the hypothesis test}
\item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
}
\description{
An object of class \code{"htest"} performing hypothesis test of equality of the expected value of the off-diagonal 
cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k=2} classes.
That is, the test performs Dixon's or Pielou's (first type of) NN symmetry test which is appropriate 
(i.e. have the appropriate asymptotic sampling distribution)
for completely mapped data and for sparsely sample data, respectively.
(See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).

The symmetry test is based on the normal approximation of the difference of the off-diagonal entries
in the NNCT and are due to \insertCite{pielou:1961,dixon:1994;textual}{nnspat}.

The \code{type="dixon"} refers to Dixon's NN symmetry test and 
\code{type="pielou"} refers to Pielou's first type of NN symmetry test.

The function yields the test statistic, \eqn{p}-value for the
corresponding alternative, the confidence interval, estimate and null value for the parameter of interest
(which is the difference of the off-diagonal entries in the NNCT), and method and name of the data set used.

The null hypothesis is that all \eqn{E(N_{12})=E(N_{21})} in the \eqn{2 \times 2} NNCT (i.e., symmetry in the 
mixed NN structure).

See also
(\insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
and the references therein.
}
\examples{
n<-20  #or try sample(1:20,1)
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))

Znnsym2cl(Y,cls)
Znnsym2cl(Y,cls,type="pielou")

Znnsym2cl(Y,cls,alt="g")
Znnsym2cl(Y,cls,type="pielou",alt="g")

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{Znnsym2cl.ss.ct}}, \code{\link{Znnsym2cl.ss}}, \code{\link{Znnsym2cl.dx.ct}},
\code{\link{Znnsym2cl.dx}}, \code{\link{Znnsym.ss.ct}}, \code{\link{Znnsym.ss}}, \code{\link{Znnsym.dx.ct}},
\code{\link{Znnsym.dx}}, \code{\link{Znnsym.dx.ct}}, \code{\link{Znnsym.dx}} and \code{\link{Znnsym}}
}
\author{
Elvan Ceyhan
}
