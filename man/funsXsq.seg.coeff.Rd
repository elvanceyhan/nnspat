% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{funsXsq.seg.coeff}
\alias{funsXsq.seg.coeff}
\alias{Xsq.seg.coeff.ct}
\alias{Xsq.seg.coeff}
\title{Chi-square Test for Segregation Coefficients}
\usage{
Xsq.seg.coeff.ct(ct, covSC)

Xsq.seg.coeff(dat, lab, ...)
}
\arguments{
\item{ct}{A nearest neighbor contingency table, used in \code{Xsq.seg.coeff.ct} only}

\item{covSC}{The covariance matrix for the segregation coefficients in the NNCT, used in \code{Xsq.seg.coeff.ct} only.
Usually output of the function \code{\link{cov.seg.coeff}}}

\item{dat}{The data set in one or higher dimensions, each row corresponds to a data point,
used in \code{Xsq.seg.coeff} only}

\item{lab}{The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.seg.coeff} only}

\item{\dots}{are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
used in \code{Xsq.seg.coeff} only}
}
\value{
A \code{list} with the elements
\item{statistic}{The chi-squared test statistic for the combined segregation coefficients}
\item{p.value}{The \eqn{p}-value for the hypothesis test}
\item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k+1)/2-1} for this function.}
\item{estimate}{The \code{vector} of estimates of the parameters, i.e., observed values of segregation coefficients 
in the NNCT.}
\item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
\item{null.value}{The null value of the parameters, i.e., expected values of segregation coefficients
in the NNCT under RL or CSR (which is 0).}
\item{method}{Description of the hypothesis test}
\item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Xsq.seg.coeff.ct} only}
\item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.seg.coeff} only}
}
\description{
Two functions: \code{Xsq.seg.coeff.ct} and \code{Xsq.seg.coeff}.

Each one performs hypothesis tests of (simultaneous) equality of the segregation coefficients in an NNCT
to the ones under RL or CSR.
That is, each performs the combined Chi-square test for segregation coefficients which is appropriate 
(i.e. have the appropriate asymptotic sampling distribution)
for completely mapped data.
(See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).

Each test is based on the Chi-square approximation of the corresponding quadratic form for the segregation
coefficients in an NNCT.
The segregation coefficients in the multi-class case are the extension of Pielou's segregation coefficient
for the two-class case.
(See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).

Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k+1)/2-1}, description of the 
alternative with the corresponding null values (i.e. expected values) of the segregation coefficients in the NNCT
(which are 0 for this function) and also the sample estimates (i.e. observed values) of the segregation
coefficients. The functions also provide names of the test statistics, the method and the data set used.

The null hypothesis for all cells \eqn{(i,j)} is that the corresponding segregation coefficients are all 
equal to the expected value (which is 0) under RL or CSR.
}
\examples{
n<-20
Y<-matrix(runif(3*n),ncol=3)
ipd<-ipd.mat(Y)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
ct<-nnct(ipd,cls)

W<-Wmat(ipd)
Qv<-Qvec(W)$q
Rv<-Rval(W)
varN<-var.nnct(ct,Qv,Rv)
covN<-cov.nnct(ct,varN,Qv,Rv)

covSC<-cov.seg.coeff(ct,covN)

Xsq.seg.coeff(Y,cls)
Xsq.seg.coeff.ct(ct,covSC)

Xsq.seg.coeff(Y,cls,method="max")

#cls as a factor
na<-floor(n/2); nb<-n-na
fcls<-rep(c("a","b"),c(na,nb))
ct<-nnct(ipd,fcls)

Xsq.seg.coeff.ct(ct,covSC)

#############
n<-40
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
ipd<-ipd.mat(Y)
ct<-nnct(ipd,cls)

W<-Wmat(ipd)
Qv<-Qvec(W)$q
Rv<-Rval(W)
varN<-var.nnct(ct,Qv,Rv)
covN<-cov.nnct(ct,varN,Qv,Rv)

covSC<-cov.seg.coeff(ct,covN)

Xsq.seg.coeff(Y,cls)
Xsq.seg.coeff.ct(ct,covSC)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{seg.coeff}}, \code{\link{Zseg.coeff.ct}} and \code{\link{Zseg.coeff}}
}
\author{
Elvan Ceyhan
}
