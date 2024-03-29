% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{nnct.boot.dis}
\alias{nnct.boot.dis}
\title{Bootstrap Nearest Neighbor Contingency Table (NNCT)}
\usage{
nnct.boot.dis(x, lab, self = TRUE, ties = TRUE, is.ipd = TRUE, ...)
}
\arguments{
\item{x}{The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
correspond to the rows (if \code{is.ipd = FALSEALSE}).}

\item{lab}{The \code{vector} of class labels (numerical or categorical)}

\item{self}{A logical argument (default=\code{TRUE}). If \code{TRUE}, for each base point, all entries in the 
row are sampled (with replacement) and if \code{FALSE} the point is excluded from the resampling (i.e. other points
are sampled with replacement).}

\item{ties}{A logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.}

\item{is.ipd}{A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
matrix, otherwise, \code{x} is taken as the data set with rows representing the data points.}

\item{...}{are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.}
}
\value{
Returns the \eqn{k \times k} NNCT where \eqn{k} is the number of classes in the data set with sampling replacement
of the rows of the IPD matrix.
}
\description{
Returns the \eqn{k \times k} NNCT with sampling replacement of the points for each base point. That is, for each base 
point, the rows in the IPD matrix are sampled with replacement and the NN counts are updated accordingly.
Row and columns of the NNCT are labeled with the corresponding class labels.

The argument self is a logical argument (default=\code{TRUE}) for including the base point in the resampling or not.
If \code{TRUE}, for each base point all entries in the row are sampled (with replacement) so the point itself can
also be sampled multiple times and if \code{FALSE} the point is excluded from the resampling (i.e. other points
are sampled with replacement).

The argument \code{ties} is a logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.

The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
with rows representing the data points.
}
\examples{
n<-20  #or try sample(1:20,1)
Y<-matrix(runif(3*n),ncol=3)
ipd<-ipd.mat(Y)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
nnct.boot.dis(ipd,cls)
nnct.boot.dis(Y,cls,is.ipd = FALSE) #may give different result from above due to random sub-sampling
nnct.boot.dis(ipd,cls,self = FALSE)
nnct.boot.dis(ipd,cls,ties = FALSE) #differences are due to ties and resampling of distances

#cls as a factor
na<-floor(n/2); nb<-n-na
fcls<-rep(c("a","b"),c(na,nb))
nnct.boot.dis(ipd,fcls)

#cls as an unsorted factor
fcls<-sample(c("a","b"),n,replace = TRUE)
nnct.boot.dis(ipd,fcls)

fcls<-sort(fcls)
nnct.boot.dis(ipd,fcls)

#############
n<-40
Y<-matrix(runif(3*n),ncol=3)
ipd<-ipd.mat(Y)
cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
nnct.boot.dis(ipd,cls)

#cls as a factor
fcls<-rep(letters[1:4],rep(10,4))
nnct.boot.dis(ipd,fcls)

#1D data points
n<-20  #or try sample(1:20,1)
X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#(i.e., a column vector), hence X<-runif(n) would not work
ipd<-ipd.mat(X)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
nnct.boot.dis(ipd,cls)

#cls as a factor
na<-floor(n/2); nb<-n-na
fcls<-rep(c("a","b"),c(na,nb))
nnct.boot.dis(ipd,fcls)

#with possible ties in the data
Y<-matrix(round(runif(3*n)*10),ncol=3)
ipd<-ipd.mat(Y)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
nnct.boot.dis(ipd,cls)
nnct.boot.dis(ipd,cls,self = FALSE)
nnct.boot.dis(ipd,cls,ties = FALSE) #differences are due to ties and resampling of distances

}
\seealso{
\code{\link{nnct}} and \code{\link{nnct.sub}}
}
\author{
Elvan Ceyhan
}
