% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{pairwise.lab}
\alias{pairwise.lab}
\title{Keeping the pair of the specified labels in the data}
\usage{
pairwise.lab(dat, lab, i, j)
}
\arguments{
\item{dat}{The data set in one or higher dimensions, each row corresponds to a data point.}

\item{lab}{The \code{vector} of class labels (numerical or categorical)}

\item{i, j}{Label of the classes that are to be retained in the post-hoc comparison.}
}
\value{
A \code{list} with two elements
\item{data.pair}{The type of the pattern from which points are to be generated}
\item{lab.pair}{The \code{"main"} title for the plot of the point pattern}
}
\description{
Keeps only the specified labels \eqn{i} and \eqn{j} and returns the data from classes with these labes and also
the corresponding label vector having class labels \eqn{i} and \eqn{j} only.

See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
}
\examples{
n<-20  #or try sample(1:20,1)
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
pairwise.lab(Y,cls,1,2)

#############
n<-40
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
pairwise.lab(Y,cls,2,3)

#cls as a factor
fcls<-rep(letters[1:4],rep(10,4))
pairwise.lab(Y,fcls,"b","c")

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{lab.onevsrest}} and \code{\link{classirest}}
}
\author{
Elvan Ceyhan
}
