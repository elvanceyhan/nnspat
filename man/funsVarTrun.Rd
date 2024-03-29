% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{funsVarTrun}
\alias{funsVarTrun}
\alias{varTrun}
\alias{varTrun.sim}
\title{Variance of Cuzick and Edwards \eqn{T_{run}} Test statistic}
\usage{
varTrun(dat, n1, ...)

varTrun.sim(dat, cc.lab, Nsim = 1000, case.lab = NULL)
}
\arguments{
\item{dat}{The data set in one or higher dimensions, each row corresponds to a data point, 
used in both functions.}

\item{n1}{Number of cases, used in \code{varTrun} only.}

\item{\dots}{are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
Used in \code{varTrun} only.}

\item{cc.lab}{Case-control labels, 1 for case, 0 for control, used in \code{varTrun.sim} only.}

\item{Nsim}{The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
variance of \eqn{T_{run}}, used in \code{varTrun.sim} only.}

\item{case.lab}{The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
such that cases are 1 and controls are 0), default is \code{NULL}, used in \code{varTrun.sim} only.}
}
\value{
The function \code{varTrun} returns the variance of Cuzick and Edwards \eqn{T_{run}} test statistic
under RL or CSR independence.
And the function \code{varTrun.sim} estimates the same variance based on simulations under the RL hypothesis.
}
\description{
Two functions: \code{varTrun} and \code{varTrun.sim}.

The function \code{varTrun} computes the (finite sample) variance of Cuzick and Edwards \eqn{T_{run}} test statistic 
which is based on the number of consecutive cases from the cases in the data under RL or CSR independence.
And the function \code{varTrun.sim} estimates this variance based on simulations under the RL hypothesis.

The only common argument for both functions is \code{dat}, the data set used in the functions.

\eqn{n_1} is an argument for \code{varTrun} and is the number of cases (denoted as \code{n1} as an argument).
The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
to match the case-control class labeling,
which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.

The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
and 1's accordingly. The argument \code{Nsim} represents the number of resamplings (without replacement) in the
RL scheme, with default being \code{1000}. \code{cc.lab}, \code{case.lab} and \code{Nsim} are arguments for \code{varTrun.sim} only.

The function \code{varTrun} might take a very long time when data size is large (even larger than 50),
hence the need for the \code{varTrun.sim} function. 

See (\insertCite{cuzick:1990;textual}{nnspat}).
}
\examples{
n<-20  #or try sample(1:20,1) #try also 40, 50, 60
set.seed(123)
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(0:1,n,replace = TRUE)
n1<-sum(cls==1)
n0<-sum(cls==0)
c(n1,n0)

varTrun(Y,n1)
varTrun(Y,n1,method="max")

n<-15  #or try sample(1:20,1) #try also 40, 50, 60
set.seed(123)
Y<-matrix(runif(3*n),ncol=3)
cls<-sample(0:1,n,replace = TRUE)
n1<-sum(cls==1)
varTrun(Y,n1) #the actual value (might take a long time if n is large)

Nmc<-1000
varTrun.sim(Y,cls,Nsim=Nmc)

#cls as a factor
na<-floor(n/2); nb<-n-na
fcls<-rep(c("a","b"),c(na,nb))
varTrun.sim(Y,fcls,Nsim=Nmc,case.lab="a")

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{ceTrun}} and \code{\link{EV.Trun}}
}
\author{
Elvan Ceyhan

Elvan Ceyhan
}
