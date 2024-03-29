% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNCTFunctions.R
\name{exact.pval2s}
\alias{exact.pval2s}
\title{\eqn{p}-value correction to the two-sided version of exact NNCT test}
\usage{
exact.pval2s(ptable, pval, type = "inc", double = FALSE)
}
\arguments{
\item{ptable}{Probability of the observed \eqn{2 \times 2} NNCT under the null hypothesis using the hypergeometric
distribution for Fisher's exact test.}

\item{pval}{Table inclusive \eqn{p}-value for Fisher's exact test on the NNCT.}

\item{type}{The type of the \eqn{p}-value correction for the two-sided exact test on the NNCT, default=\code{"inc"}.
Takes on values \code{"inc"}, \code{"exc"}, \code{"mid"}, \code{"tocher"} (or equivalently \code{1-4}, respectively) for table inclusive,
table-exclusive, mid-\eqn{p}-value, and Tocher corrected \eqn{p}-value, respectively.}

\item{double}{A logical argument (default is \code{FALSE}) to determine whether type I or II correction should be 
applied to the two-sided \eqn{p}-value. If \code{TRUE} type I correction (for doubling the minimum of the one-sided \eqn{p}-value) 
is applied, otherwise, type II correction (using the probabilities for the more extreme tables) is applied.}
}
\value{
A modified \eqn{p}-value based on the correction specified in \code{type}.
}
\description{
In using Fisher's exact test on the \eqn{2 \times 2} nearest neighbor contingency tables (NNCTs) a correction may be needed
for the \eqn{p}-value. For the one-sided alternatives, the probabilities of 
more extreme tables are summed up, including or excluding the 
probability of the table itself (or some middle way). 

There is additional complexity in \eqn{p}-values for the two-sided alternatives.
A recommended method is adding up probabilities of the same
size and smaller than the probability associated with the current table.
Alternatively, one can double the one-sided \eqn{p}-value (see (\insertCite{agresti:1992;textual}{nnspat}).

Let the probability of the contingency table itself be \eqn{p_t=f(n_{11}|n_1,n_2,c_1;\theta_0)}
where \eqn{\theta_0=(n_1-1)(n_2-1)/(n_1 n_2)} which is the odds ratio
under RL or CSR independence and
\eqn{f} is the probability mass function of the hypergeometric distribution.

**Type (I):** For double the one-sided \eqn{p}-value, we propose the following four variants:
  \itemize{
\item [(i)] twice the minimum of \eqn{p_{inc}} for the one-sided tests, which is
table-inclusive version for this type of two-sided test, and denoted as \eqn{p^I_{inc}},
\item [(ii)] twice the minimum of \eqn{p_{inc}} minus twice the table
probability \eqn{p_t}, which is table-exclusive version of this type of
two-sided test, and denoted as \eqn{p^I_{exc}},
\item [(iii)] table-exclusive version of this type of
two-sided test plus \eqn{p_t}, which is mid-\eqn{p}-value for
this test, and denoted as \eqn{p^I_{midd}},
\item [(iv)]Tocher corrected version (see \code{\link{tocher.cor}} for details).
}

**Type (II):** For summing the \eqn{p}-values of more extreme ---than that of the table--- cases
in both directions, the following variants are obtained.
The \eqn{p}-value is \eqn{p=\sum_S f(t|n_1,n_2,c_1;\theta=1)} with
\itemize{
\item [(i)] \eqn{S=\{t:\,f(t|n_1,n_2,c_1;\theta=1) \leq p_t\}}, which is
called \emph{table-inclusive version}, \eqn{p^{II}_{inc}},
\item [(ii)] the probability of the observed table is included twice, once for each side;
that is \eqn{p=p^{II}_{inc}+p_t}, which is called \emph{twice-table-inclusive version}, \eqn{p^{II}_{tinc}},
\item [(iii)] table-inclusive minus \eqn{p_t}, which is referred as \emph{table-exclusive version}, \eqn{p^{II}_{exc}},
\item [(iv)] table-exclusive plus one-half
the \eqn{p_t}, which is called \emph{mid-\eqn{p} version}, \eqn{p^{II}_{mid}} and,
\item [(v)]\emph{Tocher corrected version}, \eqn{p^{II}_{Toc}}, is obtained as before.
}

See (\insertCite{ceyhan:exact-NNCT;textual}{nnspat}) for more details.
}
\examples{
ct<-matrix(sample(20:40,4),ncol=2)
ptab<-prob.nnct(ct)
pv<-.23
exact.pval2s(ptab,pv)
exact.pval2s(ptab,pv,type="exc")
exact.pval2s(ptab,pv,type="mid")

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{exact.pval1s}} and \code{\link{tocher.cor}}
}
\author{
Elvan Ceyhan
}
