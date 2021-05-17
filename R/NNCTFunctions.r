#NNCTFunctions.r;
#Contains the ancillary functions used in NNCT calculations, such as NNCT for two classes of points

#################################################################
#####AUXILIARY FUNCTIONS#################
#################################################################
#in all these functions
#data sets are either matrices or data frames

#' @import stats
#' @import MASS
#' @import graphics
#'
#' @title Interpoint Distance Matrix
#'
#' @description 
#' This function computes and returns the distance matrix computed by using the specified distance measure to
#' compute the distances between the rows of the set of points \code{x} and \code{y} using the 
#' \code{\link[stats]{dist}} function in the \code{stats} package of the standard R distribution.
#' If \code{y} is provided (default=\code{NULL}) it yields a matrix of distances between the rows of \code{x} and rows of \code{y},
#' otherwise, it provides a square matrix with i,j-th entry being the distance between row \eqn{i} and row \eqn{j} of \code{x}.
#' This function is different from the \code{\link[stats]{dist}} function in the \code{stats} package.
#' \code{dist} returns the distance matrix in a lower triangular form, and \code{ipd.mat} returns in a full matrix.
#' \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @param x A set of points in matrix or data frame form where points correspond to the rows.
#' @param y A set of points in matrix or data frame form where points correspond to the rows (default=\code{NULL}).
#' @param \dots Additional parameters to be passed on the \code{dist} function.
#'
#' @return A distance matrix whose i,j-th entry is the distance between row \eqn{i} of \code{x} and row \eqn{j} of \code{y} if \code{y} is provided,
#' otherwise i,j-th entry is the distance between rows \eqn{i} and \eqn{j} of \code{x}.
#'
#' @seealso \code{\link[stats]{dist}}, \code{\link{ipd.mat.euc}}, \code{\link{dist.std.data}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-3
#' X<-matrix(runif(3*n),ncol=3)
#' mtd<-"euclidean" #try also "maximum", "manhattan", "canberra", "binary"
#' ipd.mat(X,method=mtd)
#' ipd.mat(X,method="minkowski",p=6)
#'
#' n<-5
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd.mat(X,Y,method=mtd)
#' ipd.mat(X[1,],Y,method=mtd)
#' ipd.mat(c(.1,.2,.3),Y,method=mtd)
#' ipd.mat(X[1,],Y[3,],method=mtd)
#'
#' #1D data points
#' X<-as.matrix(runif(3)) # need to be entered as a matrix with one column 
#' #(i.e., a column vector), hence X<-runif(3) would not work
#' ipd.mat(X)
#'
#' Y<-as.matrix(runif(5))
#' ipd.mat(X,Y)
#' ipd.mat(X[1,],Y)
#' ipd.mat(X[1,],Y[3,])
#'
#' @export
ipd.mat <- function(x,y=NULL, ...)
{
  if (!is.numeric(x) | (!is.null(y) & !is.numeric(y)))
  {stop('first argument and (if provided) second argument must be of type numeric')}
  if (is.null(y))
  {dis<-dist(x,...)
  ipd<-dist2full(dis)
  } else
  {
    ifelse(is.vector(x),x<-matrix(x,ncol=length(x)),x<-as.matrix(x))
    ifelse(is.vector(y),y<-matrix(y,ncol=length(y)),y<-as.matrix(y))
    xy<-rbind(x,y)
    nx<-nrow(x)
    ny<-nrow(y)
    dis<-dist(xy,...)
    ipd0<-dist2full(dis)
    ipd<-ipd0[1:nx,(nx+1):(nx+ny)]
  }
  ipd
} #end for the function
#' 

#################################################################

#' @title Converts a lower triangular distance matrix to a full distance matrix
#'
#' @description 
#' Converts a lower triangular distance matrix to a full distance matrix with zeroes in the diagonal.
#' The input is usually the result of the \code{\link[stats]{dist}} function in the \code{stats} package.
#' This function is adapted from Everitt's book (\insertCite{everitt:2004;textual}{nnspat})
#' 
#' @param dis A lower triangular matrix, resulting from the \code{\link[stats]{dist}} 
#' function in the \code{stats} package
#'
#' @return A square (symmetric) distance matrix with zeroes in the diagonal.
#'
#' @seealso \code{\link[stats]{dist}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-3
#' X<-matrix(runif(3*n),ncol=3)
#' dst<-dist(X)
#' dist2full(dst)
#'
#' @export
dist2full <- function(dis) {
  n<-attr(dis,"Size")
  full<-matrix(0,n,n)
  full[lower.tri(full)]<-dis
  full+t(full)
} #end for the function

#################################################################

#' @title Euclidean Interpoint Distance Matrix
#'
#' @description
#' Returns the Euclidean interpoint distance (IPD) matrix of a given the set of points \code{x} and \code{y} using 
#' two for loops with the \code{\link{euc.dist}} function of the current package.
#' If \code{y} is provided (default=\code{NULL}) it yields a matrix of Euclidean distances between the rows of \code{x} and rows of \code{y},
#' otherwise it provides a square matrix with i,j-th entry being the Euclidean distance between row \eqn{i} and row \eqn{j} 
#' of \code{x}. This function is different from the \code{\link{ipd.mat}} function in this package.
#' \code{\link{ipd.mat}} returns the full distance matrix for a variety of distance metrics (including the
#' Euclidean metric), while \code{\link{ipd.mat.euc}} uses the Euclidean distance metric only.
#' \code{ipd.mat.euc(X)} and \code{ipd.mat(X)} yield the same output for a set of points \code{X},
#' as the default metric in \code{\link{ipd.mat}} is also \code{"euclidean"}.
#' 
#' @param x A set of points in matrix or data frame form where points correspond to the rows.
#' @param y A set of points in matrix or data frame form where points correspond to the rows (default=\code{NULL}).
#'
#' @return A distance matrix whose i,j-th entry is the Euclidean distance between row \eqn{i} of \code{x} and
#' row \eqn{j} of \code{y} if \code{y} is provided, otherwise i,j-th entry is 
#' the Euclidean distance between rows \eqn{i} and \eqn{j} of \code{x}.
#'
#' @seealso \code{\link[stats]{dist}}, \code{\link{ipd.mat.euc}}, \code{\link{dist.std.data}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-3
#' X<-matrix(runif(3*n),ncol=3)
#' ipd.mat.euc(X)
#'
#' n<-5
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd.mat.euc(X,Y)
#' ipd.mat.euc(X[1,],Y)
#' ipd.mat.euc(c(.1,.2,.3),Y)
#' ipd.mat.euc(X[1,],Y[3,])
#'
#' #1D data points
#' X<-as.matrix(runif(3)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(3) would not work
#' ipd.mat.euc(X)
#'
#' Y<-as.matrix(runif(5))
#' ipd.mat.euc(X,Y)
#' ipd.mat.euc(X[1,],Y)
#' ipd.mat.euc(X[1,],Y[3,])
#'
#' @export
ipd.mat.euc <- function(x,y=NULL)
{
  ifelse(is.vector(x),x<-matrix(x,ncol=length(x)),x<-as.matrix(x))
  nx<-nrow(x)
  
  if (nx==1)
  {ipd<-euc.dist(x,x)
  } else
  {
    if (is.null(y))
    {
      ipd<-matrix(0,nx,nx)
      for (i in 1:(nx-1))
      {
        for (j in (i+1):nx)
        {
          ipd[i,j]<-euc.dist(x[i,],x[j,]);
          ipd[j,i]<-ipd[i,j]
        }
      }
      
    } else
    {
      ifelse(is.vector(y),y<-matrix(y,ncol=length(y)),y<-as.matrix(y))
      ny<-nrow(y)
      ipd<-matrix(0,nx,ny)
      for (i in 1:nx)
      {
        for (j in 1:ny)
        {
          ipd[i,j]<-euc.dist(x[i,],y[j,]);
        }
      }
      
    }
  }
  ipd
} #end for the function

#################################################################

#' @title The Euclidean distance between two vectors, matrices, or data frames
#'
#' @description Returns the Euclidean distance between \code{x} and \code{y} which can be vectors
#' or matrices or data frames of any dimension (\code{x} and \code{y} should be of same dimension).
#'
#' This function is equivalent to \code{\link[pcds]{Dist}} function in the \code{pcds} package but is 
#' different from the \code{\link[stats]{dist}} function in the \code{stats} package of the standard R 
#' distribution.
#' \code{dist} requires its argument to be a data matrix and \code{\link[stats]{dist}} computes and returns the distance matrix computed
#' by using the specified distance measure to compute the distances between the rows of a data matrix
#' (\insertCite{S-Book:1998;textual}{nnspat}),
#' while \code{euc.dist} needs two arguments to find the distances between. 
#' For two data matrices \code{A} and \code{B},
#' \code{dist(rbind(as.vector(A),as.vector(B)))} and \code{euc.dist(A,B)} yield the same result.
#'
#' @param x,y Vectors, matrices or data frames (both should be of the same type).
#'
#' @return Euclidean distance between \code{x} and \code{y}
#'
#' @seealso \code{\link[stats]{dist}} from the base package \code{stats} and
#' \code{\link[pcds]{Dist}} from the package \code{pcds}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' euc.dist(B,C);
#' euc.dist(B,B);
#' #'
#' x<-runif(10)
#' y<-runif(10)
#' euc.dist(x,y)
#' #'
#' xm<-matrix(x,ncol=2)
#' ym<-matrix(y,ncol=2)
#' euc.dist(xm,ym)
#' #'
#' euc.dist(xm,xm)
#' #'
#' dat.fr<-data.frame(b=B,c=C)
#' euc.dist(dat.fr,dat.fr)
#' euc.dist(dat.fr,cbind(B,C))
#'
#' @export
euc.dist <- function(x,y)
{
  x<-as.matrix(x)
  y<-as.matrix(y)
  dis<-sqrt(sum((x-y)^2))
  
  dis
} #end of the function
#'

#################################################################

#' @title The incidence matrix \code{W} for the NN digraph
#'
#' @description 
#' Returns the \eqn{W=(w_ij)} matrix which is used to compute \eqn{Q}, \eqn{R} and \eqn{T} values in the NN structure.
#' \eqn{w_{ij}=I(} point eqn{j} is a NN of point \eqn{i))} i.e. \eqn{w_{ij}=1} if point \eqn{j} is a NN of point \eqn{i} and 0 otherwise.
#' 
#' The argument \code{ties} is a logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} the function
#' takes ties into account by making \eqn{w_{ij}=1/m} if point \eqn{j} is a NN of point \eqn{i}
#' and there are \eqn{m} tied NNs and 0 otherwise. If \code{FALSE}, \eqn{w_{ij}=1} if point \eqn{j} is a NN of point \eqn{i} and 0 otherwise.
#' The matrix \eqn{W} is equivalent to \eqn{A=(a_{ij})} matrix with \eqn{k=1}, i.e., \code{Wmat(X)=aij.mat(X,k=1)}.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param ties A logical parameter (default=\code{FALSE}) to take ties into account in computing the \eqn{W} matrix,
#' so if it is \code{TRUE}, \eqn{w_{ij}=1/m} if point \eqn{j} is a NN of point \eqn{i} and there are \eqn{m} tied NNs and 0 otherwise
#' and if \code{FALSE}, \eqn{w_{ij}=1} if point \eqn{j} is a NN of point \eqn{i} and 0 otherwise.
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return The incidence matrix \eqn{W=(w_ij)} where \eqn{w_{ij}=I(} point eqn{j} is a NN of point \eqn{i))},
#' i.e. \eqn{w_{ij}=1} if point \eqn{j} is a NN of point \eqn{i} and 0 otherwise.
#'
#' @seealso \code{\link{aij.mat}}, \code{\link{aij.nonzero}}, and \code{\link{aij.theta}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-3
#' X<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(X)
#' Wmat(ipd)
#' Wmat(X,is.ipd = FALSE)
#'
#' n<-5
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' Wmat(ipd)
#' Wmat(Y,is.ipd = FALSE)
#' Wmat(Y,is.ipd = FALSE,method="max")
#' 
#' Wmat(Y,is.ipd = FALSE)
#' aij.mat(Y,k=1)
#'
#' #1D data points
#' X<-as.matrix(runif(5)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' Wmat(ipd)
#' Wmat(X,is.ipd = FALSE)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(15)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' Wmat(ipd,ties=TRUE)
#' Wmat(Y,ties=TRUE,is.ipd = FALSE)
#'
#' @export
Wmat <- function(x,ties=FALSE,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  W<-matrix(0,n,n)
  for (i in 1:n)
  {
    if (ties==FALSE)
    {
      W[i,ipd[i,]==sort(ipd[i,])[2]]<-1
    } else
    {
      mindis<-sort(ipd[i,])[2]
      ind.md<-which(ipd[i,]==mindis)
      lties<-sum(ind.md!=i) 
      W[i,ind.md]<-1/lties 
    }
    W[i,i]<-0
  }
  W
} #end for the function

################################################################# 

# funsQandR
#'
#' @title Functions for the Number of Shared NNs, Shared NN vector and the number of reflexive NNs
#'
#' @description
#' Four functions: \code{Qval}, \code{Qvec}, \code{sharedNN} and \code{Rval}.
#'
#' \code{Qval} returns the \eqn{Q} value, the number of points with shared nearest neighbors (NNs), which occurs when two or 
#' more points share a NN, for data in any dimension. 
#' 
#' \code{Qvec} returns the Q-value and also yields the Qv vector \eqn{Qv=(Q_0,Q_1,\ldots)} as well for data in any 
#' dimension, where \eqn{Q_j} is the number of points shared as a NN by \eqn{j} other points. 
#' 
#' \code{sharedNN} returns the \code{vector} of number of points with shared NNs, \eqn{Q=(Q_0,Q_1,\ldots)} where \eqn{Q_i} is
#' the number of points that are NN to \eqn{i} points, and if a point is a NN of \eqn{i} points, then there are \eqn{i(i-1)} 
#' points that share a NN. So \eqn{Q=\sum_{i>1} i(i-1)Q_i}.
#' 
#' \code{Rval} returns the number of reflexive NNs, R (i.e., twice the number of reflexive NN pairs).
#' 
#' These quantities are used, e.g., in computing the variances and covariances of the entries of the
#' nearest neighbor contingency tables used for Dixon's tests and other NNCT tests. 
#' The input must be the incidence matrix, \eqn{W}, of the NN digraph. 
#'
#' @param W The incidence matrix, \eqn{W}, for the NN digraph
#'
#' @return \code{Qval} returns the \eqn{Q} value
#' \code{Qvec} returns a \code{list} with two elements
#'  \item{q}{the \eqn{Q} value, the number of shared NNs}
#'  \item{qvec}{the \code{vector} of \eqn{Q_j} values} 
#' \code{sharedNN} returns a \code{matrix} with 2 rows, where first row is the \eqn{j} values and second row is
#' the corresponding vector of \eqn{Q_j} values
#' \code{Rval}{the \eqn{R} value, the number of reflexive NNs}
#' 
#' See the description above for the details of these quantities.
#' 
#' @seealso \code{\link{Tval}}, \code{\link{QRval}}, \code{\link{sharedNNmc}} and \code{\link{Ninv}}
#' 
#' @name funsQandR
NULL
#'
#' @rdname funsQandR
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for Qval
#' #3D data points
#' n<-10
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qval(W)
#'
#' #1D data points
#' X<-as.matrix(runif(10)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(10) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' Qval(W)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(15)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' Qval(W)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(15)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' Qval(W)
#'
#' @export
Qval <- function(W)
{
  W<-ceiling(W) #this is to correct for the ties in W matrix
  a<-colSums(W)
  Q<-sum(a*a-a)
  Q
} #end for the function
#'
#' @rdname funsQandR
#'
#' @examples
#' #Examples for Qvec
#' #3D data points
#' n<-10
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qvec(W)
#'
#' #2D data points
#' n<-15
#' Y<-matrix(runif(2*n),ncol=2)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qvec(W)
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(15) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' Qvec(W)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(15)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' Qvec(W)
#'
#' @export
Qvec <- function(W)
{
  W<-ceiling(W) #this is to correct for the ties in W matrix
  colsumW<-colSums(W)
  tb.csW<-table(colsumW)
  Nv<-as.numeric(tb.csW)
  
  Nlen<-length(Nv)
  Nfac<-(1:Nlen-1)*(1:Nlen-2) #this is the i*(i-1) vector
  qval<-sum(Nfac*Nv)
  
  list(q=qval,qvec=Nv) 
} #end for the function
#'
#' @rdname funsQandR
#'
#' @examples
#' #Examples for sharedNN
#' #3D data points
#' n<-10
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#'
#' #2D data points
#' n<-15
#' Y<-matrix(runif(2*n),ncol=2)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' sharedNN(W)
#'
#' @export
sharedNN <- function(W)
{
  W<-ceiling(W) #this is to correct for the ties in W matrix
  CS<-colSums(W)
  tb.cs<-table(CS)
  sh.nn1<-as.numeric(labels(tb.cs)[[1]])
  sh.nn2<-as.numeric(tb.cs)
  res<-rbind(sh.nn1,sh.nn2)
  row.names(res)<-c()
  res
} #end for the function
#'
#' @rdname funsQandR
#'
#' @examples
#' #Examples for Rval
#' #3D data points
#' n<-10
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Rval(W)
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' Rval(W)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' Rval(W)
#'
#' @export
Rval <- function(W)
{
  W<-ceiling(W) #this is to correct for the ties in W matrix
  sum(W*t(W))
} #end for the function
#'

#################################################################

#' @title Vector of Shared NNs and Number of Reflexive NNs
#'
#' @description Returns the \code{Qvec} and \code{R} where \eqn{Qvec=(Q_0,Q_1,\ldots)} with
#' \eqn{Q_j} is the number of points shared as a NN
#' by \eqn{j} other points i.e. number of points that are NN of \eqn{i} points, for \eqn{i=0,1,2,\ldots}
#' and \code{R} is the number of reflexive pairs where A and B are reflexive iff they are NN to each other.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return Returns a \code{list} with two elements
#'  \item{Qvec}{vector of \eqn{Q_j} values}
#'  \item{R}{number of reflexive points}
#'
#' @seealso \code{\link{Qval}}, \code{\link{Qvec}}, \code{\link{sharedNN}}, \code{\link{Rval}} 
#' and \code{\link{QRval}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#' Ninv(ipd)
#' Ninv(Y,is.ipd = FALSE)
#' Ninv(Y,is.ipd = FALSE,method="max")
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column 
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#' Ninv(ipd)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ny<-nrow(Y)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' sharedNN(W)
#' Qvec(W)
#' Ninv(ipd)
#'
#' @export
Ninv <- function(x,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  inv.nn <- rep(0,n); 
  Rf <- 0;
  NN.ls <- list() 
  
  for (i in 1:n)
    NN.ls[[i]] <- NN(ipd,i); #list of indices of NNs of subject i
  
  for (k in 1:n)
  {
    NN.indk<-NN.ls[[k]]
    L<-length(NN.indk)
    for (j in 1:L)
    {
      inv.nn[NN.indk[j]] <- inv.nn[NN.indk[j]]+1
      if ( sum(k==NN.ls[[NN.indk[j] ]])>0 )
        Rf<- Rf+1
    }
  }
  tb.inv.nn<-table(inv.nn)
  Nv<-as.numeric(tb.inv.nn)
  
  list(Qvec=Nv,
       R=Rf
  )
} #end for the function
#'

#################################################################

#' @title Number of Shared and Reflexive NNs
#'
#' @description Returns the \eqn{Q} and \eqn{R} values where \eqn{Q} is the number of points shared as a NN
#' by other points i.e. number of points that are NN of other points (which occurs when two or 
#' more points share a NN, for data in any dimension) and \eqn{R} is the number of reflexive pairs
#' where A and B are reflexive iff they are NN to each other.
#' 
#' These quantities are used, e.g., in computing the variances and covariances of the entries of the
#' nearest neighbor contingency tables used for Dixon's tests and other NNCT tests. 
#'
#' @param njr A \code{list} that is the output of \code{\link{Ninv}} (with first entry in the \code{list} is \code{vector} of number of shared NNs
#' and second is the \eqn{R} value, number of reflexive points)
#' 
#' @return A \code{list} with two elements
#' \item{Q}{the \eqn{Q} value, the number of shared NNs}
#' \item{R}{the \eqn{R} value, the number of reflexive NNs}
#'
#' @seealso \code{\link{Qval}}, \code{\link{Qvec}}, \code{\link{sharedNN}}, \code{\link{Rval}} 
#' and \code{\link{Ninv}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' ninv<-Ninv(ipd)
#' QRval(ninv)
#' W<-Wmat(ipd)
#' Qvec(W)$q
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column 
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' ninv<-Ninv(ipd)
#' QRval(ninv)
#' W<-Wmat(ipd)
#' Qvec(W)$q
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ny<-nrow(Y)
#' ipd<-ipd.mat(Y)
#' ninv<-Ninv(ipd)
#' QRval(ninv)
#' W<-Wmat(ipd)
#' Qvec(W)$q
#'
#' @export
QRval <- function(njr)
{
  nj<-njr[[1]]
  
  njlen<-length(nj)
  njfac<-(1:njlen-1)*(1:njlen-2) #this is the i*(i-1) vector
  q<-sum(njfac*nj)
  
  r<-njr[[2]]
  
  list(Q=q,R=r)
} #end for the function
#'

################################################################# 

#' @title \eqn{T} value in NN structure
#'
#' @description Returns the \eqn{T} value, which is the number of triplets \eqn{(z_i, z_j, z_k)} with 
#' "\eqn{NN(z_i) = NN(z_j) = z_k} and \eqn{NN(z_k) = z_j}" where \eqn{NN(\cdot)} is the nearest neighbor function.
#' Note that in the NN digraph, \eqn{T+R} is the sum of the indegrees of the points in the reflexive pairs.
#' 
#' This quantity (together with \eqn{Q} and \eqn{R}) is used in computing the variances and covariances of the entries of the
#' reflexivity contingency table. See (\insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}) for further
#' details.
#'
#' @param W The incidence matrix, \eqn{W}, for the NN digraph
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs)
#'
#' @return Returns the \eqn{T} value. See the description above for the details of this quantity.
#'
#' @seealso \code{\link{Qval}}, \code{\link{Qvec}}, \code{\link{sharedNN}} and \code{\link{Rval}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-10
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' R<-Rval(W)
#' Tval(W,R)
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' W<-Wmat(ipd)
#' R<-Rval(W)
#' Tval(W,R)
#'
#' #with ties=TRUE in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd,ties=TRUE)
#' R<-Rval(W)
#' Tval(W,R)
#'
#' @export 
Tval <- function(W,R)
{
  W<-ceiling(W) #this is to correct for the ties in W matrix
  n<-sum(W)
  ref.ind<-which(as.vector(W*t(W)>0))
  col.ref.ind<-ceiling(ref.ind/n)
  sum(W[,col.ref.ind])
  RT<-sum(W[,col.ref.ind]) #sum of the indegrees of points in the reflexive pairs
  Tval<-RT-R
  Tval
} #end for the function
#'

#################################################################

#' @title Finding the index of the NN of a given point
#'
#' @description
#' Returns the index (or indices) of the nearest neighbor(s) of subject \eqn{i} given data set or IPD matrix \code{x}.
#' It will yield a \code{vector} if there are ties, and subject indices correspond to rows (i.e. rows \code{1:n} ) if \code{x} 
#' is the data set and to rows or columns if \code{x} is the IPD matrix.  
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param i index of (i.e., row number for) the subject whose NN is to be found. 
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns the index (indices) i.e. row number(s) of the NN of subject \eqn{i}
#'
#' @seealso \code{\link{kNN}} and \code{\link{NNsub}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' NN(ipd,1)
#' NN(Y,1,is.ipd = FALSE)
#' NN(ipd,5)
#' NN(Y,5,is.ipd = FALSE)
#' NN(Y,5,is.ipd = FALSE,method="max")
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' NN(ipd,1)
#' NN(ipd,5)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ny<-nrow(Y)
#' ipd<-ipd.mat(Y)
#' for (i in 1:ny)
#'   cat(i,":",NN(ipd,i),"|",NN(Y,i,is.ipd = FALSE),"\n")
#'
#' @export 
NN <- function(x,i,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  if (n<=1 || i>n)
  {labind<-NA 
  return(labind)}
  
  D<-max(ipd) 
  ind <- 1:ncol(ipd)
  
  dis.row <- ipd[i,]
  min.dis <- min(dis.row[-i])
  dis.row[i] <- 2*D #to avoid if min distance is zero!
  labind <- ind[dis.row==min.dis]
  
  labind
} #end for the function
#'

################################################################# 

#' @title Finding the indices of the \code{k} NNs of a given point
#'
#' @description
#' Returns the indices of the \code{k} nearest neighbors of subject \eqn{i} given data set or IPD matrix \code{x}.
#' Subject indices correspond to rows (i.e. rows \code{1:n} ) if \code{x} is the data set and to rows or columns
#' if \code{x} is the IPD matrix.  
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#'
#' @inheritParams NN
#' @param k Integer specifying the number of NNs (of subject \eqn{i}).
#'
#' @return Returns the indices (i.e. row numbers) of the \code{k} NNs of subject \eqn{i}
#'
#' @seealso \code{\link{NN}}, \code{\link{NNdist}} and \code{\link{NNdist2cl}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' k<-sample(1:5,1)
#' k
#' NN(ipd,1)
#' kNN(ipd,1,k)
#' kNN(Y,1,k,is.ipd = FALSE)
#' kNN(Y,1,k,is.ipd = FALSE,method="max")
#'
#' NN(ipd,5)
#' kNN(ipd,5,k)
#' kNN(Y,5,k,is.ipd = FALSE)
#'
#' #1D data points
#' X<-as.matrix(runif(15)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' kNN(ipd,3,k)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(30)*10),ncol=3)
#' ny<-nrow(Y)
#' ipd<-ipd.mat(Y)
#' for (i in 1:ny)
#'   cat(i,":",kNN(ipd,i,k),"\n")
#'
#' @export
kNN <- function(x,i,k,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  if (n<=1 || max(i,k+1)>n)
  {knndis<-NA 
  return(knndis)}
  
  D<-max(ipd) 
  ind <- 1:ncol(ipd)
  
  dis.row <- ipd[i,]
  knndis <- order(dis.row)[2:(k+1)]
  knndis
} #end for the function
#'

#################################################################

#' @title Distances between subjects and their NNs
#'
#' @description
#' Returns the distances between subjects and their NNs. The output is an \eqn{n \times 2} matrix where \eqn{n} is the data size
#' and first column is the subject index and second column contains the corresponding distances to NN subjects.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns an \eqn{n \times 2} matrix where \eqn{n} is data size (i.e. number of subjects) and first column is the subject
#' index and second column is the NN distances.
#' 
#' @seealso \code{\link{kthNNdist}}, \code{\link{kNNdist}}, and \code{\link{NNdist2cl}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' NNdist(ipd)
#' NNdist(Y,is.ipd = FALSE)
#' NNdist(Y,is.ipd = FALSE,method="max")
#'
#' #1D data points
#' X<-as.matrix(runif(5)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' NNdist(ipd)
#' NNdist(X,is.ipd = FALSE)
#'
#' @export
NNdist <- function(x,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  if (n<=1)
  {min.dis<-NA 
  return(c(1,min.dis))}
  
  min.dis<-vector()
  for (i in 1:n)
  {
    dis.row <- ipd[i,]
    min.dis <- c(min.dis,min(dis.row[-i]))
  }
  subj.ind<-1:n
  cbind(subj.ind,min.dis)
} #end for the function
#'

################################################################# 

# funs.kNNdist
#'
#' @title Functions for the \eqn{k^{th}} and \code{k} NN distances
#'
#' @description
#' Two functions: \code{kthNNdist} and \code{kNNdist}.
#'
#' \code{kthNNdist} returns the distances between subjects and their \eqn{k^{th}} NNs. The output is an \eqn{n \times 2} matrix where 
#' \eqn{n} is the data size and first column is the subject index and second column contains the corresponding 
#' distances to \eqn{k^{th}} NN subjects. 
#' 
#' \code{kNNdist} returns the distances between subjects and their \code{k} NNs.
#' The output is an \eqn{n \times (k+1)} matrix where 
#' \eqn{n} is the data size and first column is the subject index and the remaining \code{k} columns contain the corresponding 
#' distances to \code{k} NN subjects. 
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param k Integer specifying the number of NNs (of subjects).
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return \code{kthNNdist} returns an \eqn{n \times 2} matrix where \eqn{n} is data size (i.e. number of subjects) and
#' first column is the subject index and second column is the \eqn{k^{th}} NN distances.
#' 
#' \code{kNNdist} returns an \eqn{n \times (k+1)} matrix where \eqn{n} is data size (i.e. number of subjects) and
#' first column is the subject index and the remaining \code{k} columns contain the corresponding 
#' distances to \code{k} NN subjects. 
#' 
#' @seealso \code{\link{NNdist}} and \code{\link{NNdist2cl}}
#' 
#' @name funs.kNNdist
NULL
#'
#' @rdname funs.kNNdist
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for kthNNdist
#' #3D data points, gives NAs when n<=k
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' kthNNdist(ipd,3)
#' kthNNdist(Y,3,is.ipd = FALSE)
#' kthNNdist(ipd,5)
#' kthNNdist(Y,5,is.ipd = FALSE)
#' kthNNdist(Y,3,is.ipd = FALSE,method="max")
#'
#' #1D data points
#' X<-as.matrix(runif(5)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' kthNNdist(ipd,3)
#'
#' @export
kthNNdist <- function(x,k,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  kth.dis<-vector()
  for (i in 1:n)
  {
    dis.row <- ipd[i,]
    kth.dis <- c(kth.dis,dis.row[order(dis.row)[k+1]] )
  }
  subj.ind<-1:n
  res<-list(
    kth.nndist=cbind(subj.ind,kth.dis),
    k=k)
  res
} #end for the function
#'
#' @rdname funs.kNNdist
#'
#' @examples
#' #Examples for kNNdist
#' #3D data points, gives NAs if n<=k for n,n+1,...,kNNs
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' kNNdist(ipd,3)
#' kNNdist(ipd,5)
#' kNNdist(Y,5,is.ipd = FALSE)
#'
#' kNNdist(Y,5,is.ipd = FALSE,method="max")
#'
#' kNNdist(ipd,1)
#' kthNNdist(ipd,1)
#'
#' #1D data points
#' X<-as.matrix(runif(5)) # need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(5) would not work
#' ipd<-ipd.mat(X)
#' kNNdist(ipd,3)
#'
#' @export
kNNdist <- function(x,k,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  kdist<-vector()
  for (i in 1:n)
  {
    dist.row <- ipd[i,]
    kdist <- rbind(kdist,dist.row[order(dist.row)[2:(k+1)]] )
  }
  
  res<-list(
    knndist=cbind(1:n,kdist),
    k=k)
  res
} #end for the function
#'

#################################################################

#' @title Distances between subjects from class \eqn{i} and their NNs from class \eqn{j}
#'
#' @description
#' Returns the distances between subjects from class \eqn{i} and their nearest neighbors (NNs) from class \eqn{j}. 
#' The output is a \code{list} with first entry (\code{nndist}) being an \eqn{n_i \times 3} matrix where \eqn{n_i} is the size of class \eqn{i}
#' and first column is the subject index in class \eqn{i}, second column is the subject index in NN class \eqn{j},  
#' and third column contains the corresponding distances of each class \eqn{i} subject to its NN among class \eqn{j}
#' subjects. Class \eqn{i} is labeled as base class and class \eqn{j} is labeled as NN class.
#' 
#' The argument \code{within.class.ind} is a logical argument (default=\code{FALSE}) to determine the indexing of 
#' the class \eqn{i} subjects. If \code{TRUE}, index numbering of subjects is within the class, 
#' from 1 to class size (i.e., \code{1:n_i}), according to their order in the original data;
#' otherwise, index numbering within class is just the indices in the original data.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param i,j class label of base class and NN classes, respectively.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param within.class.ind A logical parameter (default=\code{FALSE}). If \code{TRUE}, index numbering of subjects 
#' is within the class, from 1 to class size (i.e., \code{1:n_i}), according to their order in the original data;
#' otherwise, index numbering within class is just the indices in the original data.
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns a \code{list} with three elements
#'  \item{nndist}{\eqn{n_i \times 3} matrix where \eqn{n_i} is the size of class \eqn{i} and first column is the subject index in 
#'  class \eqn{i}, second column is the subject index in NN class \eqn{j}, and third column contains the corresponding
#'  distances of each class \eqn{i} subject to its NN among class \eqn{j} subjects.}
#'  \item{base.class}{label of base class} 
#'  \item{nn.class}{label of NN class} 
#' 
#' @seealso \code{\link{kthNNdist}}, \code{\link{kNNdist}}, and \code{\link{NNdist2cl}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#' NNdist2cl(ipd,1,2,clab)
#' NNdist2cl(Y,1,2,clab,is.ipd = FALSE)
#'
#' NNdist2cl(ipd,1,2,clab,within = TRUE)
#'
#' #three class case
#' clab<-sample(1:3,n,replace=TRUE) #class labels
#' table(clab)
#' NNdist2cl(ipd,2,1,clab)
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#' NNdist2cl(ipd,1,2,clab)
#' NNdist2cl(X,1,2,clab,is.ipd = FALSE)
#'
#' @export 
NNdist2cl <- function(x,i,j,lab,within.class.ind=FALSE,is.ipd=TRUE,...)
{ 
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  if (any(lab==i)*any(lab==j)==0)
  {stop('given labels i and j are not among the class labels')}
  
  if (within.class.ind==FALSE)
  {
    ns<-length(lab) #sample size of data for which IPDM is computed
    i.ind<-(1:ns)[lab==i]
    j.ind<-(1:ns)[lab==j]
    ni<-length(i.ind)
  } else
  {
    ni<-sum(lab==i)
    nj<-sum(lab==j)
    i.ind<-1:ni
    j.ind<-1:nj 
  }
  
  ipdij<-matrix(ipd[lab==i,lab==j],nrow=ni)
  min.dist<-nnj.ind<-vector()
  for (ind in 1:ni)
  {
    dist.row <- ipdij[ind,]
    min.dist.row<-min(dist.row)
    nnj.ind<-c(nnj.ind,j.ind[which(dist.row==min.dist.row)])
    min.dist <- c(min.dist,min.dist.row)
  }
  res<-list(
    nndist=cbind(i.ind,nnj.ind,min.dist),
    base.class=i,
    nn.class=j)
  res
} #end for the function
#'

#################################################################

# funs.kNNdist2cl
#'
#' @title Functions for the \eqn{k^{th}} and \code{k} NN distances
#'
#' @description
#' Two functions: \code{kthNNdist2cl} and \code{kNNdist2cl}.
#'
#' \code{kthNNdist2cl} returns the distances between subjects from class \eqn{i} and their \eqn{k^{th}} NNs from class \eqn{j}.
#' The output is a \code{list} with first entry (\code{kth.nndist}) is an \eqn{n_i \times 3} matrix where \eqn{n_i} is the size of class \eqn{i}
#' and first column is the subject index for class \eqn{i},
#' second column is the index of the \eqn{k^{th}} NN of class \eqn{i} subjects among class \eqn{j} subjects and third column 
#' contains the corresponding \eqn{k^{th}} NN distances. The other entries in the \code{list} are labels of base class and NN class
#' and the value of \code{k}, respectively.
#' 
#' \code{kNNdist2cl} returns the distances between subjects from class \eqn{i} and their \code{k} NNs from class \eqn{j}.
#' The output is a \code{list} with first entry (\code{ind.knndist}) is an \eqn{n_i \times (k+1)} matrix where \eqn{n_i} is the size of class \eqn{i},
#' first column is the indices of class \eqn{i} subjects, 2nd to \eqn{(k+1)}-st  columns are the indices of \code{k} NNs of class \eqn{i}
#' subjects among class \eqn{j} subjects. The second \code{list} entry (\code{knndist}) is an \eqn{n_i \times k} matrix where \eqn{n_i} is the 
#' size of class \eqn{i} and the columns are the \code{k}NN distances of class \eqn{i} subjects to class \eqn{j} subjects. 
#' The other entries in the \code{list} are labels of base class and NN class and the value of \code{k}, respectively. 
#' 
#' The argument \code{within.class.ind} is a logical argument (default=\code{FALSE}) to determine the indexing of the class \eqn{i}
#' subjects. If \code{TRUE}, index numbering of subjects is within the class, from 1 to class size (i.e., \code{1:n_i}), 
#' according to their order in the original data; otherwise, index numbering within class is just the indices
#' in the original data.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param k Integer specifying the number of NNs (of subjects).
#' @param i,j class label of base class and NN classes, respectively.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param within.class.ind A logical parameter (default=\code{FALSE}). If \code{TRUE}, index numbering of subjects is within the class, from 1 to class size (i.e., \code{1:n_i}), 
#' according to their order in the original data; otherwise, index numbering within class is just the indices
#' in the original data.
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return  \code{kthNNdist2cl} returns the \code{list} of elements
#'  \item{kth.nndist}{\eqn{n_i \times 3} matrix where \eqn{n_i} is the size of class \eqn{i}
#' and first column is the subject index for class \eqn{i}, second column is the index of the \code{k}-th NN of class \eqn{i} 
#' subjects among class \eqn{j} subjects and third column contains the corresponding \code{k}-th NN distances,
#' , returned by \code{Zseg.ind.ct} only} 
#'  \item{base.class}{label of base class} 
#'  \item{nn.class}{label of NN class} 
#'  \item{k}{value of \code{k} in \code{k}NN}
#'  
#' \code{kNNdist2cl} returns the \code{list} of elements
#'  \item{ind.knndist}{\eqn{n_i \times (k+1)} matrix where \eqn{n_i} is the size of class \eqn{i}, first column is the indices of class \eqn{i}
#'  subjects, 2nd to \eqn{(k+1)}-st  columns are the indices of \eqn{k} NNs of class \eqn{i} subjects among class \eqn{j} subjects.}
#' \item{knndist}{\eqn{n_i \times k} matrix where \eqn{n_i} is the size of class \eqn{i} and the columns are the \eqn{k}NN distances of class \eqn{i}
#' subjects to class \eqn{j} subjects.}
#'  \item{base.class}{label of base class} 
#'  \item{nn.class}{label of NN class} 
#'  \item{k}{value of \code{k} in \code{k}NN}  
#' 
#' @seealso \code{\link{NNdist2cl}}, \code{\link{kthNNdist}} and \code{\link{kNNdist}}
#' 
#' @name funs.kNNdist2cl
NULL
#'
#' @rdname funs.kNNdist2cl
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for kthNNdist2cl
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#' kthNNdist2cl(ipd,3,1,2,clab)
#' kthNNdist2cl(Y,3,1,2,clab,is.ipd = FALSE)
#' kthNNdist2cl(ipd,3,1,2,clab,within = TRUE)
#'
#' #three class case
#' clab<-sample(1:3,n,replace=TRUE) #class labels
#' table(clab)
#' kthNNdist2cl(ipd,3,2,3,clab)
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#' kthNNdist2cl(ipd,3,1,2,clab) # here kthNNdist2cl(ipd,3,1,12,clab) #gives an error message
#'
#' kthNNdist2cl(ipd,3,"1",2,clab)
#'
#' @export
kthNNdist2cl <- function(x,k,i,j,lab,within.class.ind=FALSE,is.ipd=TRUE,...)
{ 
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  if (any(lab==i)*any(lab==j)==0)
  {stop('given labels i and j are not among the class labels')}
  
  if (within.class.ind==FALSE)
  {
    ns<-length(lab) #sample size of data for which IPDM is computed
    i.ind<-(1:ns)[lab==i]
    j.ind<-(1:ns)[lab==j]
    ni<-length(i.ind)
  } else
  {
    ni<-sum(lab==i)
    nj<-sum(lab==j)
    i.ind<-1:ni
    j.ind<-1:nj 
  }
  
  ipdij<-matrix(ipd[lab==i,lab==j],nrow=ni)
  kth.dist<-kth.nnj.ind<-vector()
  for (ind in 1:ni)
  {
    dist.row <- ipdij[ind,]
    ord<-order(dist.row)
    kth.nnj.ind<-c(kth.nnj.ind,j.ind[ord[k]])
    kth.dist <- c(kth.dist,dist.row[ord[k]] )
  }
  
  res<-list(
    kth.nndist=cbind(i.ind,kth.nnj.ind,kth.dist),
    base.class=i,
    nn.class=j,
    k=k)
  res
} #end for the function
#'
#' @rdname funs.kNNdist2cl
#'
#' @examples
#' #Examples for kNNdist2cl
#' #3D data points
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#' kNNdist2cl(ipd,3,1,2,clab)
#' kNNdist2cl(Y,3,1,2,clab,is.ipd = FALSE)
#'
#' kNNdist2cl(ipd,3,1,2,clab,within = TRUE)
#'
#' #three class case
#' clab<-sample(1:3,n,replace=TRUE) #class labels
#' table(clab)
#' kNNdist2cl(ipd,3,1,2,clab)
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' table(clab)
#'
#' kNNdist2cl(ipd,3,1,2,clab)
#' kNNdist2cl(ipd,3,"1",2,clab) #here kNNdist2cl(ipd,3,"a",2,clab) #gives an error message
#'
#' @export
kNNdist2cl <- function(x,k,i,j,lab,within.class.ind=FALSE,is.ipd=TRUE,...)
{ 
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  if (any(lab==i)*any(lab==j)==0)
  {stop('given labels for i and j are not among the class labels')}
  
  if (within.class.ind==FALSE)
  {
    ns<-length(lab) #sample size of data for which IPDM is computed
    i.ind<-(1:ns)[lab==i]
    j.ind<-(1:ns)[lab==j]
    ni<-length(i.ind)
  } else
  {
    ni<-sum(lab==i)
    nj<-sum(lab==j)
    i.ind<-1:ni
    j.ind<-1:nj 
  }
  
  ipdij<-matrix(ipd[lab==i,lab==j],nrow=ni)
  n<-nrow(ipdij)
  kdist<-nnj.ind<-vector()
  for (ind in 1:n)
  {
    dist.row <- ipdij[ind,]
    ord<-order(dist.row)
    nnj.ind<-rbind(nnj.ind,j.ind[ord[1:k]])
    kdist <- rbind(kdist,dist.row[ord[1:k]] )
  }
  
  res<-list(
    ind.knndist=cbind(i.ind,nnj.ind),
    knndist=kdist,
    base.class=i,
    nn.class=j,
    k=k)
  res
} #end for the function
#'

#################################################################

#' @title Finding the index of the NN of a given point among a subset of points
#'
#' @description
#' Returns the index (indices) of the nearest neighbor(s) of subject \eqn{i} (other than subject \eqn{i}) among the indices of points 
#' provided in the subsample \code{ss} using the given data set or IPD matrix \code{x}. The indices in \code{ss} determine the
#' columns of the IPD matrix to be used in this function. 
#' It will yield a \code{vector} if there are ties, and subject indices correspond to rows (i.e. rows \code{1:n} ) if \code{x} 
#' is the data set and to rows or columns if \code{x} is the IPD matrix.  
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param ss indices of subjects (i.e., row indices in the data set) among with the NN of subject is to be found
#' @inheritParams NN
#'
#' @return Returns a \code{list} with the elements
#'  \item{base.ind}{index of the base subject}
#'  \item{ss.ind}{the index (indices) i.e. row number(s) of the NN of subject \eqn{i} among the subjects with indices
#' provided in \code{ss}}
#'  \item{ss.dis}{distance from subject \eqn{i} to its NN among the subjects in \code{ss}}
#'
#' @seealso \code{\link{NN}} and \code{\link{kNN}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' #3D data points bura
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' #indices of the subsample ss
#' ss<-sample(1:n,floor(n/2),replace=FALSE)
#' NNsub(ss,ipd,2)
#' NNsub(ss,Y,2,is.ipd = FALSE)
#' NNsub(ss,ipd,5)
#'
#' #1D data points
#' n<-15
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' #two class case
#' clab<-sample(1:2,n,replace=TRUE) #class labels
#' #indices of the subsample ss
#' ss<-sample(1:n,floor(n/2),replace=FALSE)
#' NNsub(ss,ipd,2)
#' NNsub(ss,ipd,5)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(60)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' ss<-sample(1:20,10,replace=FALSE) #class labels
#' NNsub(ss,ipd,2)
#' NNsub(ss,ipd,5)
#'
#' @export
NNsub <- function(ss,x,i,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  if (n<=1 || i>n)
  {
  res<-list(base.ind=NA,
         ss.ind=NA,
         ss.dis=NA)
    return(res)
  }
  if (sum(ss==i)==1)
  { 
    ssi<-ss[ss!=i]
    dis.vec<-ipd[i,ssi] #distances from i to points whose indices are in \code{ss}
    min.dis<-min(dis.vec)
    dlen<-sum(dis.vec==min.dis)
    
    ord<-order(dis.vec)
    ss.ind<-ss[ord[1:dlen]]
    D<-dis.vec[ord[1:dlen]]
  } else
  {
    dis.vec<-ipd[i,ss] #distances from i to points whose indices are in \code{ss}
    min.dis<-min(dis.vec)
    dlen<-sum(dis.vec==min.dis)
    
    ord<-order(dis.vec)
    ss.ind<-ss[ord[1:dlen]]
    D<-dis.vec[ord[1:dlen]]
  }
  
  list(base.ind=i,
       ss.ind=ss.ind,
       ss.dis=D)
} #end for the function
#'

################################################################# 

#' @title Nearest Neighbor Contingency Table (NNCT)
#'
#' @description
#' Returns the \eqn{k \times k} NNCT given the IPD matrix or data set \code{x} where \eqn{k} is 
#' the number of classes in the data set.
#' Rows and columns of the NNCT are labeled with the corresponding class labels.
#' 
#' The argument \code{ties} is a logical argument (default=\code{FALSE}) to take ties into account or not.
#' If \code{TRUE} a NN 
#' contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#'
#' See also (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat})
#' and the references therein.
#'
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param ties A logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
#' contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns the \eqn{k \times k} NNCT where \eqn{k} is the number of classes in the data set.
#'
#' @seealso \code{\link{nnct.sub}}, \code{\link{scct}}, \code{\link{rct}}, and \code{\link{tct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#' nnct(ipd,cls,ties = TRUE)
#'
#' nnct(Y,cls,is.ipd = FALSE)
#' nnct(Y,cls,is.ipd = FALSE,method="max")
#' nnct(Y,cls,is.ipd = FALSE,method="mink",p=6)
#'
#' #with one class, it works but really uninformative
#' cls<-rep(1,n)
#' nnct(ipd,cls)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nnct(ipd,fcls)
#'
#' #cls as an unsorted factor
#' fcls1<-sample(c("a","b"),n,replace = TRUE)
#' nnct(ipd,fcls1)
#'
#' fcls2<-sort(fcls1)
#' nnct(ipd,fcls2) #ipd needs to be sorted as well, otherwise this result will not agree with fcls1
#'
#' nnct(Y,fcls1,ties = TRUE,is.ipd = FALSE)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#' nnct(Y,cls,is.ipd = FALSE)
#'
#' #cls as a factor
#' fcls<-rep(letters[1:4],rep(10,4))
#' nnct(ipd,fcls)
#'
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nnct(ipd,fcls)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(3*n)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#' nnct(ipd,cls,ties = TRUE)
#'
#' @export
nnct <- function(x,lab,ties=FALSE,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  ord<-order(lab)#ordering the class label lab first 
  #(to be consistent with row and column labeling in the NNCT)
  lab<-sort(lab) 
  ipd<-ipd[ord,ord]
  flab<-as.factor(lab) #converting class labels to factors
  clnames<-levels(flab)
  k<-length(clnames)
  
  ct<-matrix(0,k,k)
  rownames(ct)<-colnames(ct)<-clnames #row and column names for the NNCT
  if (n<=1)
  {return(ct)}
  
  nlab<-as.numeric(flab)  #converting class labels to numbers
  for(i in 1:n)
  {
    ind <- NN(ipd,i);
    lind<-length(ind)
    for (j in 1:lind)
    {
      addend<-ifelse(ties==FALSE,1,1/lind)
      ct[nlab[i],nlab[ind[j]]]<- ct[nlab[i],nlab[ind[j]]]  + addend
    }
  }
  ct
} #end for the function
#'

#################################################################

#' @title Nearest Neighbor Contingency Table (NNCT) with (only) base points restricted to a subsample
#'
#' @description
#' Returns the \eqn{k \times k} NNCT with (only) base points are restricted to be in the subset of indices \code{ss} using
#' the IPD matrix or data set \code{x} where \eqn{k} is the number of classes in the data set. That is, the base points
#' are the points with indices in \code{ss} but for the NNs the function checks all the points in the data set 
#' (including the points in \code{ss}). 
#' Row and columns of the NNCT are labeled with the corresponding class labels.
#' 
#' The argument \code{ties} is a logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
#' contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' @param ss indices of subjects (i.e., row indices in the data set) chosen to be the base points
#' @inheritParams nnct
#'
#' @return Returns the \eqn{k \times k} NNCT where \eqn{k} is the number of classes in the data set with (only) base points
#' restricted to a subsample \code{ss}.
#'
#' @seealso \code{\link{nnct}} and \code{\link{nnct.boot.dis}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#'
#' #subsampling indices
#' ss<-sample(1:n,floor(n/2))
#' nnct.sub(ss,ipd,cls)
#' nnct.sub(ss,Y,cls,is.ipd = FALSE)
#' nnct.sub(ss,ipd,cls,ties = TRUE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nnct.sub(ss,ipd,fcls)
#'
#' #cls as an unsorted factor
#' fcls<-sample(c("a","b"),n,replace = TRUE)
#' nnct(ipd,fcls)
#' nnct.sub(ss,ipd,fcls)
#'
#' fcls<-sort(fcls)
#' nnct.sub(ss,ipd,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ss<-sample(1:40,30)
#' nnct.sub(ss,ipd,cls)
#'
#' #cls as a factor
#' fcls<-rep(letters[1:4],rep(10,4))
#' nnct.sub(ss,ipd,cls)
#'
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct(ipd,cls)
#'
#' #subsampling indices
#' ss<-sample(1:n,floor(n/2))
#' nnct.sub(ss,ipd,cls)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(120)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ss<-sample(1:40,30)
#' nnct.sub(ss,ipd,cls)
#' nnct.sub(ss,ipd,cls,ties = TRUE)
#'
#' @export
nnct.sub <- function(ss,x,lab,ties=FALSE,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  N<-nrow(ipd) #sample size
  ssvec<-rep(0,N)
  ssvec[ss]<-ss
  
  ord<-order(lab)
  lab<-sort(lab) #ordering the class label lab first 
  #(to be consistent with row and column labeling in the NNCT)
  ipd<-ipd[ord,ord]
  ssvec<-ssvec[ord]
  ss<-ssvec[ssvec>0]
  
  flab<-as.factor(lab) #converting class labels to factors
  clnames<-levels(flab)
  k<-length(clnames)
  
  nss<-length(ss) #length of the indices to select for NNCT construction
  
  ct<-matrix(0,k,k)
  rownames(ct)<-colnames(ct)<-clnames #row and column names for the NNCT
  if (nss<=1)
  {return(ct)}
  
  nlab<-as.numeric(flab)  #converting class labels to numbers
  for(i in 1:nss)
  {
    ind <- NN(ipd,ss[i]);
    lind<-length(ind)
    for (j in 1:lind)
    {
      addend<-ifelse(ties==FALSE,1,1/lind)
      ct[nlab[ss[i]],nlab[ind[j]]]<- ct[nlab[ss[i]],nlab[ind[j]]]  + addend
    }
  }
  ct
} #end for the function
#'

#################################################################

#' @title Bootstrap Nearest Neighbor Contingency Table (NNCT)
#' 
#' @description
#' Returns the \eqn{k \times k} NNCT with sampling replacement of the points for each base point. That is, for each base 
#' point, the rows in the IPD matrix are sampled with replacement and the NN counts are updated accordingly.
#' Row and columns of the NNCT are labeled with the corresponding class labels.
#' 
#' The argument self is a logical argument (default=\code{TRUE}) for including the base point in the resampling or not.
#' If \code{TRUE}, for each base point all entries in the row are sampled (with replacement) so the point itself can
#' also be sampled multiple times and if \code{FALSE} the point is excluded from the resampling (i.e. other points
#' are sampled with replacement).
#' 
#' The argument \code{ties} is a logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
#' contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#'
#' @inheritParams nnct
#' @param self A logical argument (default=\code{TRUE}). If \code{TRUE}, for each base point, all entries in the 
#' row are sampled (with replacement) and if \code{FALSE} the point is excluded from the resampling (i.e. other points
#' are sampled with replacement).
#'
#' @return Returns the \eqn{k \times k} NNCT where \eqn{k} is the number of classes in the data set with sampling replacement
#' of the rows of the IPD matrix.
#'
#' @seealso \code{\link{nnct}} and \code{\link{nnct.sub}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct.boot.dis(ipd,cls)
#' nnct.boot.dis(Y,cls,is.ipd = FALSE) #may give different result from above due to random sub-sampling
#' nnct.boot.dis(ipd,cls,self = FALSE)
#' nnct.boot.dis(ipd,cls,ties = FALSE) #differences are due to ties and resampling of distances
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nnct.boot.dis(ipd,fcls)
#'
#' #cls as an unsorted factor
#' fcls<-sample(c("a","b"),n,replace = TRUE)
#' nnct.boot.dis(ipd,fcls)
#'
#' fcls<-sort(fcls)
#' nnct.boot.dis(ipd,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct.boot.dis(ipd,cls)
#'
#' #cls as a factor
#' fcls<-rep(letters[1:4],rep(10,4))
#' nnct.boot.dis(ipd,fcls)
#'
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct.boot.dis(ipd,cls)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nnct.boot.dis(ipd,fcls)
#'
#' #with possible ties in the data
#' Y<-matrix(round(runif(3*n)*10),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' nnct.boot.dis(ipd,cls)
#' nnct.boot.dis(ipd,cls,self = FALSE)
#' nnct.boot.dis(ipd,cls,ties = FALSE) #differences are due to ties and resampling of distances
#'
#' @export 
nnct.boot.dis <- function(x,lab,self=TRUE,ties=TRUE,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  
  n<-nrow(ipd)
  ord<-order(lab)#ordering the class label lab first 
  #(to be consistent with row and column labeling in the NNCT)
  
  lab<-sort(lab) 
  ipd<-ipd[ord,ord]
  
  flab<-as.factor(lab) #converting class labels to factors
  clnames<-levels(flab)
  k<-length(clnames)
  
  ct<-matrix(0,k,k)
  rownames(ct)<-colnames(ct)<-clnames #row and column names for the NNCT
  if (n<=1)
  {return(ct)}
  
 # pr<-rep(1,n)
  
  nlab<-as.numeric(flab)  #converting class labels to numbers
  for(i in 1:n)
  {
    pr<-rep(1,n)
    if (self==FALSE) { pr[i]<-0 }
    ss<-sample(1:n,replace=TRUE,prob=pr) # sampling with replacement of column indices of the ipd
    ind <- NNsub(ss,ipd,i)$ss.ind;
    lind<-length(ind)
    for (j in 1:lind)
    {
      addend<-ifelse(ties==FALSE,1,1/lind)
      ct[nlab[i],nlab[ind[j]]]<- ct[nlab[i],nlab[ind[j]]]  + addend
    } 
  }
  ct
} #end for the function
#'

#################################################################

# funsOnevsRest
#'
#' @title Functions for one versus rest type labeling
#'
#' @description
#' Two functions: \code{lab.onevsrest} and \code{classirest}.
#'
#' Both functions relabel the points, keeping class \eqn{i} label as is and relabeling the other classes as "rest".
#' Used in the one-vs-rest type comparisons after the overall segregation test is found to be significant.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#' 
#' @param i label of the class that is to be retained in the post-hoc comparison. 
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' 
#' @return Both functions return the data relabeled as class \eqn{i} label is retained and the remaining is
#' relabeled as "rest".
#' 
#' @seealso \code{\link{pairwise.lab}}
#' 
#' @name funsOnevsRest
NULL
#'
#' @rdname funsOnevsRest
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' lab.onevsrest(1,cls)
#' classirest(2,cls)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' lab.onevsrest("a",fcls)
#' lab.onevsrest("b",fcls)
#' classirest("b",fcls)
#'
#' #cls as a factor
#' fcls<-rep(letters[1:4],rep(10,4))
#' lab.onevsrest("b",fcls)
#' classirest("b",fcls)
#'
#' @export
lab.onevsrest <- function(i,lab)
{
  if (any(lab==i)==FALSE)
  {stop('given label i is not among the class labels')}
  
  lab[lab!=i]<-"rest"
  lab
} #end for the function
#' 
#' @rdname funsOnevsRest
#'
#' @export
classirest <- function(i,lab)
{
  if (any(lab==i)==FALSE)
  {stop('given label i is not among the class labels')}
  
  lab[lab!=i]<-"rest"
  lab
} #end for the function
#'

#################################################################

#' @title Keeping the pair of the specified labels in the data
#'
#' @description
#' Keeps only the specified labels \eqn{i} and \eqn{j} and returns the data from classes with these labes and also
#' the corresponding label vector having class labels \eqn{i} and \eqn{j} only.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param i,j Label of the classes that are to be retained in the post-hoc comparison.
#'
#' @return A \code{list} with two elements
#' \item{data.pair}{The type of the pattern from which points are to be generated}
#' \item{lab.pair}{The \code{"main"} title for the plot of the point pattern}
#' 
#' @seealso \code{\link{lab.onevsrest}} and \code{\link{classirest}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' pairwise.lab(Y,cls,1,2)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' pairwise.lab(Y,cls,2,3)
#'
#' #cls as a factor
#' fcls<-rep(letters[1:4],rep(10,4))
#' pairwise.lab(Y,fcls,"b","c")
#'
#' @export
pairwise.lab <- function(dat,lab,i,j)
{
  if (any(lab==i)*any(lab==j)==0)
  {stop('at least one of the given labels for i and j are not among the class labels')}
  
  ind<-(lab==i | lab == j)
  pair.dat <-dat[ind,];
  cl<-lab[ind]
  
  list(data.pair=pair.dat,
       lab.pair=cl)
} #end for the function
#'

#################################################################

# funsRowColSums
#'
#' @title Functions for row and column sums of a matrix
#'
#' @description
#' Two functions: \code{row.sum} and \code{col.sum}.
#'
#' \code{row.sum} returns the row sums of a given matrix (in particular a contingency table) as a vector and 
#' \code{col.sum} returns the column sums of a given matrix as a vector. \code{row.sum} is equivalent to 
#' \code{\link[base]{rowSums}} function and \code{col.sum} is equivalent to \code{\link[base]{colSums}}
#' function in the \code{base} package.
#'   
#' @param ct A matrix, in particular a contingency table
#' 
#' @return 
#' \code{row.sum} returns the row sums of \code{ct} as a vector
#' \code{col.sum} returns the column sums of \code{ct} as a vector
#' 
#' @seealso \code{\link[base]{rowSums}} and \code{\link[base]{colSums}}
#' 
#' @name funsRowColSums
NULL
#'
#' @rdname funsRowColSums
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' row.sum(ct)
#' rowSums(ct)
#'
#' col.sum(ct)
#' colSums(ct)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' row.sum(ct)
#' rowSums(ct)
#'
#' col.sum(ct)
#' colSums(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' row.sum(ct)
#' rowSums(ct)
#'
#' col.sum(ct)
#' colSums(ct)
#'
#' @export
row.sum <- function(ct)
{
  apply(ct,1,sum)
} #end for the function
#' 
#' @rdname funsRowColSums
#'
#' @export
col.sum <- function(ct)
{
  apply(ct,2,sum)
} #end for the function
#'

#################################################################

#' @title Dixon's Segregation Indices for NNCTs
#'
#' @description
#' Returns Dixon's segregation indices in matrix form based on entries of the NNCT, \code{ct}. 
#' Segregation index for cell \eqn{i,j} is defined as \eqn{log(N_{ii}(n-n_i)/((n_i-N_{ii})(n_i-1))} if \eqn{i=j}
#' and
#' as \eqn{log(N_{ij}(n-n_j-1)/((n_i-N_{ij})(n_j))} if \eqn{i \ne j}. 
#' See (\insertCite{dixon:NNCTEco2002,ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' 
#' The argument \code{inf.corr} is a logical argument (default=\code{FALSE}) to avoid \eqn{\pm \infty} for the segregation
#' indices. If \code{TRUE} indices are modified so that they are finite and if \code{FALSE} the above definition is used. 
#' (See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).
#'
#' @param ct A contingency table, in particular an NNCT
#' @param inf.corr A logical argument (default=\code{FALSE}). If \code{TRUE}, indices are modified so that 
#' they are finite and if \code{FALSE} the above definition in the description is used.
#' 
#' @return Returns a \code{matrix} of segregation indices which is of the same dimension as \code{ct}.
#'
#' @seealso \code{\link{Pseg.coeff}}, \code{\link{seg.coeff}}, \code{\link{Zseg.ind}}
#' and \code{\link{Zseg.ind.ct}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#' seg.ind(ct)
#' seg.ind(ct,inf.corr = TRUE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' seg.ind(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' seg.ind(ct)
#' seg.ind(ct,inf.corr = TRUE)
#'
#' ct<-matrix(c(0,10,5,5),ncol=2)
#' seg.ind(ct)
#'
#' seg.ind(ct,inf.corr = TRUE)
#'
#' @export
seg.ind <- function(ct,inf.corr=FALSE)
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct)
  
  Sdix<- matrix(0,k,k);
  if (inf.corr==FALSE) #without infinity index correction
  {
    for (i in 1:k)
      for (j in 1:k)
      {
        if (i == j)
          Sdix[i,i]<- log((ct[i,i]/(rs[i]-ct[i,i]))/((rs[i]-1)/(n-rs[i]))) 
        else 
          Sdix[i,j]<- log((ct[i,j]/(rs[i]-ct[i,j]))/(rs[j]/(n-rs[j]-1)));
      }
  } else  #with infinity index correction
  {
    for (i in 1:k)
      for (j in 1:k)
      {
        if (i == j)
          Sdix[i,i]<- log(((ct[i,i]+1)/(rs[i]-ct[i,i]+1))/( (rs[i]*(rs[i]-1)+n-1)/(rs[i]*(n-rs[i])+n-1) )) 
        else 
          Sdix[i,j]<- log(((ct[i,j]+1)/(rs[i]-ct[i,j]+1))/( (rs[i]*rs[j]+n-1)/(rs[i]*(n-rs[j]-1)+n-1) ));
      }  
  }
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(Sdix)<-colnames(Sdix)<-clnames #row and column names for the matrix
  
  if (any(is.na(Sdix)))
  { warning("some of the segregation indices are NaN, 
  so use inf.corr=TRUE in the function call!")
    return(Sdix)}
  
  if( (any(abs(Sdix)==Inf)))
  {warning("some of the segregation indices are +- infinity, 
  so use inf.corr=TRUE in the function call!")}
  
  Sdix
} #end for the function
#'

#################################################################

# funsZsegind
#'
#' @title Z Tests for Segregation Indices
#'
#' @description
#' Two functions: \code{Zseg.ind.ct} and \code{Zseg.ind}.
#'
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' segregation indices from their expected values under RL or CSR for each segregation index in the NNCT.
#' The test for each cell \eqn{i,j} is based on the normal approximation of the corresponding segregation index.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, lower and upper confidence levels, sample estimates (i.e. observed values) and null value(s) (i.e. expected values) for the segregation indices
#' and also names of the test statistics, estimates, null value and the method and the data set used.
#'
#' The null hypothesis for each cell \eqn{i,j} is that the corresponding segregation index equal to the expected value
#' (which is 0) under RL or CSR.
#'
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Zseg.ind.ct} only 
#' @param varN The variance matrix for cell counts in the NNCT, \code{ct} ; used in \code{Zseg.ind.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zseg.ind} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zseg.ind} only
#' @param inf.corr A logical argument (default=\code{FALSE}). If \code{TRUE}, indices are modified so that 
#' they are finite and if \code{FALSE} the above definition in the description is used.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, for the segregation indices
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Zseg.ind} only
#'
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of test statistics for the segregation indices}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the segregation indices at the given confidence
#' level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits of the segregation indices,
#' provided in \code{conf.level}.}
#' \item{estimate}{Estimate of the parameter, i.e., matrix of the observed segregation indices}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function}
#' \item{null.value}{Hypothesized values for the parameters, i.e. the null values of the segregation indices, 
#' which are all 0 under RL or CSR.}
#' \item{null.name}{Name of the null value}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zseg.ind.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zseg.ind} only}
#' 
#' @seealso \code{\link{seg.ind}} and \code{\link{Zseg.coeff}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZsegind
NULL
#'
#' @rdname funsZsegind
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' seg.ind(ct)
#' seg.ind(ct,inf.corr=TRUE)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' varN
#'
#' Zseg.ind(Y,cls)
#' Zseg.ind(Y,cls,inf.corr=TRUE)
#' Zseg.ind.ct(ct,varN)
#'
#' Zseg.ind(Y,cls,alt="g")
#' Zseg.ind.ct(ct,varN,alt="g")
#'
#' Zseg.ind(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Zseg.ind(Y,cls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' varN
#'
#' Zseg.ind(Y,cls)
#' Zseg.ind(Y,cls,inf.corr = TRUE)
#'
#' Zseg.ind.ct(ct,varN)
#' Zseg.ind.ct(ct,varN,inf.corr = TRUE)
#'
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' Zseg.ind(X,cls)
#' Zseg.ind.ct(ct,varN)
#' Zseg.ind.ct(ct,varN,inf.corr=TRUE)
#'
#' @export
Zseg.ind.ct <- function(ct,varN,inf.corr=FALSE,
                        alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-dsi<-seg.ind(ct,inf.corr=inf.corr) #inf.corr takes care of infinity correction in the indices
  estimate.name <-"Dixon's segregation indices"
  
  if (all(is.na(dsi)) || all(is.na(varN)))
  {stop('All of the segregation indices or all entries of the variance, varN, are NaN, so the tests are not defined')}
  
  null.seg.ind<-0
  names.null <-"(expected value of each of the) Dixon's segregation indices"
  
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(rs)
  
  ts<-stderrS<- matrix(0,k,k);
  if (inf.corr==FALSE) #without infinity index correction
  {
    for (i in 1:k)
      for (j in 1:k)
      {
        if (i == j)
        { stderrS[i,i]<-(sqrt(varN[i,i])*(n-1)^2)/(rs[i]*(n-rs[i])*(rs[i]-1))
        ts[i,i]<- dsi[i,i]/stderrS[i,i] 
        } else 
        {
          stderrS[i,j]<-(sqrt(varN[i,j])*(n-1)^2)/(rs[i]*rs[j]*(n-rs[j]-1))
          ts[i,j]<- dsi[i,j]/stderrS[i,j] 
        }
      }
    #for +infinity correction
    for (i in 1:k)
      for (j in 1:k)
      {
        if (!is.na(dsi[i,j]) & dsi[i,j]==Inf)
          ts[i,j]<- 10^10
      }
    #for -infinity correction
    for (i in 1:k)
      for (j in 1:k)
      {
        if (!is.na(dsi[i,j]) & dsi[i,j]==-Inf)
          ts[i,j]<- -10^10
      }
  } else  #with infinity index correction
  {
    for (i in 1:k)
      for (j in 1:k)
      {
        if (i == j)
        { stderrS[i,i]<-(sqrt(varN[i,i])*(rs[i]+2)*(n-1)^2)/((rs[i]*(rs[i]-1)+n-1)*(rs[i]*(n-rs[i])+n-1))
        ts[i,i]<- dsi[i,i]/stderrS[i,i] 
        } else 
        {
          stderrS[i,j]<-(sqrt(varN[i,j])*(rs[i]+2)*(n-1)^2)/((rs[i]*rs[j]+n-1)*(rs[i]*(n-rs[j]-1)+n-1))
          ts[i,j]<- dsi[i,j]/stderrS[i,j] 
        }
      }
  }
  
  if (all(is.na(ts)))
  {stop('All of the test stat statistics are NaN, the test are not well defined')}
  
  if (alternative == "less") {
    pval <-pnorm(ts)
    lcl <-NULL
    ucl <-estimate+qnorm(conf.level)*stderrS
  }
  else if (alternative == "greater") {
    pval <-pnorm(ts, lower.tail = FALSE)
    ucl <-NULL
    lcl <-estimate-qnorm(conf.level)*stderrS
  }
  else {
    pval <-2 * pnorm(-abs(ts))
    alpha <-1 - conf.level
    crit.val <-qnorm(1-alpha/2)
    lcl <-estimate-crit.val*stderrS
    ucl <-estimate+crit.val*stderrS
  }
  
  cnf.lvl<-conf.level
  
  ifelse(inf.corr==FALSE,
         method <-c("Large Sample z-Test for Dixon's Segregation Indices"),
         method <-c("Large Sample z-Test for Dixon's Segregation Indices\n
                    with correction for zero cell counts"))
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  if (!is.null(lcl)) {rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) {rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"standardized segregation indices (i.e., z-scores)"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name, #this is for other functions to have a different description for the sample estimates
    null.value = null.seg.ind,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZsegind
#'
#' @export
Zseg.ind <- function(dat,lab,inf.corr=FALSE,
                     alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  
  res<- Zseg.ind.ct(ct,varN,inf.corr=inf.corr,alternative=alternative,conf.level=conf.level)
  
  dname <-deparse(substitute(dat))
  res$data.name<-dname
  
  return(res)
} #end for the function
#'

#################################################################

#' @title Expected Values of the Cell Counts in NNCT
#'
#' @description Returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the expected cell counts of
#' the NNCT under RL or CSR. The class sizes given as the row sums of \code{ct} and the row and column names are
#' inherited from \code{ct}.
#' 
#' See also (\insertCite{dixon:1994,ceyhan:eest-2010;textual}{nnspat}).
#'
#' @param ct A nearest neighbor contingency table
#'
#' @return A \code{matrix} of the expected values of cell counts in the NNCT.
#'
#' @seealso \code{\link{nnct}} and \code{\link{EV.tct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.nnct(ct)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#' EV.nnct(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.nnct(ct)
#'
#' ct<-matrix(c(0,10,5,5),ncol=2)
#' EV.nnct(ct)
#'
#' @export
EV.nnct <- function(ct)
{
  rs <- row.sum(ct)
  
  k<-length(rs);
  n<-sum(rs) 
  
  EN<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    {
      if (i == j)
        EN[i,j]<- rs[i]*(rs[i]-1)/(n-1)  
      else 
        EN[i,j]<- rs[i]*rs[j]/(n-1)
    }
  clnames<-colnames(ct)
  rownames(EN)<-colnames(EN)<-clnames #row and column names from the NNCT
  EN
} #end for the function
#'

#################################################################

#' @title Expected Values of the Type I cell-specific tests
#'
#' @description Returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the expected values
#' of the Type I cell-specific test statistics, \eqn{T^I_{ij}}. 
#' The row and column names are inherited from \code{ct}. 
#' These expected values are valid under RL or CSR.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}) and the references therein.
#' 
#' @param ct A nearest neighbor contingency table
#' 
#' @return A \code{matrix} of the expected values of Type I cell-specific tests.
#'
#' @seealso \code{\link{EV.tct}}, \code{\link{tct}} and \code{\link{EV.nnct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.tctI(ct)
#' 
#' @export
EV.tctI <- function(ct)
{
  rs <- row.sum(ct); 
  k<-length(rs);
  n<-sum(rs) 
  
  expT<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    {
      if (i == j)
        expT[i,j]<- rs[i]*(rs[i]-n)/(n*(n-1))  
      else 
        expT[i,j]<- rs[i]*rs[j]/(n*(n-1))
    }
  clnames<-colnames(ct)
  rownames(expT)<-colnames(expT)<-clnames #row and column names from the NNCT
  expT
} #end for the function
#'

#################################################################

#' @title Expected Values of the Types I-IV cell-specific tests
#'
#' @description Returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the expected values
#' of the \eqn{T_{ij}} values which are the Types I-IV cell-specific test statistics (i.e., \eqn{T^I_{ij}-T^{IV}_{ij}})
#' under RL or CSR. 
#' The row and column names are inherited from \code{ct}. The type argument specifies the type
#' of the cell-specific test among the types I-IV tests.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}) and the references therein.
#' 
#' @param ct A nearest neighbor contingency table
#' @param type The type of the cell-specific test, default=\code{"III"}. Takes on values \code{"I"}-\code{"IV"} (or 
#' equivalently \code{1-4}, respectively.
#' 
#' @return A \code{matrix} of the expected values of Type I-IV cell-specific tests.
#'
#' @seealso \code{\link{EV.tctI}}, \code{\link{tct}} and \code{\link{EV.nnct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.tct(ct,2)
#' EV.tct(ct,"II")
#' EV.tctI(ct)
#' 
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#' EV.tct(ct,2)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.tct(ct,2)
#'
#' ct<-matrix(c(0,10,5,5),ncol=2)
#' EV.tct(ct,2)
#'
#' @export
EV.tct <- function(ct,type="III")
{
  ET<-  switch(type,
               I = { ET<- EV.tctI(ct) },
               II = { ET<- EV.tctI(ct) },
               III = { k<-nrow(ct)
               ET<- matrix(0,k,k) },
               IV = { k<-nrow(ct)
               ET<- matrix(0,k,k) }
  )
  if (is.null(ET)) stop("Enter numbers 1-4 or I-IV in quotes for type")

  clnames<-colnames(ct)
  rownames(ET)<-colnames(ET)<-clnames #row and column names from the nnct
  ET
} #end for the function
#'

#################################################################

#' @title \eqn{T} Contingency Table (TCT)
#'
#' @description Returns the \code{T} contingency table (TCT), which is a matrix of same dimension as, \code{ct}, 
#' whose entries are the values of the Types I-IV cell-specific test statistics, \eqn{T^I_{ij}-T^{IV}_{ij}}. 
#' The row and column names are inherited from \code{ct}. The type argument specifies the type
#' of the cell-specific test among the types I-IV tests. 
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}) and the references therein.
#' 
#' @inheritParams EV.tct
#' 
#' @return A \code{matrix} of the values of Type I-IV cell-specific tests
#'
#' @seealso \code{\link{cellsTij}} and \code{\link{nnct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' type.lab<-c("I","II","III","IV")
#' for (i in 1:4)
#' { print(paste("T_ij values for cell specific tests for type",type.lab[i]))
#'   print(tct(ct,i))
#' }
#'
#' tct(ct,"II")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#' tct(ct,2)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' tct(ct,2)
#'
#' ct<-matrix(c(0,10,5,5),ncol=2)
#' tct(ct,2)
#'
#' @export 
tct <- function(ct,type="III") 
{
  rs <- row.sum(ct); 
  cs <- col.sum(ct); 
  
  k<-nrow(ct);
  n<-sum(ct) 
  
  cells<- matrix(0,k,k);
  cells <- switch(type,
         I = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               cells[i,j]<- ct[i,j]-(rs[i]*cs[j])/n 
             } 
           cells},
         II = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               cells[i,j]<- ct[i,j]-(rs[i]*rs[j])/n 
             } 
           cells},
         III = {   
           for (i in 1:k)
             for (j in 1:k)
             {
               if (i == j)
                 cells[i,j]<- ct[i,j]-(rs[i]-1)*cs[j]/(n-1)  
               else 
                 cells[i,j]<- ct[i,j]-rs[i]*cs[j]/(n-1)
             } 
           cells},
         IV = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               if (i == j)
                 cells[i,j]<- (rs[i]/n)* ((n-1)*ct[i,j]/(rs[i]-1)-cs[j])  
               else 
                 cells[i,j]<- (1/n)*((n-1)*ct[i,j]-rs[i]*cs[j])
             } 
           cells}
  )
  
  if (is.null(cells)) stop("Enter numbers 1-4 or I-IV in quotes for type")
  
  clnames<-colnames(ct)
  rownames(cells)<-colnames(cells)<-clnames #row and column names from the NNCT 
  
  cells
} #end for the function
#'

################################################################# 

#' @title Entries for the Types I-IV cell-specific tests
#'
#' @description Returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the values
#' of the Types I-IV cell-specific test statistics, \eqn{T^I_{ij}-T^{IV}_{ij}}. 
#' The row and column names are inherited from \code{ct}. The type argument specifies the type
#' of the cell-specific test among the types I-IV tests. 
#' Equivalent to the function \code{\link{tct}} in this package.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}) and the references therein.
#' 
#' @inheritParams EV.tct
#' 
#' @return A \code{matrix} of the values of Type I-IV cell-specific tests
#'
#' @seealso \code{\link{tct}} and \code{\link{nnct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' type.lab<-c("I","II","III","IV")
#' for (i in 1:4)
#' { print(paste("T_ij values for cell specific tests for type",type.lab[i]))
#'   print(cellsTij(ct,i))
#' }
#'
#' cellsTij(ct,"II")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#' cellsTij(ct,2)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' cellsTij(ct,2)
#'
#' ct<-matrix(c(0,10,5,5),ncol=2)
#' cellsTij(ct,2)
#'
#' @export 
cellsTij <- function(ct,type="III")  
{
  k<-nrow(ct);
  rs<-row.sum(ct)
  cs<- col.sum(ct)
  n<-sum(ct)
  
  cells<- matrix(0,k,k);
  cells <- switch(type,
         I = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               cells[i,j]<- ct[i,j]-(rs[i]*cs[j])/n 
             } 
           cells},
         II = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               cells[i,j]<- ct[i,j]-(rs[i]*rs[j])/n 
             } 
           cells},
         III = {   
           for (i in 1:k)
             for (j in 1:k)
             {
               if (i == j)
                 cells[i,j]<- ct[i,j]-(rs[i]-1)*cs[j]/(n-1)  
               else 
                 cells[i,j]<- ct[i,j]-rs[i]*cs[j]/(n-1)
             } 
           cells},
         IV = { 
           for (i in 1:k)
             for (j in 1:k)
             {
               if (i == j)
                 cells[i,j]<- (rs[i]/n)* ((n-1)*ct[i,j]/(rs[i]-1)-cs[j])  
               else 
                 cells[i,j]<- (1/n)*((n-1)*ct[i,j]-rs[i]*cs[j])
             } 
           cells}
  )
  
  if (is.null(cells)) stop("Enter numbers 1-4 or I-IV in quotes for type")
  
  clnames<-colnames(ct)
  rownames(cells)<-colnames(cells)<-clnames #row and column names from the NNCT 
  
  cells
} #end for the function
#'

#################################################################

# funs.pijPij
#'
#' @title The functions for probability of selecting a number of points from respective classes
#' 
#' @description
#' The ancillary probability functions used in computation of the variance-covariance matrices
#' of various NN spatial tests such as NNCT tests and tests based on other contingency tables.
#' These functions can be classified as \code{pij} and \code{Pij} type functions. The \code{pij} functions are for individual 
#' probabilities and the corresponding \code{Pij} functions are the summed \code{pij} values. For example \eqn{p_{iijk}}
#' is the probability of any 4 points with 2 from class \eqn{i}, and others are from classes \eqn{j} and \eqn{k}. 
#' These probabilities are for data from RL or CSR.
#'
#' @param n A positive integer representing the size of the data set (i.e., number of observations in the data
#' set).
#' @param k,l,m,p Positive integers, usually representing the class sizes, used in \code{pij} type functions only.
#' Number of these arguments required depends on the number of distinct indices of \eqn{p}, e.g. \eqn{p_{ij}} requires
#' \code{k,l,n} and \eqn{p_{iijk}} requires \code{k,l,m,n} as input.
#' @param nvec A \code{vector} ofpositive integers representing the sizes of classes in the data set, used in 
#' \code{Pij} type functions only.
#'
#' @return Probability values for the selected points being from the indicated classes.
#'
#' @name funs.pijPij
NULL
#'
#' @seealso \code{\link{pk}}
#'
#' @rdname funs.pijPij
#'
p11 <- function(k,n)
{
  k*(k-1)/(n*(n-1))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P11 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp11<-0
  for (i in 1:k)
  {
    sp11<-sp11+p11(nvec[i],n)
  }
  sp11
} #end for the function
#'
#' @rdname funs.pijPij
#'
p12 <- function(k,l,n)
{
  k*l/(n*(n-1))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P12 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp12<-0
  for (i in 1:(k-1))
    for (j in (i+1):k)
    {
      sp12<-sp12+p12(nvec[i],nvec[j],n)
    }
  2*sp12
} #end for the function
#'
#' @rdname funs.pijPij
#'
p111 <- function(k,n)
{
  k*(k-1)*(k-2)/(n*(n-1)*(n-2))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P111 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp111<-0
  for (i in 1:k)
  {
    sp111<-sp111+p111(nvec[i],n)
  }
  sp111
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1111 <- function(k,n)
{
  k*(k-1)*(k-2)*(k-3)/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P1111 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp1111<-0
  for (i in 1:k)
  {
    sp1111<-sp1111+p1111(nvec[i],n)
  }
  sp1111
} #end for the function
#'
#' @rdname funs.pijPij
#'
p112 <- function(k,l,n)
{
  k*(k-1)*l/(n*(n-1)*(n-2))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P112 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp112<-0
  for (i in 1:k)
  {
    jind<-(1:k)[-i]
    for (j in jind)
    {
      sp112<-sp112+p112(nvec[i],nvec[j],n)
    }
  }
  sp112
} #end for the function
#'
#' @rdname funs.pijPij
#'
p122 <- function(k,l,n)
{
  k*l*(l-1)/(n*(n-1)*(n-2))
} #end for the function
#'
#' @rdname funs.pijPij
#'
p123 <- function(k,l,m,n)
{
  k*l*m/(n*(n-1)*(n-2))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P123 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp123<-0
  for (i in 1:(k-2))
  {
    for (j in (i+1):(k-1))
    {
      for (l in (j+1):k)
      {
        sp123<-sp123+p123(nvec[i],nvec[j],nvec[l],n)
      }
    }
  }
  6*sp123
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1234 <- function(k,l,m,p,n)
{
  k*l*m*p/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P1234 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp1234<-0
  for (i in 1:(k-3))
  {
    for (j in (i+1):(k-2))
    {
      for (l in (j+1):(k-1))
      {
        for (m in (l+1):k)
        {
          sp1234<-sp1234+p1234(nvec[i],nvec[j],nvec[l],nvec[m],n)
        }
      }
    }
  }
  24*sp1234
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1223 <- function(k,l,m,n)
{
  k*l*(l-1)*m/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1123 <- function(k,l,m,n)
{
  k*(k-1)*l*m/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P1123 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp1123<-0
  for (i in 1:k)
  {
    jind<-(1:k)[-i]
    for (j in jind)
    {
      lind<-(1:k)[-c(i,j)]
      for (l in lind)
      {
        sp1123<-sp1123+p1123(nvec[i],nvec[j],nvec[l],n)
      }
    }
  }
  sp1123
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1122 <- function(k,l,n)
{
  k*(k-1)*l*(l-1)/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P1122 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp1122<-0
  for (i in 1:(k-1))
    for (j in (i+1):k)
    {
      sp1122<-sp1122+p1122(nvec[i],nvec[j],n)
    }
  2*sp1122
} #end for the function
#'
#' @rdname funs.pijPij
#'
p1112 <- function(k,l,n)
{
  k*(k-1)*(k-2)*l/(n*(n-1)*(n-2)*(n-3))
} #end for the function
#'
#' @rdname funs.pijPij
#'
P1112 <- function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  sp1112<-0
  for (i in 1:k)
  {
    jind<-(1:k)[-i]
    for (j in jind)
    {
      sp1112<-sp1112+p1112(nvec[i],nvec[j],n)
    }
  }
  sp1112
} #end for the function
#'

#################################################################

#' @title Probability of \code{k} items selected from the class with size \eqn{n_1}
#'
#' @description
#' Returns the ratio \eqn{n_1(n_1-1) \cdots (n_1-(k-1))/(n(n-1) \cdots (n-(k-1))}, 
#' which is the probability that the \code{k} selected
#' objects are from class 1 with size \eqn{n_1} (denoted as \code{n1} as an argument)
#' and the total data size is \code{n}.
#' This probability is valid under RL or CSR.
#' 
#' @param n A positive integer representing the size of the data set 
#' (i.e., number of observations in the data set).
#' @param n1 A positive integer, representing the class size of interest which has size \eqn{n_1}
#' @param k Number of items selected from the data set
#' 
#' @return Returns the probability of \code{k} items selected from \code{n} items are from the class of interest
#'(i.e., from the class whose size is \eqn{n_1})
#'
#' @seealso \code{\link{p11}} and \code{\link{p12}} etc.
#'
pk <- function(n,n1,k)
{
  r<-prod(n1:(n1-(k-1)))/prod(n:(n-(k-1)))
  r
}
#'

#################################################################

#' @title Variances of Cell Counts in an NNCT
#'
#' @description Returns the variances of cell counts \eqn{N_{ij}} for \eqn{i,j=1,\ldots,k} in the NNCT, \code{ct} in matrix form which
#' is of the same dimension as \code{ct}. These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @param ct A nearest neighbor contingency table
#' @param Q The number of shared NNs
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs)
#'
#' @return A \code{matrix} of same dimension as, \code{ct}, whose entries are the variances of the cell counts 
#' in the NNCT with class sizes given as the row sums of \code{ct}. The row and column names are inherited from \code{ct}.
#'
#' @seealso \code{\link{var.tct}}, \code{\link{var.nnsym}} and \code{\link{cov.nnct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' var.nnct(ct,Qv,Rv)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#' var.nnct(ct,Qv,Rv)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' var.nnct(ct,Qv,Rv)
#'
#' @export
var.nnct <- function(ct,Q,R) 
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct) 
  
  var<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    {
      if (i == j)
      {
        Paa <- p11(rs[i],n); Paaa <- p111(rs[i],n); Paaaa <- p1111(rs[i],n)
        var[i,j]<- (n+R)*Paa+(2*n-2*R+Q)*Paaa+(n^2-3*n-Q+R)*Paaaa-(n*Paa)^2
      }  
      else 
      {
        Pab <- p12(rs[i],rs[j],n); Paab <- p112(rs[i],rs[j],n); Paabb <- p1122(rs[i],rs[j],n)
        var[i,j]<- n*Pab+Q*Paab+(n^2-3*n-Q+R)*Paabb-(n*Pab)^2
      }
    }
  
  clnames<-colnames(ct)
  rownames(var)<-colnames(var)<-clnames #row and column names from the NNCT
  var
} #end for the function
#'

#################################################################

# funsZcell.nnct
#'
#' @title Dixon's Cell-specific Z Tests of Segregation for NNCT
#'
#' @description
#' Two functions: \code{Zcell.nnct.ct} and \code{Zcell.nnct}.
#'
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' cell counts from the expected values under RL or CSR for each cell (i.e., entry) in the NNCT.
#' The test for each cell \eqn{i,j} is based on the normal approximation of the corresponding cell count, \eqn{N_{ij}}
#' and are due to \insertCite{dixon:1994,dixon:NNCTEco2002;textual}{nnspat}.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e., null value(s)), lower and upper confidence levels, sample estimates (i.e. observed values)
#' for the cell counts and also names of the test statistics, estimates, null values and the method and
#' the data set used.
#' 
#' The null hypothesis for each cell \eqn{i,j} is that the corresponding cell count is equal to the expected value
#' under RL or CSR, that is \eqn{E[N_{ii}] = n_i(n_i - 1)/(n - 1)} and \eqn{E[N_{ij}] = n_i n_j/(n - 1)}
#' where \eqn{n_i} is the size of 
#' class \eqn{i} and \eqn{n} is the size of the data set.
#'
#' See also (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table, used in \code{\link{Zcell.nnct.ct}} only 
#' @param varN The variance matrix for cell counts in the NNCT, \code{ct} ; used in \code{\link{Zcell.nnct.ct}} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{\link{Zcell.nnct}} only 
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{\link{Zcell.nnct}} only 
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, for the cell counts, i.e.
#' \eqn{N_{ij}} values
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{\link{Zcell.nnct}} only 
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of Dixon's cell-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the cell counts at the given confidence
#' level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits of the cell counts,
#' provided in \code{conf.level}.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the observed cell counts which is the NNCT}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of 
#' the cell counts.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{\link{Zcell.nnct.ct}} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{\link{Zcell.nnct}} only }
#' 
#' @seealso \code{\link{Zcell.nnct.2s}}, \code{\link{Zcell.nnct.rs}}, \code{\link{Zcell.nnct.ls}},
#' \code{\link{Zcell.nnct.pval}} and \code{\link{Zcell.tct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZcell.nnct
NULL
#'
#' @rdname funsZcell.nnct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' varN
#'
#' Zcell.nnct(Y,cls)
#' Zcell.nnct(Y,cls,alt="g")
#'
#' Zcell.nnct.ct(ct,varN)
#' Zcell.nnct.ct(ct,varN,alt="g")
#'
#' Zcell.nnct(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Zcell.nnct(Y,cls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' Zcell.nnct(Y,cls)
#' Zcell.nnct.ct(ct,varN)
#'
#' @export
Zcell.nnct.ct <- function(ct,varN,alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  if (all(is.na(varN)))
  {stop('All entries of the variance, varN, are NaN (due to data size <=3), so the tests are not defined')}
  
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-ct
  estimate.name <-c("NNCT entries")
  
  EV.<-EV.nnct(ct)
  nullNij<-EV.
  names.null <-"NNCT entries"
  
  k<-nrow(ct);
  
  ts<-stderrN<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    { stderrN[i,j]<-sqrt(varN[i,j])
    ts[i,j] <- (ct[i,j]-EV.[i,j])/stderrN[i,j]
    }
  
  if (all(is.na(ts)))
  {stop('All of the test stat statistics are NaN, these cell-specific tests are not defined')}
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderrN
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderrN
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderrN
           ucl <-estimate+crit.val*stderrN
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("Alternative must be one of less, greater, or two.sided")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  rownames(nullNij)<-colnames(nullNij)<-clnames
  if (!is.null(lcl)) {rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) {rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"Dixon's NNCT cell-specific tests of segregation, Z"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name, #this is for other functions to have a different description for the sample estimates
    null.value = nullNij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZcell.nnct
#'
#' @export
Zcell.nnct <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  EV<-EV.nnct(ct)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  
  res<- Zcell.nnct.ct(ct,varN,alternative=alternative,conf.level=conf.level)
  
  dname <-deparse(substitute(dat))
  
  res$data.name<-dname
  return(res)
} #end for the function
#'

#################################################################

# funsZcell.nnct.pval
#'
#' @title \eqn{p}-values for Cell-specific Z Test Statistics for NNCT
#'
#' @description
#' Four functions: \code{Zcell.nnct.2s}, \code{Zcell.nnct.rs}, \code{Zcell.nnct.ls} and \code{Zcell.nnct.pval}.
#'
#' These functions yield a contingency table (i.e., a matrix) of the \eqn{p}-values for the cell-specific Z
#' test statistics for the NNCT and take the cell-specific \eqn{Z} test statistics in matrix form as their argument.
#' \code{Zcell.nnct.pval} yields an array of size \eqn{k \times k \times 3} where 1st entry of the array is the matrix of \eqn{p}-values for the
#' two-sided alternative, 2nd entry of the array is the matrix of \eqn{p}-values for the left-sided alternative
#' and 3rd entry of the array is the matrix of \eqn{p}-values for the right-sided alternative.
#' And each of \code{Zcell.nnct.2s}, \code{Zcell.nnct.rs} and \code{Zcell.nnct.ls} yield a \eqn{k \times k} matrix of \eqn{p}-values
#' for the two-sided, right-sided and left-sided alternative, respectively.
#' 
#' The functions \code{Zcell.nnct.2s}, \code{Zcell.nnct.rs} and \code{Zcell.nnct.ls} are equivalent to
#' \code{\link{Zcell.nnct}(...,alt)$p.val} where \code{alt="two-sided"}, \code{"greater"} and \code{"less"}, respectively, with the appropriate
#' arguments for the function \code{\link{Zcell.nnct}} (see the examples below).
#'
#' See also (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010;textual}{nnspat}).
#' 
#' @param zt A \eqn{k \times k} matrix of the cell-specific \eqn{Z} test statistics
#' 
#' @return 
#' \code{Zcell.nnct.pval} returns a \eqn{k \times k \times 3} array whose 1st entry is the matrix of \eqn{p}-values for the
#' two-sided alternative, 2nd entry is the matrix of \eqn{p}-values for the left-sided alternative
#' and 3rd entry is the matrix of \eqn{p}-values for the right-sided alternative
#' \code{Zcell.nnct.2s} returns a \eqn{k \times k} matrix of \eqn{p}-values for the two-sided alternative
#' \code{Zcell.nnct.rs} returns a \eqn{k \times k} matrix of \eqn{p}-values for the right-sided alternative
#' \code{Zcell.nnct.ls} returns a \eqn{k \times k} matrix of \eqn{p}-values for the left-sided alternative
#'  
#' @seealso \code{\link{Zcell.nnct}} and \code{\link{Zcell.nnct.ct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZcell.nnct.pval
NULL
#'
#' @rdname funsZcell.nnct.pval
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' TS<-Zcell.nnct(Y,cls)$statistic
#' TS
#' pv<-Zcell.nnct.pval(TS)
#' pv
#'
#' Zcell.nnct(Y,cls,alt="t")$p.val
#' Zcell.nnct(Y,cls,alt="l")$p.val
#' Zcell.nnct(Y,cls,alt="g")$p.val
#'
#' Zcell.nnct.2s(TS)
#'
#' Zcell.nnct.ls(TS)
#'
#' Zcell.nnct.rs(TS)
#'
#' @export
Zcell.nnct.pval <- function(zt) 
{
  k<-length(zt[1,]);
  
  pval<- array(0,dim=c(k,k,3));
  for (i in 1:k)
    for (j in 1:k)
    { pls<- pnorm(zt[i,j]); prs<-1-pls; p2s<-2*min(pls,prs)
    pval[i,j,] <- c(p2s,pls,prs)
    }
  clnames<-rownames(zt)
  ct.names <- c("p-values for two-sided tests","p-values for left-sided tests","p-values for right-sided tests")
  dimnames(pval)<-list(clnames,clnames,ct.names)
  
  pval
} #end for the function
#'
#' @rdname funsZcell.nnct.pval
#'
#' @export 
Zcell.nnct.2s <- function(zt) 
{
  k<-nrow(zt);
  
  pval<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    { pls<-pnorm(zt[i,j])
    pval[i,j] <- 2*min(pls,1-pls)
    }
  
  clnames<-rownames(zt)
  rownames(pval)<-colnames(pval)<-clnames #row and column names from the NNCT
  pval
} #end for the function
#'
#' @rdname funsZcell.nnct.pval
#'
#' @export
Zcell.nnct.ls <- function(zt) 
{
  k<-nrow(zt);
  
  pval<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    {
      pval[i,j] <- pnorm(zt[i,j])
    }
  clnames<-rownames(zt)
  rownames(pval)<-colnames(pval)<-clnames #row and column names from the NNCT
  pval
} #end for the function
#'
#' @rdname funsZcell.nnct.pval
#'
#' @export
Zcell.nnct.rs <- function(zt) 
{
  k<-nrow(zt);
  
  pval<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    {
      pval[i,j] <- 1-pnorm(zt[i,j])
    }
  clnames<-rownames(zt)
  rownames(pval)<-colnames(pval)<-clnames #row and column names from the NNCT
  pval
} #end for the function
#'

#################################################################

#' @title Conversion of a Matrix to a Vector
#'
#' @description
#' Converts the contingency table (or any matrix) \code{ct} to a \code{vector} by default row-wise (i.e., by appending
#' each row one after the other) or column-wise, and also returns the entry indices (in the original matrix \code{ct})
#' in a \eqn{k^2 \times 2} matrix
#' 
#' @param ct A matrix, in particular a contingency table
#' @param byrow A logical argument (default=\code{TRUE}). If \code{TRUE}, rows of \code{ct} are appended to obtain the vector
#' and if \code{FALSE} columns of \code{ct} are appended to obtain the vector.
#'
#' @return A \code{list} with two elements
#' \item{vec}{The \code{vector}ized form the matrix \code{ct}, by default appending the rows of \code{ct}}
#' \item{ind}{The \eqn{k^2 \times 2} matrix of entry indices (in the original matrix \code{ct}) whose i-th row corresponds
#' to the i-th entry in \code{vec}.}
#'
#' @seealso \code{\link{ind.nnsym}} and \code{\link{ind.seg.coeff}},
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#' mat2vec(ct)
#' mat2vec(ct,byrow=FALSE)
#'
#' #an arbitrary 3x3 matrix
#' M<-matrix(sample(10:20,9),ncol=3)
#' M
#' mat2vec(M)
#' mat2vec(M,byrow=FALSE)
#'
#' @export
mat2vec <- function(ct,byrow=TRUE)
{
  k<-nrow(ct)
  
  ifelse(byrow==TRUE,vec<-as.vector(t(ct)),vec<-as.vector(ct))
  
  if (k==1)
  {lt<-rbind(c(1,1))
  return(list(vec=vec,
              ind=lt))}
  
  if (byrow==TRUE)
  {lt<-cbind(rep(1,k),c(1:k))
  for (j in 2:k)
    lt<-rbind(lt,cbind(rep(j,k),c(1:k)));
  } else
  {
    lt<-cbind(c(1:k),rep(1,k))
    for (j in 2:k)
      lt<-rbind(lt,cbind(c(1:k),rep(j,k))); 
  }
  list(vec=vec,
       ind=lt)
} #end for the function
#'

#################################################################

#' @title Covariance Matrix of the Cell Counts in an NNCT
#'
#' @description Returns the covariance matrix of cell counts \eqn{N_{ij}} for \eqn{i,j=1,\ldots,k} in the NNCT, \code{ct}. 
#' The covariance matrix is of dimension \eqn{k^2 \times k^2} and its entries are \eqn{cov(N_{ij},N_{kl})} when \eqn{N_{ij}} values are
#' by default corresponding to the row-wise vectorization of \code{ct}. If \code{byrow=FALSE}, the column-wise 
#' vectorization of \code{ct} is used.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#'
#' See also (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @param ct A nearest neighbor contingency table
#' @param varN The \eqn{k \times k} variance matrix of cell counts of NNCT, \code{ct}.
#' @param Q The number of shared NNs
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs)
#' @param byrow A logical argument (default=\code{TRUE}). If \code{TRUE}, rows of \code{ct} are appended to obtain the vector
#' and if \code{FALSE} columns of \code{ct} are appended to obtain the vector.
#'
#' @return The \eqn{k^2 \times k^2} covariance matrix of cell counts \eqn{N_{ij}} for \eqn{i,j=1,\ldots,k} in the NNCT, \code{ct} 
#'
#' @seealso \code{\link{covNrow2col}}, \code{\link{cov.tct}} and \code{\link{cov.nnsym}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' cov.nnct(ct,varN,Qv,Rv)
#' cov.nnct(ct,varN,Qv,Rv,byrow=FALSE)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' cov.nnct(ct,varN,Qv,Rv)
#' cov.nnct(ct,varN,Qv,Rv,byrow=FALSE)
#' 
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' cov.nnct(ct,varN,Qv,Rv)
#'
#' @export
cov.nnct <- function(ct,varN,Q,R,byrow=TRUE)
{
  varvec<-mat2vec(varN)
  Vsq<-varvec$vec
  ind<-varvec$ind
  
  rs <- row.sum(ct); 
  n<-sum(ct)
  m<-nrow(ct); 
  cov<-matrix(0,m^2,m^2);
  
  for (i in 1:(m^2))
    for (j in i:(m^2))
      if (i==j)
        cov[i,j] <- Vsq[i]
  else 
  {if (ind[i,1]==ind[i,2] & ind[j,1]==ind[j,2])
  {
    k <- ind[i,1]; l<-ind[j,1]
    Paa<- p11(rs[k],n); Pbb<- p11(rs[l],n); Paabb<- p1122(rs[k],rs[l],n)
    cov[i,j] <- (n^2-3*n-Q+R)*Paabb-n^2*Paa*Pbb
    cov[j,i] <- cov[i,j]
  }
    else
    {if (ind[i,1]==ind[i,2] & ind[i,1]==ind[j,1])
    {
      k <- ind[i,1]; l<-ind[j,2]
      Paa<- p11(rs[k],n); Pbb<- p11(rs[l],n); Paaab<- p1112(rs[k],rs[l],n);
      Pab<- p12(rs[k],rs[l],n); Paab<- p112(rs[k],rs[l],n);
      cov[i,j] <- (n-R)*Paab+(n^2-3*n-Q+R)*Paaab-n^2*Paa*Pab
      cov[j,i] <- cov[i,j]
    }
      else
      {if (ind[j,1]==ind[j,2] & ind[j,1]==ind[i,1])
      {
        k <- ind[j,1]; l<-ind[i,2]
        Paa<- p11(rs[k],n); Pbb<- p11(rs[l],n); Paaab<- p1112(rs[k],rs[l],n);
        Pab<- p12(rs[k],rs[l],n); Paab<- p112(rs[k],rs[l],n);
        cov[i,j] <- (n-R)*Paab+(n^2-3*n-Q+R)*Paaab-n^2*Paa*Pab
        cov[j,i] <- cov[i,j]
      }
        else
        {if (ind[i,1]==ind[i,2] & ind[i,1]==ind[j,2])
        {
          k <- ind[i,1]; l<-ind[j,1]
          Paa<- p11(rs[k],n); Paaab<- p1112(rs[k],rs[l],n);
          Pab<- p12(rs[k],rs[l],n); Paab<- p112(rs[k],rs[l],n);
          cov[i,j] <- (n-R+Q)*Paab+(n^2-3*n-Q+R)*Paaab-n^2*Paa*Pab
          cov[j,i] <- cov[i,j]
        }
          else
          {if (ind[j,1]==ind[j,2] & ind[j,1]==ind[i,2])
          {
            k <- ind[j,1]; l<-ind[i,1]
            Paa<- p11(rs[k],n); Paaab<- p1112(rs[k],rs[l],n);
            Pab<- p12(rs[k],rs[l],n); Paab<- p112(rs[k],rs[l],n);
            cov[i,j] <- (n-R+Q)*Paab+(n^2-3*n-Q+R)*Paaab-n^2*Paa*Pab
            cov[j,i] <- cov[i,j]
          }
            else
            {if (ind[i,1]==ind[i,2])
            {
              k <- ind[i,1]; l<-ind[j,1]; q<-ind[j,2] 
              Paa<- p11(rs[k],n); Paabc<- p1123(rs[k],rs[l],rs[q],n); Pbc<- p12(rs[l],rs[q],n); 
              cov[i,j] <- (n^2-3*n-Q+R)*Paabc-n^2*Paa*Pbc
              cov[j,i] <- cov[i,j]
            }
              else
              {if (ind[j,1]==ind[j,2])
              {
                k <- ind[j,1]; l<-ind[i,1]; q<-ind[i,2] 
                Paa<- p11(rs[k],n); Paabc<- p1123(rs[k],rs[l],rs[q],n); Pbc<- p12(rs[l],rs[q],n); 
                cov[i,j] <- (n^2-3*n-Q+R)*Paabc-n^2*Paa*Pbc
                cov[j,i] <- cov[i,j]
              }
                else
                {if (ind[i,1]==ind[j,1])
                {
                  k <- ind[i,1]; l<-ind[i,2]; q<-ind[j,2] 
                  Pab<- p12(rs[k],rs[l],n); Paabc<- p1123(rs[k],rs[l],rs[q],n); Pac<- p12(rs[k],rs[q],n); 
                  cov[i,j] <- (n^2-3*n-Q+R)*Paabc-n^2*Pab*Pac
                  cov[j,i] <- cov[i,j]
                }
                  else
                  {if (ind[i,1]==ind[j,2] & ind[i,2]==ind[j,1] )
                  {
                    k <- ind[i,1]; l<-ind[j,1]; 
                    Pba<-Pab<- p12(rs[k],rs[l],n); Paabb<- p1122(rs[k],rs[l],n); 
                    Paab<- p112(rs[k],rs[l],n); Pabb<- p122(rs[k],rs[l],n)
                    cov[i,j] <- R*Pab+(n-R)*(Paab+Pabb)+(n^2-3*n-Q+R)*Paabb-n^2*Pab*Pba
                    cov[j,i] <- cov[i,j]
                  }
                    else
                    {if (ind[i,2]==ind[j,1])
                    {
                      k <- ind[i,1]; l<-ind[i,2]; q<-ind[j,2]
                      Pab<- p12(rs[k],rs[l],n); Pbc<- p12(rs[l],rs[q],n); Pabc<- p123(rs[k],rs[l],rs[q],n); 
                      Pabbc<- p1223(rs[k],rs[l],rs[q],n)
                      cov[i,j] <- (n-R)*Pabc+(n^2-3*n-Q+R)*Pabbc-n^2*Pab*Pbc
                      cov[j,i] <- cov[i,j]
                    }
                      else
                      {if (ind[j,2]==ind[i,1])
                      {
                        k <- ind[j,1]; l<-ind[j,2]; q<-ind[i,2]
                        Pab<- p12(rs[k],rs[l],n); Pbc<- p12(rs[l],rs[q],n); Pabc<- p123(rs[k],rs[l],rs[q],n); 
                        Pabbc<- p1223(rs[k],rs[l],rs[q],n)
                        cov[i,j] <- (n-R)*Pabc+(n^2-3*n-Q+R)*Pabbc-n^2*Pab*Pbc
                        cov[j,i] <- cov[i,j]
                      }
                        else
                        {if (ind[i,2]==ind[j,2])
                        {
                          k <- ind[i,1]; l<-ind[i,2]; q<-ind[j,1]
                          Pab<- p12(rs[k],rs[l],n); Pbc<- p12(rs[l],rs[q],n); Pabc<- p123(rs[k],rs[l],rs[q],n); 
                          Pabbc<- p1223(rs[k],rs[l],rs[q],n)
                          cov[i,j] <- Q*Pabc+(n^2-3*n-Q+R)*Pabbc-n^2*Pab*Pbc
                          cov[j,i] <- cov[i,j]
                        }
                          else
                          {
                            k <- ind[i,1]; l<-ind[i,2]; q<-ind[j,1]; p<-ind[j,2]
                            Pab<- p12(rs[k],rs[l],n); Pcd<- p12(rs[q],rs[p],n); 
                            Pabcd<- p1234(rs[k],rs[l],rs[q],rs[p],n)
                            cov[i,j] <- (n^2-3*n-Q+R)*Pabcd-n^2*Pab*Pcd
                            cov[j,i] <- cov[i,j]
                          }     
                        }}}}}}}}}}}}
  
  if (byrow==TRUE)
  {
    rcov<-cov #rcov is the result for covariance
  } else # if byrow==FALSE
  { 
    #the part for the covariance of the colum-wise vectorization of N vector
    q<-nrow(ct)
    ind<-vector()
    for (j in 1:q)
    {
      ind<- c(ind,j+seq(0,q-1)*q)
    }
    rcov<-cov[ind,ind] 
  }
  rcov
} #end for the function
#'

#################################################################

#' @title Conversion of the Covariance Matrix of the Row-wise Vectorized Cell Counts to Column-wise
#' Vectorized Cell Counts in an NNCT
#'
#' @description Converts the \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts \eqn{N_{ij}} for \eqn{i,j=1,\ldots,k} 
#' in the NNCT, \code{ct} to the covariance matrix of column-wise vectorized cell counts.
#' In the output, the covariance matrix entries are \eqn{cov(N_{ij},N_{kl})} when \eqn{N_{ij}} values are
#' corresponding to the column-wise vectorization of \code{ct}.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#'
#' See also 
#' (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.
#'
#' @return The \eqn{k^2 \times k^2} covariance matrix of column-wise vectorized cell counts \eqn{N_{ij}} for 
#' \eqn{i,j=1,\ldots,k} in the NNCT, \code{ct}.
#'
#' @seealso \code{\link{cov.nnct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' covNrow<-cov.nnct(ct,varN,Qv,Rv)
#' covNcol1<-cov.nnct(ct,varN,Qv,Rv,byrow=FALSE)
#' covNcol2<-covNrow2col(covNrow)
#'
#' covNrow
#' covNcol1
#' covNcol2
#'
#' all.equal(covNcol1,covNcol2)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' covNrow<-cov.nnct(ct,varN,Qv,Rv)
#' covNcol1<-cov.nnct(ct,varN,Qv,Rv,byrow=FALSE)
#' covNcol2<-covNrow2col(covNrow)
#'
#' covNrow
#' covNcol1
#' covNcol2
#'
#' all.equal(covNcol1,covNcol2)
#'
#' #1D data points
#' n<-20  #or try sample(1:20,1)
#' X<-as.matrix(runif(n))# need to be entered as a matrix with one column
#' #(i.e., a column vector), hence X<-runif(n) would not work
#' ipd<-ipd.mat(X)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' cov.nnct(ct,varN,Qv,Rv)
#'
#' @export
covNrow2col <- function(covN)
{
  k<-sqrt(nrow(covN))
  
  ind<-vector()
  for (j in 1:k)
  {
    ind<- c(ind,j+seq(0,k-1)*k)
  }
  rcov<-covN[ind,ind] 
  rcov
} #end for the function
#'

#################################################################

# funs.vartct
#'
#' @title Functions for Variances of Cell Counts in the Types I, III and IV TCTs
#'
#' @description
#' Three functions: \code{var.tctI}, \code{var.tctIII} and \code{var.tctIV}.
#'
#' These functions return the variances of \eqn{T_{ij}} values for \eqn{i,j=1,\ldots,k} in the TCT in matrix form which
#' is of the same dimension as TCT for types I, III and IV tests. 
#' The argument \code{covN} must be the covariance between \eqn{N_{ij}} values which are obtained from the NNCT by row-wise
#' vectorization.
#' These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.
#' 
#' @return Each of these functions returns a \code{matrix} of same dimension as, \code{ct}, whose entries are the variances of
#' the entries in the TCT for the corresponding type of cell-specific test.
#' The row and column names are inherited from \code{ct}.
#' 
#' @seealso \code{\link{var.tct}} and \code{\link{var.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.vartct
NULL
#'
#' @rdname funs.vartct
#' 
var.tctI <- function(ct,covN) 
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct) 
  
  var<- matrix(0,k,k);
  for (j in 1:k)
  {
    ind<-j+seq(0,k-1)*k 
    VarCj<-sum(covN[ind,ind])
    for (i in 1:k)
    {
      covNijCj<-sum(covN[(i-1)*k+j,ind])
      var[i,j]<- covN[(i-1)*k+j,(i-1)*k+j]+(rs[i]/n)^2*VarCj-2*(rs[i]/n)*covNijCj
    }
  }
  clnames<-rownames(ct)
  rownames(var)<-colnames(var)<-clnames #row and column names from the NNCT 
  var
} #end for the function
#' 
#' @rdname funs.vartct
#'
var.tctIII <- function(ct,covN) 
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct) 
  
  var<- matrix(0,k,k);
  for (j in 1:k)
  {
    ind<-j+seq(0,k-1)*k 
    VarCj<-sum(covN[ind,ind])
    for (i in 1:k)
    {
      covNijCj<-sum(covN[(i-1)*k+j,ind])
      if (i == j)
      {
        var[i,j]<- covN[(i-1)*k+j,(i-1)*k+j]+(rs[i]-1)^2/(n-1)^2*VarCj-2*(rs[i]-1)/(n-1)*covNijCj
      }  
      else 
      {
        var[i,j]<- covN[(i-1)*k+j,(i-1)*k+j]+rs[i]^2/(n-1)^2*VarCj-2*rs[i]/(n-1)*covNijCj
      }
    }
  }
  clnames<-rownames(ct)
  rownames(var)<-colnames(var)<-clnames #row and column names from the NNCT 
  var
} #end for the function
#'
#' @rdname funs.vartct
#'
var.tctIV <- function(ct,covN) 
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct) 
  
  var<- matrix(0,k,k);
  for (j in 1:k)
  {
    ind<-j+seq(0,k-1)*k 
    VarCj<-sum(covN[ind,ind])
    for (i in 1:k)
    {
      covNijCj<-sum(covN[(i-1)*k+j,ind])
      if (i == j)
      {
        var[i,j]<- (rs[i]/n)^2*((n-1)^2/(rs[i]-1)^2*covN[(i-1)*k+j,(i-1)*k+j]+VarCj-2*(n-1)/(rs[i]-1)*covNijCj)
      }  
      else 
      {
        var[i,j]<- (1/n)^2*((n-1)^2*covN[(i-1)*k+j,(i-1)*k+j]+rs[i]^2*VarCj-2*rs[i]*(n-1)*covNijCj)
      }
    }
  }
  
  clnames<-rownames(ct)
  rownames(var)<-colnames(var)<-clnames #row and column names from the NNCT 
  var
} #end for the function
#'
 
#################################################################

#' @title Variances of Entries in a TCT
#'
#' @description Returns the variances of \eqn{T_{ij}} values for \eqn{i,j=1,\ldots,k} in the TCT in matrix form which
#' is of the same dimension as TCT for types I-IV tests. 
#' The argument \code{covN} must be the covariance between \eqn{N_{ij}} values which are obtained from the NNCT by row-wise
#' vectorization. type determines the type of the test for which variances are to be computed, with default=\code{"III"}.
#' These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{ceyhan:SJScorrected2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.
#' @param type The type of the cell-specific test, default=\code{"III"}. Takes on values \code{"I"}-\code{"IV"} (or 
#' equivalently \code{1-4}, respectively.
#'
#' @return A \code{matrix} of same dimension as, \code{ct}, whose entries are the variances of
#' the entries in the TCT for the corresponding type of cell-specific test.
#' The row and column names are inherited from \code{ct}.
#'
#' @seealso \code{\link{var.nnct}}, \code{\link{var.tctI}}, \code{\link{var.tctIII}}, \code{\link{var.tctIV}} 
#' and \code{\link{cov.tct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' var.tct(ct,covN,"I")
#' var.tct(ct,covN,2)
#' var.tct(ct,covN,"III")
#' var.tct(ct,covN,"IV")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' var.tct(ct,covN,"I")
#' var.tct(ct,covN,2)
#'
#' @export
var.tct <- function(ct,covN,type="III") 
{
  
  var <-  switch(type,
         I = { var <- var.tctI(ct,covN) },
         II = { k<-nrow(ct)
         var <- matrix(diag(covN),nrow=k,byrow = TRUE) },
         III = { var <- var.tctIII(ct,covN) },
         IV = { var <- var.tctIII(ct,covN)  }
  )
  
  if (is.null(var)) stop("Enter numbers 1-4 or I-IV in quotes for type")
  
  clnames<-rownames(ct)
  rownames(var)<-colnames(var)<-clnames #row and column names from the NNCT 
  var
} #end for the function
#'

#################################################################

# funs.auxcovtct
#'
#' @title Auxiliary Functions for Computing Covariances Between Cell Counts in the TCT
#'
#' @description
#' Five functions: \code{cov.2cells}, \code{cov.cell.col}, \code{covNijCk}, \code{cov2cols}
#' and \code{covCiCj}
#'
#' These are auxiliary functions for computing covariances between entries in the TCT for the types I-IV
#' cell-specific tests. The covariances between \eqn{T_{ij}} values for \eqn{i,j=1,\ldots,k} in the TCT require covariances
#' between two cells in the NNCT, between a cell and column sum, and between two column sums in the NNCT.
#' \code{cov.2cells} computes the covariance between two cell counts \eqn{N_{ij}} and \eqn{N_{kl}} in an NNCT,
#' \code{cov.cell.col} and \code{covNijCk} are equivalent and they compute the covariance between cell count \eqn{N_{ij}}
#' and sum of column \eqn{k}, \eqn{C_k},
#' \code{cov2cols} and \code{covCiCj} are equivalent and they compute the covariance between sums of two columns, 
#' \eqn{C_i} and \eqn{C_j}.
#' The index arguments refer to which entry or column sum is intended in the NNCT. 
#' The argument \code{covN} must be the covariance between \eqn{N_{ij}} values which are obtained from NNCT by row-wise vectorization.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' @param i,j,k,l Indices of the cell counts or column sums whose covariance is to be computed. All four are
#' needed for \code{cov.2cells} referring to cells \eqn{(i,j)} and \eqn{(k,l)}; only three indices \eqn{i,j,k} are needed for
#' \code{cov.cell.col} and \code{covNijCk} referring to cell \eqn{(i,j)} and column \eqn{k};
#' only two indices \eqn{i,j} are needed for \code{cov2cols} and \code{covCiCj} referring to columns \eqn{i} and \eqn{j}. 
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.
#' 
#' @return 
#' \code{cov.2cells} returns the covariance between two cell counts \eqn{N_{ij}} and \eqn{N_{kl}} in an NNCT,
#' \code{cov.cell.col} and \code{covNijCk} return the covariance between cell count \eqn{N_{ij}}
#' and sum of column \eqn{k}, \eqn{C_k},
#' \code{cov2cols} and \code{covCiCj} return the covariance between sums of two columns, 
#' \eqn{C_i} and \eqn{C_j}.
#' 
#' @seealso \code{\link{cov.tct}} and \code{\link{cov.nnct}}
#' 
#' @name funs.auxcovtct
NULL
#'
#' @rdname funs.auxcovtct
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.2cells(1,1,1,2,ct,covN)
#'
#' cov.cell.col(2,2,1,ct,covN)
#' covNijCk(2,2,1,ct,covN)
#'
#' cov.2cols(2,1,ct,covN)
#' covCiCj(2,1,ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.2cells(2,3,1,2,ct,covN)
#'
#' cov.cell.col(1,1,2,ct,covN)
#' covNijCk(1,1,2,ct,covN)
#'
#' cov.2cols(3,4,ct,covN)
#' covCiCj(3,4,ct,covN)
#' 
#' @export
cov.2cells <- function(i,j,k,l,ct,covN)  
{
  nr<-nrow(ct);
  covN[(i-1)*nr+j,(k-1)*nr+l]
} #end for the function
#' 
#' @rdname funs.auxcovtct
#'
#' @export
cov.cell.col <- function(i,j,k,ct,covN)
{
  nr<-nrow(ct);
  covNijCk<-0
  for (r in 1:nr)
  {
    covNijCk<-covNijCk+covN[(i-1)*nr+j,(r-1)*nr+k]
  }
  covNijCk
} #end for the function
#'
#' @rdname funs.auxcovtct
#'
#' @export
covNijCk <- function(i,j,k,ct,covN)
{
  nr<-nrow(ct);
  
  if (any(c(i,j,k)>nr))
  {stop('Subscript(s) out of bounds of the ct')}
  
  ind<-cbind(rep(1:nr^2,rep(nr^2,nr^2)),rep(1:nr^2),rep(1:nr,rep(nr^3,nr)),rep(1:nr,rep(nr^2,nr)),
             rep(rep(1:nr,rep(nr,nr)),nr),rep(rep(1:nr,nr),nr))
  
  indlis<-which(ind[,3]==i & ind[,4]==j & ind[,6]==k)
  crs<-ind[indlis,][,1:2]
  sum(covN[crs])
} #end for the function
#'
#' @rdname funs.auxcovtct
#'
#' @export
cov.2cols <- function(i,j,ct,covN)
{
  nr<-nrow(ct);
  covCiCj<-0
  for (s in 1:nr)
  {
    for (r in 1:nr)
    {
      covCiCj<-covCiCj+covN[(s-1)*nr+i,(r-1)*nr+j]
    }
  }
  covCiCj
} #end for the function
#'
#' @rdname funs.auxcovtct
#'
#' @export
covCiCj <- function(i,j,ct,covN)
{
  nr<-nrow(ct);
  
  if (any(c(i,j)>nr))
  {stop('Subscript(s) out of bounds of the ct')}
  
  ind<-cbind(rep(1:nr^2,rep(nr^2,nr^2)),rep(1:nr^2),rep(1:nr,rep(nr^3,nr)),rep(1:nr,rep(nr^2,nr)),
             rep(rep(1:nr,rep(nr,nr)),nr),rep(rep(1:nr,nr),nr))
  
  indlis<-which(ind[,4]==i & ind[,6]==j)
  crs<-ind[indlis,][,1:2]
  sum(covN[crs])
} #end for the function
#'

#################################################################

# funs.covtct
#'
#' @title Functions for Covariances of the Entries of the Types I, III and IV TCTs
#'
#' @description
#' Four functions: \code{cov.tctI}, \code{cov.tctIII}, \code{cov.tct3} and \code{cov.tctIV}.
#'
#' These functions return the covariances between between entries in the TCT for the types I, III, and IV
#' cell-specific tests in matrix form which is of dimension \eqn{k^2 \times k^2}.
#' The covariance matrix entries are \eqn{cov(T_{ij},T_{kl})} when \eqn{T_{ij}} values are by default corresponding to 
#' the row-wise vectorization of TCT. 
#' The argument \code{CovN} must be the covariance between \eqn{N_{ij}} values which are obtained from the NNCT by row-wise
#' vectorization.
#' The functions \code{cov.tctIII} and \code{cov.tct3} are equivalent.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table
#' @param CovN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct}.
#' 
#' @return Each of these functions returns a \eqn{k^2 \times k^2} covariance matrix, whose entries are the covariances of
#' the entries in the TCTs for the corresponding type I-IV cell-specific test.
#' The row and column names are inherited from \code{ct}.
#' 
#' @seealso \code{\link{cov.tct}} and \code{\link{cov.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.covtct
NULL
#' 
#' @author Elvan Ceyhan
#'
#' @examples 
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.tctI(ct,covN)
#' cov.tctIII(ct,covN)
#' cov.tctIV(ct,covN)
#' 
#' @rdname funs.covtct
#' 
#' @export
cov.tctI <- function(ct,CovN) 
{
  rs <- row.sum(ct); 
  q<-nrow(ct);
  n<-sum(ct) 
  
  cov<- matrix(0,q^2,q^2);
  
  for (r in 1:(q^2))
    for (s in 1:(q^2))
    {
      i<-ceiling(r/q)
      j<-r-q*floor(r/q-10^(-q))
      
      k<-ceiling(s/q)
      l<-s-q*floor(s/q-10^(-q))
      
      indl<-l+seq(0,q-1)*q
      indj<-j+seq(0,q-1)*q
      
      covNijCl<-sum(CovN[r,indl])
      covNklCj<-sum(CovN[s,indj])
      covCjCl<-sum(CovN[indj,indl])
      
      cov[r,s]<-CovN[r,s]-(rs[k]/n)*covNijCl-(rs[i]/n)*covNklCj+(rs[i]*rs[k]/(n^2))*covCjCl
    }
  cov
} #end for the function
#' 
#' @rdname funs.covtct
#'
#' @export
cov.tctIII <- function(ct,CovN) 
{
  rs <- row.sum(ct); 
  q<-nrow(ct);
  qsq<-q^2
  n<-sum(ct) 
  
  cov<- matrix(0,qsq,qsq);
  for (r in 1:(qsq))
    for (s in 1:(qsq))
    {
      i<-ceiling(r/q) 
      j<-r-q*floor(r/q-10^(-q)) 
      
      k<-ceiling(s/q) 
      l<-s-q*floor(s/q-10^(-q)) 
      
      indl<-l+seq(0,q-1)*q
      indj<-j+seq(0,q-1)*q
      
      covNijCl<-sum(CovN[r,indl])
      covNklCj<-sum(CovN[s,indj])
      covCjCl<-sum(CovN[indj,indl])
      
      if (i==j & k==l)
      {cov[r,s]<-CovN[r,s]-(rs[k]-1)/(n-1)*covNijCl-(rs[i]-1)/(n-1)*covNklCj+(rs[i]-1)*(rs[k]-1)/(n-1)^2*covCjCl}
      else
      {if (i==j & k!=l)
      {cov[r,s]<-CovN[r,s]-rs[k]/(n-1)*covNijCl-(rs[i]-1)/(n-1)*covNklCj+(rs[i]-1)*rs[k]/(n-1)^2*covCjCl}
        else
        {if (i!=j & k==l)
        {cov[r,s]<-CovN[r,s]-(rs[k]-1)/(n-1)*covNijCl-rs[i]/(n-1)*covNklCj+rs[i]*(rs[k]-1)/(n-1)^2*covCjCl}
          else
          {cov[r,s]<-CovN[r,s]-rs[k]/(n-1)*covNijCl-rs[i]/(n-1)*covNklCj+rs[i]*rs[k]/(n-1)^2*covCjCl}          
        }
      }
    }
  cov
} #end for the function
#'
#' @rdname funs.covtct
#' 
#' @export
cov.tct3 <- function(ct,CovN)
{ 
  q<-nrow(ct)
  qsq<-q^2;
  rs <- row.sum(ct);
  n <- sum(ct)
  
  CM<-matrix(0,qsq,qsq)
  for (t in 1:qsq)
  {
    for (u in 1:qsq)
    {
      j<- ifelse(t %% q==0,q,t %% q) #%% is the modulus or mod operator, like 5%%2 is 1
      i<- (t-j)/q+1
      l<- ifelse(u %% q==0,q,u %% q)
      k<- (u-l)/q+1
      if (i==j & k==l)
      {
        CM[t,u]<-cov.2cells(i,j,k,l,ct,CovN)-(rs[i]-1)/(n-1)*cov.cell.col(k,l,i,ct,CovN)-
          (rs[k]-1)/(n-1)*cov.cell.col(i,j,k,ct,CovN)+(rs[i]-1)*(rs[k]-1)/(n-1)^2*cov.2cols(i,k,ct,CovN)  
      }
      
      if (i==j & k!=l)
      {
        CM[t,u]<-cov.2cells(i,j,k,l,ct,CovN)-(rs[i]-1)/(n-1)*cov.cell.col(k,l,i,ct,CovN)-
          rs[k]/(n-1)*cov.cell.col(i,j,l,ct,CovN)+(rs[i]-1)*rs[k]/(n-1)^2*cov.2cols(i,l,ct,CovN)  
      }
      
      if (i!=j & k==l)
      {
        CM[t,u]<-cov.2cells(i,j,k,l,ct,CovN)-(rs[k]-1)/(n-1)*cov.cell.col(i,j,k,ct,CovN)-
          rs[i]/(n-1)*cov.cell.col(k,l,j,ct,CovN)+rs[i]*(rs[k]-1)/(n-1)^2*cov.2cols(j,k,ct,CovN)  
      }
      
      if (i!=j & k!=l)
      {
        CM[t,u]<-cov.2cells(i,j,k,l,ct,CovN)-rs[k]/(n-1)*cov.cell.col(i,j,l,ct,CovN)-
          rs[i]/(n-1)*cov.cell.col(k,l,j,ct,CovN)+rs[i]*rs[k]/(n-1)^2*cov.2cols(j,l,ct,CovN)  
      }
    }
  }
  CM
} #end for the function
#'
#' @rdname funs.covtct
#' 
#' @export
cov.tctIV <- function(ct,CovN) 
{
  rs <- row.sum(ct); 
  q<-nrow(ct);
  n<-sum(ct) 
  
  cov<- matrix(0,q^2,q^2);
  
  for (r in 1:(q^2))
    for (s in 1:(q^2))
    {
      i<-ceiling(r/q)
      j<-r-q*floor(r/q-10^(-q))
      
      k<-ceiling(s/q)
      l<-s-q*floor(s/q-10^(-q))
      
      indl<-l+seq(0,q-1)*q
      indj<-j+seq(0,q-1)*q
      
      covNijCl<-sum(CovN[r,indl])
      covNklCj<-sum(CovN[s,indj])
      covCjCl<-sum(CovN[indj,indl])
      
      if (i==j & k==l)
      {cov[r,s]<-rs[i]*rs[k]/(n^2)*((n-1)^2/((rs[i]-1)*(rs[k]-1))*CovN[r,s]-(n-1)/(rs[i]-1)*covNijCl-(n-1)/(rs[k]-1)*covNklCj+covCjCl)}
      else
      {if (i==j & k!=l)
      {cov[r,s]<-rs[i]/(n^2)*((n-1)^2/(rs[i]-1)*CovN[r,s]-(n-1)*rs[k]/(rs[i]-1)*covNijCl-(n-1)*covNklCj+rs[k]*covCjCl)}
        else
        {if (i!=j & k==l)
        {cov[r,s]<-rs[k]/(n^2)*((n-1)^2/(rs[k]-1)*CovN[r,s]-(n-1)*covNijCl-(n-1)*rs[i]/(rs[k]-1)*covNklCj+rs[i]*covCjCl)}
          else
          {cov[r,s]<-1/(n^2)*((n-1)^2*CovN[r,s]-(n-1)*rs[k]*covNijCl-(n-1)*rs[i]*covNklCj+rs[i]*rs[k]*covCjCl)}          
        }
      }
    }
  cov
} #end for the function
#'

#################################################################

#' @title Covariance Matrix of the Entries of the Type I-IV TCTs
#'
#' @description Returns the covariance matrix of the entries \eqn{T_{ij}} for \eqn{i,j=1,\ldots,k} in the TCT for the types I, III, 
#' and IV cell-specific tests. The covariance matrix is of dimension \eqn{k^2 \times k^2} and its entries are \eqn{cov(T_{ij},T_{kl})}
#' when \eqn{T_{ij}} values are by default corresponding to the row-wise vectorization of TCT. 
#' The argument \code{covN} must be the covariance matrix of \eqn{N_{ij}} values which are obtained from the NNCT by row-wise
#' vectorization.
#' The functions \code{cov.tctIII} and \code{cov.tct3} are equivalent.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#'
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}).
#'
#' @inheritParams var.tct
#'
#' @return The \eqn{k^2 \times k^2} covariance matrix of the entries \eqn{T_{ij}} for \eqn{i,j=1,\ldots,k} in the Type I-IV TCTs
#'
#' @seealso \code{\link{cov.nnct}} and \code{\link{cov.nnsym}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.tct(ct,covN,type=1)
#' cov.tct(ct,covN,type="I")
#' cov.tct(ct,covN,type="II")
#' cov.tct(ct,covN,type="III")
#' cov.tct(ct,covN,type="IV")
#' cov.tctI(ct,covN)
#'
#' cov.tct(ct,covN)
#' cov.tctIII(ct,covN)
#' cov.tct3(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#'
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.tct(ct,covN,type=3)
#' cov.tct(ct,covN,type="III")
#'
#' cov.tctIII(ct,covN)
#' cov.tct3(ct,covN)
#'
#' @export
cov.tct <- function(ct,covN,type="III") 
{
  
  cov <- switch(type,
         I = { cov <- cov.tctI(ct,covN) },
         II = { cov <- covN },
         III = { cov <- cov.tctIII(ct,covN) },
         IV = { cov <- cov.tctIII(ct,covN)  }
  )
  
  if (is.null(cov)) stop("Enter numbers 1-4 or I-IV in quotes for type")  
  
  cov
} #end for the function
#' 

################################################################# 

# funsZcell.tct
#'
#' @title Types I-IV Cell-specific Z Tests of Segregation based on NNCTs
#'
#' @description
#' Two functions: \code{Zcell.tct.ct} and \code{Zcell.tct}.
#'
#' All functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' entries of types I-IV TCT, \eqn{T_{ij}}, from their expected values under RL or CSR for each entry.
#' The test for each entry \eqn{i,j} is based on the normal approximation of the corresponding \eqn{T_{ij}} value
#' and are due to \insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e. null value(s)), lower and upper confidence levels and sample estimates (i.e. observed values)
#' for the \eqn{T_{ij}} values and also names of the test statistics, estimates, null values and the method and the data
#' set used.
#' 
#' The null hypothesis for each entry \eqn{i,j} is that the corresponding value \eqn{T_{ij}} is equal to the expected value
#' under RL or CSR, see \insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat} for more detail.
#'
#' See also (\insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}) and references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Zcell.tct.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized cell counts of NNCT, \code{ct} ;
#' used in \code{Zcell.tct.ct} only.
#' @param type The type of the cell-specific test, default=\code{"III"}. Takes on values \code{"I"}-\code{"IV"} (or 
#' equivalently \code{1-4}, respectively.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zcell.tct} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zcell.tct} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, for the \eqn{T_{ij}} values
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Zcell.tct} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of Types I-IV cell-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the \eqn{T_{ij}} values at the given confidence
#' level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits of the entries, provided in \code{conf.level}.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the observed \eqn{T_{ij}} values which is the TCT}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of 
#' \eqn{T_{ij}} values in the TCT.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zcell.tct.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zcell.tct} only}
#' 
#' @seealso \code{\link{Zcell.nnct.ct}} and \code{\link{Zcell.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZcell.tct
NULL
#'
#' @rdname funsZcell.tct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' type<-"I" #try also "II", "III", and "IV"
#' Zcell.tct(Y,cls,type)
#' Zcell.tct(Y,cls,type,alt="g")
#' Zcell.tct(Y,cls,type,method="max")
#'
#' Zcell.tct.ct(ct,covN)
#' Zcell.tct.ct(ct,covN,type)
#' Zcell.tct.ct(ct,covN,type,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Zcell.tct(Y,cls,type)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' Zcell.tct(Y,cls,type)
#' Zcell.tct.ct(ct,covN,type)
#'
#' @export
Zcell.tct.ct <- function(ct,covN,type="III",alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  Tct<-tct(ct,type=type) 
  estimate<-Tct
  
  ET<-EV.tct(ct,type)
  varT<-var.tct(ct,covN,type)
  
  if (all(is.na(varT)))
  {stop('All entries of the variance, varT, are NaN (due to data size <= 3), so the tests are not defined')}
  
  nullTij<-ET
  estimate.name <- names.null <-paste("TCT ",type," entries",sep="")
  
  k<-nrow(ct);
  
  ts<-stderrT<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    { stderrT[i,j]<-sqrt(varT[i,j])
    ts[i,j] <- (Tct[i,j]-ET[i,j])/stderrT[i,j]
    }
  
  if (all(is.na(ts)))
  {stop('All of the test stat statistics are NaN, the test are not well defined')}
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderrT
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderrT
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderrT
           ucl <-estimate+crit.val*stderrT
         }
  )
 
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-paste("Type ",type," Contingency Table (TCT ",type,") Cell-Specific Tests",sep="")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  rownames(nullTij)<-colnames(nullTij)<-clnames
  if (!is.null(lcl)) {rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) {rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-paste("Type ",type," cell-specific tests of segregation, Z",sep="")
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name, #this is for other functions to have a different description for the sample estimates
    null.value = nullTij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZcell.tct
#'
#' @export
Zcell.tct <- function(dat,lab,type="III",alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv)
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  res<- Zcell.tct.ct(ct,covN,type,alternative=alternative,conf.level=conf.level)
  
  dname <-deparse(substitute(dat))
  res$data.name<-dname
  
  return(res)
} #end for the function
#'

#################################################################

# funsZcell.spec
#'
#' @title Cell-specific Z Tests of Segregation for NNCTs
#'
#' @description
#' Two functions: \code{Zcell.spec.ct} and \code{Zcell.spec}.
#'
#' All functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' entries of NNCT or types I-IV TCTs from the expected values under RL or CSR for each entry.
#' The test for each entry \eqn{i,j} is based on the normal approximation of the corresponding \eqn{T_{ij}} value
#' and are due to \insertCite{dixon:NNCTEco2002;textual}{nnspat}
#' and \insertCite{ceyhan:jkss-posthoc-2017;textual}{nnspat}, respectively.
#' 
#' The \code{type="dixon"} or \code{"nnct"} refers to Dixon's cell-specific test of segregation, and
#' \code{type="I"}-\code{"IV"} refers to types I-IV cell-specific tests, respectively.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e. null value(s)), lower and upper confidence levels and sample estimates (i.e. observed values)
#' for the \eqn{N_{ij}} or \eqn{T_{ij}} values and also names of the test statistics, estimates, null values and the method and
#' the data set used.
#' 
#' The null hypothesis for each entry \eqn{i,j} is that the corresponding value \eqn{N_{ij}} or \eqn{T_{ij}} is equal to the 
#' expected value under RL or CSR.
#' 
#' See also
#' (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Zcell.spec.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{Zcell.spec.ct} only.
#' @param type The type of the cell-specific test with no default.
#' Takes on values \code{"dixon"} or \code{"nnct"} for Dixon's cell-specific tests and \code{"I"}-\code{"IV"} for types I-IV cell-specific
#' tests (or equivalently \code{1-6}, respectively).
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zcell.spec} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zcell.spec} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, for the \eqn{N_{ij}} or \eqn{T_{ij}} values
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Zcell.spec} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of cell-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the \eqn{N_{ij}} or \eqn{T_{ij}} values at the given confidence
#' level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits of the entries, provided in \code{conf.level}.}
#' \item{estimate}{Estimates of the parameters, NNCT or TCT, i.e., matrix of the observed \eqn{N_{ij}} or \eqn{T_{ij}} values
#' which is NNCT or TCT, respectively.}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of the 
#' the null \eqn{N_{ij}} values in an NNCT or \eqn{T_{ij}} values in an TCT.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zcell.spec.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zcell.spec} only}
#' 
#' @seealso \code{\link{Zcell.nnct.ct}}, \code{\link{Zcell.nnct}}, \code{\link{Zcell.tct.ct}} 
#' and \code{\link{Zcell.tct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZcell.spec
NULL
#'
#' @rdname funsZcell.spec
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' type<-"IV" #"dixon" #try also "nnct", "I", "II", "III", and "IV"
#' cell.spec(Y,cls,type)
#' cell.spec(Y,cls,type,alt="g")
#'
#' cell.spec.ct(ct,covN,type)
#' cell.spec.ct(ct,covN,type="II",alt="g")
#'
#' cell.spec(Y,cls,type,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' cell.spec(Y,cls,type="I")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cell.spec(Y,cls,type)
#' cell.spec.ct(ct,covN,type)
#'
#' @export
cell.spec.ct <- function(ct,covN,type,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  k<-nrow(ct)
  if ((type %in% c("dixon","nnct","I","II","III","IV"))==FALSE)
  {stop("type is misspecified, should be one of dixon, nnct, I, II, III or IV (in quotes)")}
  
  if (type %in% c("dixon","nnct"))
  {
    varN<-matrix(diag(covN),byrow=T,ncol=k)
    res<- Zcell.nnct.ct(ct,varN,alternative=alternative,conf.level=conf.level)
  } else
  {
    res<- Zcell.tct.ct(ct,covN,type=type,alternative=alternative,conf.level=conf.level) 
  }
  return(res)
} #end for the function
#'
#' @rdname funsZcell.spec
#'
#' @export
cell.spec <- function(dat,lab,type,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  if ((type %in% c("dixon","nnct","I","II","III","IV"))==FALSE)
  {stop("type is misspecified, should be one of dixon, nnct, I, II, III or IV (in quotes)")}
  
  if (type %in% c("dixon","nnct"))
  {
    res<- Zcell.nnct(dat,lab,alternative=alternative,conf.level=conf.level,...)
  } else
  {
    res<- Zcell.tct(dat,lab,type=type,alternative=alternative,conf.level=conf.level,...) 
  }
  
  dname <-deparse(substitute(dat))
  res$data.name<-dname
  
  return(res)
} #end for the function
#'

#################################################################

# funsC_MI_II
#'
#' @title Correction Matrices for the Covariance Matrix of NNCT entries
#'
#' @description
#' Two functions: \code{correct.cf1} and \code{correct.cf1}.
#'
#' Each function yields matrices which are used in obtaining covariance matrices of \eqn{T_{ij}} values for 
#' types I and II tests from the usual Chi-Square test of contingency tables (i.e. Pielou's test) applied
#' on NNCTs.
#' The output matrices are to be term-by-term multiplied with the covariance matrix of 
#' the entries of NNCT. See Sections 3.1 and 3.2 in 
#' (\insertCite{ceyhan:SJScorrected2010;textual}{nnspat})
#' or
#' Sections 3.5.1 and 3.5.2 in 
#' (\insertCite{ECarXivCorrected:2008;textual}{nnspat}) for more details.
#' 
#' @param ct A nearest neighbor contingency table
#'
#' @return 
#' Both functions return a correction matrix which is to be multiplied with the covariance matrix of
#' entries of the NNCT so as to obtain types I and II overall tests from Pielou's test of segregation.
#' See the description above for further detail.
#' 
#' @seealso \code{\link{nnct.cr1}} and \code{\link{nnct.cr2}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsC_MI_II
NULL
#'
#' @rdname funsC_MI_II
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' #correction type 1
#' CM1<-correct.cf1(ct)
#' CovN.cf1<-covN*CM1
#'
#' #correction type 2
#' CM2<-correct.cf2(ct)
#' CovN.cf2<-covN*CM2
#'
#' covN
#' CovN.cf1
#' CovN.cf2
#'
#' @export
correct.cf1 <- function(ct)
{
  rs <- row.sum(ct);
  cs <- col.sum(ct); 
  
  n<-sum(rs)
  m<-nrow(ct); 
  
  ctvec<-mat2vec(ct)
  lst<-ctvec$ind
  
  cm<-matrix(0,m^2,m^2);
  for (i in 1:(m^2))
    for (j in i:(m^2))
    {
      cm[i,j]<-rs[lst[i,1]]*cs[lst[i,2]]*rs[lst[j,1]]*cs[lst[j,2]]
      cm[j,i]<-cm[i,j]
    }
  n/sqrt(cm)
} #end for the function
#'
#' @rdname funsC_MI_II
#'
#' @export
correct.cf2 <- function(ct)
{
  rs <- row.sum(ct);
  
  n<-sum(rs)
  m<-nrow(ct); 
  
  ctvec<-mat2vec(ct)
  lst<-ctvec$ind
  
  cm<-matrix(0,m^2,m^2);
  for (i in 1:(m^2))
    for (j in i:(m^2))
    {
      cm[i,j]<-rs[lst[i,1]]*rs[lst[i,2]]*rs[lst[j,1]]*rs[lst[j,2]]
      cm[j,i]<-cm[i,j]
    }
  n/sqrt(cm)
} #end for the function
#'

#################################################################

# funsN_I_II
#'
#' @title Correction Matrices for the NNCT entries
#'
#' @description
#' Two functions: \code{nnct.cr1} and \code{nnct.cr1}.
#'
#' Each function yields matrices which are used in obtaining the correction term to be added to 
#' the usual Chi-Square test of contingency tables (i.e. Pielou's test) applied
#' on NNCTs to obtain types I and II overall tests.
#' The output contingency tables are to be row-wise vectorized to obtain \eqn{N_I} and \eqn{N_{II}} vectors. 
#' See Sections 3.1 and 3.2 in 
#' (\insertCite{ceyhan:SJScorrected2010;textual}{nnspat})
#' or
#' Sections 3.5.1 and 3.5.2 in 
#' (\insertCite{ECarXivCorrected:2008;textual}{nnspat}) for more details.
#' 
#' @param ct A nearest neighbor contingency table
#'
#' @return 
#' Both functions return a \eqn{k \times k} contingency table which is to be row-wise vectorized to obtain \eqn{N_I} and \eqn{N_{II}}
#' vectors which are used in the correction summands to obtain types I and II overall tests from Pielou's 
#' test of segregation.
#' See the description above for further detail.
#' 
#' @seealso \code{\link{correct.cf1}} and \code{\link{correct.cf2}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsN_I_II
NULL
#'
#' @rdname funsN_I_II
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' #correction type 1
#' ct1<-nnct.cr1(ct)
#'
#' #correction type 2
#' ct2<-nnct.cr2(ct)
#'
#' ct
#' ct1
#' ct2
#'
#' @export
nnct.cr1 <- function(ct)
{
  rs <- row.sum(ct);
  cs <- col.sum(ct); 
  
  n<-sum(rs)
  m<-nrow(ct); 
  til<-matrix(0,m,m);
  
  for (i in 1:m)
    for (j in 1:m)
    {
      til[i,j]<-(ct[i,j]-rs[i]*cs[j]/n)/sqrt(rs[i]*cs[j]/n)
    }
  til
} #end for the function
#'
#' @rdname funsN_I_II
#'
#' @export
nnct.cr2 <- function(ct)
{
  rs <- row.sum(ct);
  
  n<-sum(rs)
  m<-nrow(ct); 
  til<-matrix(0,m,m);
  
  for (i in 1:m)
    for (j in 1:m)
    {
      til[i,j]<-(ct[i,j]-rs[i]*rs[j]/n)/sqrt(rs[i]*rs[j]/n)
    }
  til
} #end for the function
#'

#################################################################

# funs.base.class.spec
#'
#' @title Base Class-specific Chi-square Tests based on NNCTs
#'
#' @description
#' Two functions: \code{base.class.spec.ct} and \code{base.class.spec}.
#'
#' Both functions are objects of class \code{"classhtest"} but with different arguments (see the parameter list below).
#' Each one performs class specific segregation tests due to Dixon for \eqn{k \ge 2} classes. That is,
#' each one performs hypothesis tests of deviations of 
#' entries in each row of NNCT from the expected values under RL or CSR for each row. 
#' Recall that row labels in the NNCT are base class labels.
#' The test for each row \eqn{i} is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{dixon:NNCTEco2002;textual}{nnspat}.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} for each base class \eqn{i}, description of the 
#' alternative with the corresponding null values (i.e. expected values) for the row \eqn{i}, estimates for the entries in row \eqn{i}
#' for \eqn{i=1,\ldots,k}. The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis for each row is that the corresponding \eqn{N_{ij}} entries in row \eqn{i} are equal to their 
#' expected values under RL or CSR.
#' 
#' See also (\insertCite{dixon:NNCTEco2002,ceyhan:stat-neer-class2009;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{base.class.spec.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{base.class.spec.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{base.class.spec} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{base.class.spec} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{base.class.spec} only
#'    
#' @return A \code{list} with the elements
#' \item{type}{Type of the class-specific test, which is \code{"base"} for this function}
#' \item{statistic}{The \code{vector} of base class-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{vector} of \eqn{p}-values for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k-1} for this function.}
#' \item{estimate}{Estimates of the parameters, NNCT, i.e., matrix of the observed \eqn{N_{ij}} values
#' which is the NNCT.}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of 
#' the \eqn{N_{ij}} values in the NNCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{base.class.spec.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{base.class.spec} only}
#' 
#' @seealso \code{\link{NN.class.spec.ct}}, \code{\link{NN.class.spec}}, \code{\link{class.spec.ct}} 
#' and \code{\link{class.spec}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.base.class.spec
NULL
#'
#' @rdname funs.base.class.spec
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' base.class.spec(Y,cls)
#' base.class.spec.ct(ct,covN)
#' base.class.spec(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' base.class.spec(Y,fcls)
#' base.class.spec.ct(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' base.class.spec(Y,cls)
#' base.class.spec.ct(ct,covN)
#'
#' @export
base.class.spec.ct <- function(ct,covN)
{
  typ="base" #type of the class-specific test, it is base here (for the columns) in the NNCT
  k<-nrow(ct);
  nsq<-mat2vec(ct)$vec #row-wise vectorization
  
  EN<-EV.nnct(ct)
  Ensq <- mat2vec(EN)$vec
  
  ts<-vector()
  for (i in 1:k)
  {
    Ns <- nsq[(k*(i-1)+1):(k*i)];
    ENs<-Ensq[(k*(i-1)+1):(k*i)];
    CovN<-covN[(k*(i-1)+1):(k*i),(k*(i-1)+1):(k*i)];
    ts<- c(ts,t(Ns-ENs) %*% ginv(CovN) %*% (Ns-ENs))
  }
  pval<- 1-pchisq(ts,df=k-1)
  
  method <-c("Dixon's Base Class-Specific Tests")
  ts.names <-"Dixon's base class-specific tests of segregation"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    type=typ,
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    df=k-1,
    estimate = ct,
    null.value = EN,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"classhtest"
  return(rval)
} #end for the function
#'
#' @rdname funs.base.class.spec
#'
#' @export
base.class.spec <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q; Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  rval<-base.class.spec.ct(ct,covN)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsNNclass.spec
#'
#' @title NN Class-specific Chi-square Tests based on NNCTs
#'
#' @description
#' Two functions: \code{NN.class.spec.ct} and \code{NN.class.spec}.
#'
#' Both functions are objects of class \code{"classhtest"} but with different arguments (see the parameter list below).
#' Each one performs class specific segregation tests for the columns, i.e., NN categories for \eqn{k \ge 2} classes.
#' That is,
#' each one performs hypothesis tests of deviations of 
#' entries in each column of NNCT from the expected values under RL or CSR for each column. 
#' Recall that column labels in the NNCT are NN class labels.
#' The test for each column \eqn{i} is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{ceyhan:stat-neer-class2009;textual}{nnspat}.
#' 
#' The argument \code{covN} must be covariance of column-wise vectorization of NNCT if the logical argument \code{byrow=FALSE}
#' otherwise the function converts \code{covN} (which is done row-wise) to columnwise version with \code{\link{covNrow2col}}
#' function.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} for each base class \eqn{i}, description of the 
#' alternative with the corresponding null values (i.e. expected values) for the column \eqn{i}, estimates for the entries in column \eqn{i}
#' for \eqn{i=1,\ldots,k}. The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis for each column is that the corresponding \eqn{N_{ij}} entries in column \eqn{i} are equal to their 
#' expected values under RL or CSR.
#' 
#' See also (\insertCite{dixon:NNCTEco2002,ceyhan:stat-neer-class2009;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{NN.class.spec.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of column-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{NN.class.spec.ct} only.
#' @param byrow A logical argument (default=\code{TRUE}). If \code{TRUE}, rows of \code{ct} are appended to obtain the vector
#' of \eqn{N_{ij}} values and \code{covN} is the covariance matrix for this vector and if \code{FALSE} columns of \code{ct} are appended 
#' to obtain the \eqn{N_{ij}} vector and \code{covN} is converted to the row-wise version by covNrow2col function;used in \code{NN.class.spec.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{NN.class.spec} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{NN.class.spec} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{NN.class.spec} only
#'    
#' @return A \code{list} with the elements
#' \item{type}{Type of the class-specific test, which is \code{"NN"} for this function}
#' \item{statistic}{The \code{vector} of NN class-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{vector} of \eqn{p}-values for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k} for this function.}
#' \item{estimate}{Estimates of the parameters, transpose of the NNCT, i.e., transpose of the matrix of the 
#' observed \eqn{N_{ij}} values which is the transpose of NNCT.}
#' \item{null.value}{Transpose of the matrix of hypothesized null values for the parameters which are expected
#' values of the \eqn{N_{ij}} values in the NNCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{NN.class.spec.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{NN.class.spec} only}
#' 
#' @seealso \code{\link{base.class.spec.ct}}, \code{\link{base.class.spec}}, \code{\link{class.spec.ct}} 
#' and \code{\link{class.spec}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsNNclass.spec
NULL
#'
#' @rdname funsNNclass.spec
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covNrow<-cov.nnct(ct,varN,Qv,Rv)
#' covNcol<-covNrow2col(covNrow)
#'
#' NN.class.spec(Y,cls)
#' NN.class.spec(Y,cls,method="max")
#'
#' NN.class.spec.ct(ct,covNrow)
#' NN.class.spec.ct(ct,covNcol,byrow = FALSE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' NN.class.spec(Y,fcls)
#' NN.class.spec.ct(ct,covNrow)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covNrow<-cov.nnct(ct,varN,Qv,Rv)
#' covNcol<-covNrow2col(covNrow)
#'
#' NN.class.spec(Y,cls)
#'
#' NN.class.spec.ct(ct,covNrow)
#' NN.class.spec.ct(ct,covNcol,byrow = FALSE)
#'
#' @export
NN.class.spec.ct <- function(ct,covN,byrow=TRUE)
{
  typ="NN" #type of the class-specific test, it is NN here (for the columns) in the NNCT
  k<-nrow(ct);
  nsq<-mat2vec(ct,byrow=FALSE)$vec #row-wise vectorization
  
  EN<-EV.nnct(ct)
  Ensq <- mat2vec(EN,byrow=FALSE)$vec
  
  ifelse(byrow==TRUE,covN<-covNrow2col(covN),covN<-covN)
  
  if (all(is.na(covN)))
  {stop('All entries of the covariance matrix, covN, are NaN, so the tests are not defined')}
  
  ts<-vector()
  for (i in 1:k)
  {
    Ns <- nsq[(k*(i-1)+1):(k*i)];
    ENs<-Ensq[(k*(i-1)+1):(k*i)];
    CovN<-covN[(k*(i-1)+1):(k*i),(k*(i-1)+1):(k*i)];
    ts<- c(ts,t(Ns-ENs) %*% solve(CovN) %*% (Ns-ENs))
  }
  pval<- 1-pchisq(ts,df=k)
  
  method <-c("NN Class-Specific Tests")
  ts.names <-"NN class-specific tests of segregation"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    type=typ,
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    df=k,
    estimate = t(ct), # to have columns as rows in the print function
    null.value = t(EN), # to have columns as rows in the print function
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"classhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsNNclass.spec
#'
#' @export
NN.class.spec <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q; Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv) #row-wise vectorization
  
  rval<-NN.class.spec.ct(ct,covN,byrow=TRUE)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funs.class.spec
#'
#' @title Class-specific Chi-square Tests based on NNCTs
#'
#' @description
#' Two functions: \code{class.spec.ct} and \code{class.spec}.
#'
#' Both functions are objects of class \code{"classhtest"} but with different arguments (see the parameter list below).
#' Each one performs class specific segregation tests for the rows if \code{type="base"} and 
#' columns if \code{type="NN"} for \eqn{k \ge 2} classes.
#' That is,
#' each one performs hypothesis tests of deviations of 
#' entries in each row (column) of NNCT from the expected values under RL or CSR for each row (column)
#' if \code{type="base"} (\code{"NN"}). 
#' Recall that row labels of the NNCT are base class labels and
#' column labels in the NNCT are NN class labels.
#' The test for each row (column) \eqn{i} is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{dixon:NNCTEco2002;textual}{nnspat} 
#' (\insertCite{ceyhan:stat-neer-class2009;textual}{nnspat}).
#' 
#' The argument \code{covN} must be covariance of row-wise (column-wise) vectorization of NNCT if \code{type="base"}
#' (\code{type="NN"}).
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} for each base class \eqn{i}, description of the 
#' alternative with the corresponding null values (i.e. expected values) for the row (column) \eqn{i}, estimates for the entries in 
#' row (column) \eqn{i} for \eqn{i=1,\ldots,k} if \code{type="base"} (\code{type="NN"}).
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis for each row (column) is that the corresponding \eqn{N_{ij}} entries in row (column) \eqn{i} are 
#' equal to their expected values under RL or CSR.
#' 
#' See also (\insertCite{dixon:NNCTEco2002,ceyhan:stat-neer-class2009;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{class.spec.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{class.spec.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{class.spec} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{class.spec} only
#' @param type The type of the class-specific tests with default=\code{"base"}. 
#' Takes on values\code{"base"} for (Dixon's) base class-specific test
#' and\code{"NN"} for NN class-specific test.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{class.spec} only
#'    
#' @return A \code{list} with the elements
#' \item{type}{Type of the class-specific test, which is \code{"base"} or \code{"NN"} for this function}
#' \item{statistic}{The \code{vector} of class-specific test statistics}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{vector} of \eqn{p}-values for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k-1} for base class-specific test
#' and \eqn{k} for NN class-specific test.}
#' \item{estimate}{Estimates of the parameters, NNCT, i.e., the matrix of the 
#' observed \eqn{N_{ij}} values for base class-specific test and transpose of the NNCT for
#' the NN class-specific test.}
#' \item{null.value}{The \code{matrix} of hypothesized null values for the parameters which are expected values
#' of the \eqn{N_{ij}} values for the base class-specific test and transpose of this
#' matrix for the NN-class specific test.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{class.spec.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{class.spec} only}
#' 
#' @seealso \code{\link{base.class.spec.ct}}, \code{\link{base.class.spec}}, \code{\link{NN.class.spec.ct}} 
#' and \code{\link{NN.class.spec}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.class.spec
NULL
#'
#' @rdname funs.class.spec
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' class.spec(Y,cls)
#' class.spec(Y,cls,type="NN")
#'
#' class.spec.ct(ct,covN)
#' class.spec.ct(ct,covN,type="NN")
#'
#' class.spec(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' class.spec(Y,fcls)
#' class.spec(Y,fcls,type="NN")
#'
#' class.spec.ct(ct,covN)
#' class.spec.ct(ct,covN,type="NN")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' class.spec(Y,cls)
#' class.spec(Y,cls,type="NN")
#'
#' class.spec.ct(ct,covN)
#' class.spec.ct(ct,covN,type="NN")
#'
#' @export
class.spec.ct <- function(ct,covN,type="base")
{
  if ((type %in% c("NN","base"))==FALSE)
  {stop("type is misspecified, should be one of base or NN (in quotes)")}
  
  if (type == "base")
  {
    res<- base.class.spec.ct(ct,covN)
  } else
  {
    res<- NN.class.spec.ct(ct,covN) 
  }
  
  return(res)
} #end for the function
#'
#' @rdname funs.class.spec
#'
#' @export
class.spec <- function(dat,lab,type="base",...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q; Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  res<-class.spec.ct(ct,covN,type)
  
  dname <-deparse(substitute(dat))
  res$data.name<-dname
  
  return(res)
} #end for the function
#'

#################################################################

# funs.overall.nnct
#'
#' @title Dixon's Overall Test of Segregation for NNCT
#'
#' @description
#' Two functions: \code{overall.nnct.ct} and \code{overall.nnct}.
#'
#' Both functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' cell counts from the expected values under RL or CSR for all cells (i.e., entries) combined in the NNCT.
#' That is, each test is Dixon's overall test of segregation based on NNCTs for \eqn{k \ge 2} classes.
#' This overall test is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{dixon:1994,dixon:NNCTEco2002;textual}{nnspat}.
#' Both functions exclude the last column of the NNCT (in fact any column will do and last column
#' is chosen without loss of generality), to avoid ill-conditioning of the covariance matrix (for its inversion
#' in the quadratic form).
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of NNCT entries, sample estimates (i.e. observed values) of the entries in NNCT.
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all \eqn{N_{ij}} entries are equal to their expected values under RL or CSR.
#'
#' See also 
#' (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:jkss-posthoc-2017;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{overall.nnct.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{overall.nnct.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{overall.nnct} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{overall.nnct} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{overall.nnct} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The overall chi-squared statistic}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)} for this function.}
#' \item{estimate}{Estimates of the parameters, NNCT, i.e., matrix of the observed \eqn{N_{ij}} values
#' which is the NNCT.}
#' \item{est.name,est.name2}{Names of the estimates, former is a longer description of the estimates
#' than the latter.}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of the
#' the \eqn{N_{ij}} values in the NNCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{overall.nnct.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{overall.nnct} only}
#' 
#' @seealso \code{\link{overall.seg.ct}}, \code{\link{overall.seg}}, \code{\link{overall.tct.ct}}
#' and \code{\link{overall.tct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.overall.nnct
NULL
#'
#' @rdname funs.overall.nnct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' overall.nnct(Y,cls)
#' overall.nnct.ct(ct,covN)
#'
#' overall.nnct(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' overall.nnct(Y,fcls)
#' overall.nnct.ct(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' overall.nnct(Y,cls)
#' overall.nnct.ct(ct,covN)
#'
#' @export
overall.nnct.ct <- function(ct,covN)
{
  k<-nrow(ct);
  nsq<-mat2vec(ct)$vec #row-wise vectorization
  
  EN<-EV.nnct(ct)
  Ensq <- mat2vec(EN)$vec
  
  ind<-(1:k)*k #removes the last column in the NNCT, to avoid the rank deficiency
  Nsq<-nsq[-ind]
  ENsq<-Ensq[-ind]
  Covar<-covN[-ind,-ind]
  ts<- t(Nsq-ENsq) %*% ginv(Covar,tol=1.490116e-20) %*% (Nsq-ENsq)
  nu<-k*(k-1)
  
  pval<- 1-pchisq(ts,df=nu)
  
  method <-"Dixon's Overall Test of Segregation"
  
  dname <-deparse(substitute(ct))
  
  estimate<-ct
  estimate.name <-"Nearest Neighbor Contingency Table (NNCT) entries"
  estimate.name2 <-"NNCT entries"
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = EN,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funs.overall.nnct
#'
#' @export
overall.nnct <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q; Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  rval<-overall.nnct.ct(ct,covN)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funs.overall.tct
#'
#' @title Types I-IV Overall Tests of Segregation for NNCT
#'
#' @description
#' Two functions: \code{overall.tct.ct} and \code{overall.tct}.
#'
#' All functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' cell counts from the expected values under RL or CSR for all cells (i.e., entries) combined in the TCT.
#' That is, each test is one of Types I-IV overall test of segregation based on TCTs for \eqn{k \ge 2} classes.
#' This overall test is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{ceyhan:SJScorrected2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}.
#' Both functions exclude some row and/or column of the TCT, to avoid ill-conditioning of the covariance matrix
#' of the NNCT (for its inversion in the quadratic form).
#' In particular, type-II removes the last column, and all other types remove the last row and column.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)} for type II test and \eqn{(k-1)^2} 
#' for the other types, description of the 
#' alternative with the corresponding null values (i.e. expected values) of TCT entries, sample estimates (i.e. observed values) of the entries in TCT.
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all Tij entries for the specified type are equal to their expected values
#' under RL or CSR.
#'
#' See also 
#' (\insertCite{ceyhan:SJScorrected2010,ceyhan:jkss-posthoc-2017;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{overall.tct.ct} only
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct} ;
#' used in \code{overall.tct.ct} only.
#' @param type The type of the overall segregation test, default=\code{"III"}.
#' Takes on values \code{"I"}-\code{"IV"} (or equivalently \code{1-4}, respectively.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{overall.tct} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{overall.tct} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{overall.tct} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The overall chi-squared statistic for the specified type}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)} for type=\code{"II"} and \eqn{(k-1)^2} for others.}
#' \item{estimate}{Estimates of the parameters, TCT, i.e., matrix of the observed \eqn{T_{ij}} values
#' which is the TCT.}
#' \item{est.name,est.name2}{Names of the estimates, former is a longer description of the estimates
#' than the latter.}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of the
#' the \eqn{T_{ij}} values in the TCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{overall.tct.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{overall.tct} only}
#' 
#' @seealso \code{\link{overall.seg.ct}}, \code{\link{overall.seg}}, \code{\link{overall.nnct.ct}}
#' and \code{\link{overall.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.overall.tct
NULL
#'
#' @rdname funs.overall.tct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' overall.tct(Y,cls)
#' overall.tct(Y,cls,type="I")
#' overall.tct(Y,cls,type="II")
#' overall.tct(Y,cls,type="III")
#' overall.tct(Y,cls,type="IV")
#'
#' overall.tct(Y,cls,method="max")
#'
#' overall.tct.ct(ct,covN)
#' overall.tct.ct(ct,covN,type="I")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' overall.tct(Y,fcls)
#' overall.tct.ct(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' overall.tct(Y,cls)
#' overall.tct.ct(ct,covN)
#'
#' @export
overall.tct.ct <- function(ct,covN,type="III")
{
  Tct<-tct(ct,type) 
  estimate<-Tct
  estimate.name <-paste("Type ",type," Contingency Table (TCT ",type,") entries",sep="")
  estimate.name2 <-paste("TCT ",type," entries",sep="")
  
  k<-nrow(ct);
  tsq<-mat2vec(Tct)$vec #row-wise vectorization
  
  ET<-EV.tct(ct,type)
  Etsq <- mat2vec(ET)$vec
  
  covT<-cov.tct(ct,covN,type)
  
  if (type == "II")
  {nu<-k*(k-1)
  ind<-(1:k)*k #removes the last column in the ct
  } else
  {
    nu<-(k-1)^2
    ind<-c((1:(k-1))*k,(k^2-k)+(1:k))#removes the last row and column in the ct
    #this makes type IV equivalent to type III
    #ind<-(k^2-k)+(1:k)#removes last row (makes type IV equivalent to Dixon's)
  }
  
  Tsq<-tsq[-ind]
  ETsq<-Etsq[-ind]
  Covar<-covT[-ind,-ind]
  ts<- t(Tsq-ETsq) %*% ginv(Covar,tol=1.490116e-20) %*% (Tsq-ETsq)
  
  pval<- 1-pchisq(ts,df=nu)
  
  method <-paste("Type ",type," Overall Segregation Test",sep="")
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = ET,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funs.overall.tct
#'
#' @export
overall.tct <- function(dat,lab,type="III",...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q; Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  rval<-overall.tct.ct(ct,covN,type)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funs.overall.seg
#'
#' @title Overall Segregation Tests for NNCTs
#'
#' @description
#' Two functions: \code{overall.seg.ct} and \code{overall.seg}.
#'
#' All functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' cell counts from the expected values under RL or CSR for all cells (i.e., entries) combined in the NNCT or TCT.
#' That is, each test is one of Dixon's or Types I-IV overall test of segregation based on NNCTs or TCTs
#' for \eqn{k \ge 2} classes.
#' Each overall test is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{dixon:1994,dixon:NNCTEco2002;textual}{nnspat}
#' and to \insertCite{ceyhan:SJScorrected2010,ceyhan:jkss-posthoc-2017;textual}{nnspat}, respectively.
#' All functions exclude some row and/or column of the TCT, to avoid ill-conditioning of the covariance matrix
#' of the NNCT (for its inversion in the quadratic form), see the relevant functions under See also section below.
#' 
#' The \code{type="dixon"} or \code{"nnct"} refers to Dixon's overall test of segregation, and
#' \code{type="I"}-\code{"IV"} refers to types I-IV overall tests, respectively.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)} for type II and Dixon's test 
#' and \eqn{(k-1)^2} for the other types, description of the 
#' alternative with the corresponding null values (i.e. expected values) of TCT entries, sample estimates (i.e. observed values) of the entries in TCT.
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all \eqn{N_{ij}} or \eqn{T_{ij}} entries for the specified type are equal to their expected values
#' under RL or CSR, respectively.
#' 
#' See also
#' (\insertCite{dixon:1994,dixon:NNCTEco2002,ceyhan:eest-2010,ceyhan:SJScorrected2010,ceyhan:jkss-posthoc-2017;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{overall.seg.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT, \code{ct};
#' used in \code{overall.seg.ct} only.
#' @param type The type of the overall test with no default.
#' Takes on values \code{"dixon"} or \code{"nnct"} for Dixon's overall test and \code{"I"}-\code{"IV"} for types I-IV cell-specific
#' test (or equivalently \code{1-6}, respectively).
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{overall.seg} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{overall.seg} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{overall.seg} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The overall chi-squared statistic for the specified type}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)} for type II and Dixon's tests
#' and \eqn{(k-1)^2} for others.}
#' \item{estimate}{Estimates of the parameters, NNCT for Dixon's test and type I-IV TCT for others.}
#' \item{est.name,est.name2}{Names of the estimates, former is a longer description of the estimates
#' than the latter.}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of the
#' the \eqn{N_{ij}} values in the NNCT or \eqn{T_{ij}} values in the TCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{overall.seg.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{overall.seg} only}
#'  
#' @seealso \code{\link{overall.nnct.ct}}, \code{\link{overall.nnct}}, \code{\link{overall.tct.ct}}
#' and \code{\link{overall.tct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.overall.seg
NULL
#'
#' @rdname funs.overall.seg
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' type<-"dixon" #try also "nnct", I", "II", "III", and "IV"
#' overall.seg(Y,cls,type)
#' overall.seg(Y,cls,type,method="max")
#' overall.seg(Y,cls,type="I")
#'
#' overall.seg.ct(ct,covN,type)
#' overall.seg.ct(ct,covN,type="I")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' overall.seg(Y,fcls,type="I")
#' overall.seg.ct(ct,covN,type)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' overall.seg(Y,cls,type="I")
#' overall.seg.ct(ct,covN,type)
#'
#' @export
overall.seg.ct <- function(ct,covN,type)
{
  if ((type %in% c("dixon","nnct","I","II","III","IV"))==FALSE)
  {stop("type is misspecified, should be one of dixon, nnct, I, II, III or IV (in quotes)")}
  
  if (type %in% c("dixon","nnct"))
  {
    res<- overall.nnct.ct(ct,covN)
  } else
  {
    res<- overall.tct.ct(ct,covN,type)
  }
  
  return(res)
} #end for the function
#'
#' @rdname funs.overall.seg
#'
#' @export
overall.seg <- function(dat,lab,type,...)
{
  if ((type %in% c("dixon","nnct","I","II","III","IV"))==FALSE)
  {stop("type is misspecified, should be one of dixon, nnct, I, II, III or IV (in quotes)")}
  
  if (type %in% c("dixon","nnct"))
  {
    res<- overall.nnct(dat,lab,...)
  } else
  {
    res<- overall.tct(dat,lab,type,...)
  }
  
  return(res)
} #end for the function
#'

################################################ 
#functions for symmetry, reflexivity, and niche specificity
################################################

# funsXsq.nnsym.ss
#'
#' @title Pielou's First Type of NN Symmetry Test with Chi-square Approximation for multiple classes
#' (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{Xsq.nnsym.ss.ct} and \code{Xsq.nnsym.ss}.
#' 
#' Both functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs the hypothesis test of equality of the expected value of the off-diagonal 
#' cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k \ge 2} classes.
#' That is, each performs Pielou's first type of NN symmetry test which is also equivalent to McNemar's
#' test on the NNCT. The test is appropriate (i.e. have the appropriate asymptotic sampling distribution)
#' provided that data is obtained by sparse sampling.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the chi-squared approximation of the corresponding quadratic form
#' and are due to \insertCite{pielou:1961;textual}{nnspat}.
#' 
#' The argument cont.corr is a logical argument (default=\code{TRUE}) for continuity correction to this test.
#' If \code{TRUE} the continuity correction to McNemar's test is implemented, 
#' and if \code{FALSE} such a correction is not implemented. 
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)/2}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of differences of the off-diagonal entries,(which is
#' 0 for this function) and also the sample estimates (i.e. observed values) of absolute differences of th
#' off-diagonal entries of NNCT (in the upper-triangular form).
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{ij})=E(N_{ji})} for all entries for \eqn{i \ne j} (i.e., symmetry in the 
#' mixed NN structure).
#' In the output, the test statistic, \eqn{p}-value and \code{df} are valid only for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Xsq.nnsym.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Xsq.nnsym.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.nnsym.ss} only
#' @param cont.corr A logical argument (default=\code{TRUE}). 
#' If \code{TRUE} the continuity correction to McNemar's test is implemented, 
#' and if \code{FALSE} such a correction is not implemented. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Xsq.nnsym.ss} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for Pielou's first type of NN symmetry test}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)/2} for this function.}
#' \item{estimate}{Estimates, i.e., absolute differences of the off-diagonal entries of 
#' NNCT (in the upper-triangular form).}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null values for the differences between the expected values of the off-diagonal 
#' entries, which is 0 for this function.}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Xsq.nnsym.ss.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.nnsym.ss} only}
#'  
#' @seealso \code{\link{Znnsym2cl.ss.ct}}, \code{\link{Znnsym2cl.ss}}, \code{\link{Znnsym.ss.ct}},
#' \code{\link{Znnsym.ss}}, \code{\link{Xsq.nnsym.dx.ct}}, \code{\link{Xsq.nnsym.dx}}
#' and \code{\link{Qsym.test}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsXsq.nnsym.ss
NULL
#'
#' @rdname funsXsq.nnsym.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' Xsq.nnsym.ss(Y,cls)
#' Xsq.nnsym.ss.ct(ct)
#'
#' Xsq.nnsym.ss(Y,cls,method="max")
#'
#' Xsq.nnsym.ss(Y,cls,cont.corr=FALSE)
#' Xsq.nnsym.ss.ct(ct,cont.corr=FALSE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Xsq.nnsym.ss(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' Xsq.nnsym.ss(Y,cls)
#' Xsq.nnsym.ss.ct(ct)
#' Xsq.nnsym.ss.ct(ct,cont.corr = FALSE)
#'
#' @export
Xsq.nnsym.ss.ct <- function(ct,cont.corr=TRUE) 
{
  k<-nrow(ct)
  
  if (k==1)
  {stop('contingency table must be kxk with k>=2')}
  
  ts<-0;
  for (i in 1:(k-1))
    for (j in (i+1):k)
    { if (ct[i,j]+ct[j,i]!=0)
    {
      if (cont.corr==TRUE)
      {ts <- ts +(abs(ct[i,j]-ct[j,i])-1)^2/(ct[i,j]+ct[j,i])
      } else
      {ts <- ts +(ct[i,j]-ct[j,i])^2/(ct[i,j]+ct[j,i])}
    }
    }
  nu<-k*(k-1)/2
  pval<- 1-pchisq(ts,df=nu)
  
  ifelse(cont.corr==TRUE,
         method <-"Pielou's Test of NN symmetry for Sparse Sampling - Type I (with Continuity Correction)",
         method <-"Pielou's Test of NN symmetry for Sparse Sampling - Type I (without Continuity Correction)")
  
  dname <-deparse(substitute(ct))
  
  diff.mat<-abs(ct[upper.tri(ct)]- t(ct)[upper.tri(ct)])
  estimate<-matrix(0,k,k)
  estimate[upper.tri(estimate, diag=FALSE)]<-diff.mat 
  
  estimate.name <-c("absolute differences of the off-diagonal entries of NNCT")
  estimate.name2 <-c("absolute differences of the off-diagonal entries of NNCT (in the upper-triangular form)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = 0,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsXsq.nnsym.ss
#'
#' @export
Xsq.nnsym.ss <- function(dat,lab,cont.corr=TRUE,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-Xsq.nnsym.ss.ct(ct,cont.corr) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZnnsym2cl.ss
#'
#' @title Pielou's First Type of NN Symmetry Test with Normal Approximation for Two Classes
#' (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{Znnsym2cl.ss.ct} and \code{Znnsym2cl.ss}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected value of the off-diagonal 
#' cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k=2} classes.
#' That is, each performs Pielou's first type of NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' provided that data is obtained by sparse sampling.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the normal approximation of the difference of the off-diagonal entries
#' in the NNCT and are due to \insertCite{pielou:1961;textual}{nnspat}.
#'
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, estimate and null value for the parameter of interest
#' (which is the difference of the off-diagonal entries in the NNCT), and method and name of the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{12})=E(N_{21})} in the \eqn{2 \times 2} NNCT (i.e., symmetry in the 
#' mixed NN structure).
#' In the output, the test statistic, \eqn{p}-value and the confidence interval are valid only 
#' for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Znnsym2cl.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnsym2cl.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnsym2cl.ss} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnsym2cl.ss} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for Pielou's first type of NN symmetry test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}} in the \eqn{2 \times 2} NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate, i.e., the difference of the off-diagonal entries of the \eqn{2 \times 2} NNCT, \eqn{N_{12}-N_{21}}.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{12})-E(N_{21})} in the \eqn{2 \times 2} NNCT, which is 0 for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link{Xsq.nnsym.ss.ct}}, \code{\link{Xsq.nnsym.ss}}, \code{\link{Znnsym.ss.ct}} and
#' \code{\link{Znnsym.ss}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnsym2cl.ss
NULL
#'
#' @rdname funsZnnsym2cl.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' Znnsym2cl.ss(Y,cls)
#' Znnsym2cl.ss.ct(ct)
#'
#' Znnsym2cl.ss(Y,cls,method="max")
#'
#' Znnsym.ss.ct(ct)
#'
#' Znnsym2cl.ss(Y,cls,alt="g")
#' Znnsym2cl.ss.ct(ct,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Znnsym2cl.ss(Y,fcls)
#' 
#' #############
#' ct<-matrix(sample(1:20,4),ncol=2)
#' Znnsym2cl.ss.ct(ct) #gives an error message if ct<-matrix(sample(1:20,9),ncol=3)
#'
#' @export
Znnsym2cl.ss.ct <- function(ct,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  ts<-stderr<- 0
  {stderr <-sqrt(ct[1,2]+ct[2,1])
    if (stderr!=0)
    { ts <-(ct[1,2]-ct[2,1])/stderr}
  }
  names(ts) <-"standardized difference between # of mixed NNs, Z"
  estimate<-ct[1,2]-ct[2,1]
  
  names(estimate) <-c("difference between # of mixed NNs")
  method <-c("Pielou's NN Symmetry Test for Two Classes (for Sparse Sampling)")
  
  nullij<-0
  names(nullij) <-"(expected) difference between # of mixed NNs"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = nullij,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZnnsym2cl.ss
#'
#' @export
Znnsym2cl.ss <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-Znnsym2cl.ss.ct(ct,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZnnsym2cl.dx
#'
#' @title Dixon's NN Symmetry Test with Normal Approximation for Two Classes
#'
#' @description
#' Two functions: \code{Znnsym2cl.dx.ct} and \code{Znnsym2cl.dx}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs the hypothesis test of equality of the expected value of the off-diagonal 
#' cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k=2} classes.
#' That is, each performs Dixon's NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the normal approximation of the difference of the off-diagonal entries
#' in the NNCT and are due to \insertCite{dixon:1994;textual}{nnspat}.
#'
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, estimate and null value for the parameter of interest
#' (which is the difference of the off-diagonal entries in the NNCT), and method and name of the data set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{12})=E(N_{21})} in the \eqn{2 \times 2} NNCT (i.e., symmetry in the 
#' mixed NN structure).
#' 
#' See also
#' (\insertCite{dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Znnsym2cl.dx.ct} only 
#' @param Q The number of shared NNs, used in \code{Znnsym2cl.dx.ct} only
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs),
#' used in \code{Znnsym2cl.dx.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnsym2cl.dx} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnsym2cl.dx} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnsym2cl.dx} only
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for Pielou's first type of NN symmetry test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}} in the \eqn{2 \times 2} NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate, i.e., the difference of the off-diagonal entries of the \eqn{2 \times 2} NNCT, \eqn{N_{12}-N_{21}}.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{12})-E(N_{21})} in the \eqn{2 \times 2} NNCT, which is 0 for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link{Znnsym2cl.ss.ct}}, \code{\link{Znnsym2cl.ss}}, \code{\link{Znnsym.dx.ct}},
#' \code{\link{Znnsym.dx}}, \code{\link{Xsq.nnsym.dx.ct}} and \code{\link{Xsq.nnsym.dx}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnsym2cl.dx
NULL
#'
#' @rdname funsZnnsym2cl.dx
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' Znnsym2cl.dx(Y,cls)
#' Znnsym2cl.dx.ct(ct,Qv,Rv)
#'
#' Znnsym2cl.dx(Y,cls,method="max")
#'
#' Znnsym2cl.dx(Y,cls,alt="g")
#' Znnsym2cl.dx.ct(ct,Qv,Rv,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Znnsym2cl.dx(Y,fcls)
#'
#' #############
#' ct<-matrix(sample(1:20,4),ncol=2)
#' Znnsym2cl.dx.ct(ct,Qv,Rv) #gives an error message if ct<-matrix(sample(1:20,9),ncol=3)
#' #here, Qv and Rv values are borrowed from above, to highlight a point
#'
#' @export
Znnsym2cl.dx.ct <- function(ct,Q,R,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  rs <- row.sum(ct); 
  n<-sum(rs)
  n1<-rs[1]; n2<-rs[2];
  n12<-ct[1,2]; n21<-ct[2,1]
  
  VarN12<-n*p12(n1,n2,n)+Q*p112(n1,n2,n)+(n^2-3*n-Q+R)*p1122(n1,n2,n)-(n*p12(n1,n2,n))^2
  VarN21<-n*p12(n2,n1,n)+Q*p112(n2,n1,n)+(n^2-3*n-Q+R)*p1122(n2,n1,n)-(n*p12(n2,n1,n))^2
  CovN12N21<-R*p12(n1,n2,n)+(n-R)*(p112(n1,n2,n)+p122(n1,n2,n))+(n^2-3*n-Q+R)*p1122(n1,n2,n)-n^2*p12(n1,n2,n)*p12(n2,n1,n)
  stderr <-sqrt(VarN12+VarN21-2*CovN12N21)
  ts<-(n12-n21)/stderr
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the NN symmetry test is not defined')}
  
  names(ts) <-"standardized difference between # of mixed NNs, Z"
  estimate<-n12-n21
  
  names(estimate) <-c("difference between # of mixed NNs")
  method <-c("Dixon's NN Symmetry Test for Two Classes")
  
  nullij<-0
  names(nullij) <-"(expected) difference between # of mixed NNs"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = nullij,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZnnsym2cl.dx
#'
#' @export
Znnsym2cl.dx <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  
  rval<-Znnsym2cl.dx.ct(ct,Qv,Rv,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

#' @title NN Symmetry Test with Normal Approximation for Two Classes
#'
#' @description
#' An object of class \code{"htest"} performing hypothesis test of equality of the expected value of the off-diagonal 
#' cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k=2} classes.
#' That is, the test performs Dixon's or Pielou's (first type of) NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data and for sparsely sample data, respectively.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' The symmetry test is based on the normal approximation of the difference of the off-diagonal entries
#' in the NNCT and are due to \insertCite{pielou:1961,dixon:1994;textual}{nnspat}.
#' 
#' The \code{type="dixon"} refers to Dixon's NN symmetry test and 
#' \code{type="pielou"} refers to Pielou's first type of NN symmetry test.
#'
#' The function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, estimate and null value for the parameter of interest
#' (which is the difference of the off-diagonal entries in the NNCT), and method and name of the data set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{12})=E(N_{21})} in the \eqn{2 \times 2} NNCT (i.e., symmetry in the 
#' mixed NN structure).
#' 
#' See also
#' (\insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param type The type of the NN symmetry test with default=\code{"dixon"}.
#' Takes on values \code{"dixon"} and \code{"pielou"} for Dixon's and Pielou's (first type) NN symmetry test
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#'   
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for Pielou's first type of NN symmetry test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}} in the \eqn{2 \times 2} NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate, i.e., the difference of the off-diagonal entries of the \eqn{2 \times 2} NNCT, \eqn{N_{12}-N_{21}}.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{12})-E(N_{21})} in the \eqn{2 \times 2} NNCT, which is 0 for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link{Znnsym2cl.ss.ct}}, \code{\link{Znnsym2cl.ss}}, \code{\link{Znnsym2cl.dx.ct}},
#' \code{\link{Znnsym2cl.dx}}, \code{\link{Znnsym.ss.ct}}, \code{\link{Znnsym.ss}}, \code{\link{Znnsym.dx.ct}},
#' \code{\link{Znnsym.dx}}, \code{\link{Znnsym.dx.ct}}, \code{\link{Znnsym.dx}} and \code{\link{Znnsym}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Znnsym2cl(Y,cls)
#' Znnsym2cl(Y,cls,type="pielou")
#'
#' Znnsym2cl(Y,cls,alt="g")
#' Znnsym2cl(Y,cls,type="pielou",alt="g")
#'
#' @export
Znnsym2cl <- function(dat,lab,type="dixon",alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  if ((type %in% c("dixon","pielou"))==FALSE)
  {stop("type is misspecified, should be one of dixon or pielou (in quotes)")}
  
  if (type == "dixon")
  {
    res<- Znnsym2cl.dx(dat,lab,alternative=alternative,conf.level = conf.level) 
  } else
  {
    res<- Znnsym2cl.ss(dat,lab,alternative=alternative,conf.level = conf.level) 
  }
  
  return(res)
} #end for the function
#'

#################################################################

# funsZnnsym.ss
#'
#' @title Pielou's Pairwise NN Symmetry Test with Normal Approximation (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{Znnsym.ss.ct} and \code{Znnsym.ss}.
#' 
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the off-diagonal 
#' cell counts (i.e., entries) for each pair \eqn{i,j} of classes under RL or CSR in the NNCT for \eqn{k \ge 2} classes.
#' That is, each performs Pielou's first type of NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' provided that data is obtained by sparse sampling.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the normal approximation of the differences of the off-diagonal entries
#' in the NNCT and are due to \insertCite{pielou:1961;textual}{nnspat}.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values, lower and upper confidence levels, sample estimates (i.e. observed values)
#' and null value(s) (i.e. expected values) for the
#' \eqn{N_{ij}-N_{ji}} values for \eqn{i \ne j} (all in the upper-triangular form except for the null value, which is 0 for all
#' pairs) and also names of the test statistics, estimates, null values and the method and the data
#' set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT (i.e., symmetry in the 
#' mixed NN structure) for \eqn{k \ge 2}.
#' In the output, the test statistic, \eqn{p}-value and the lower and upper confidence limits are valid only 
#' for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Znnsym.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnsym.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnsym.ss} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{ij}-N_{ji}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnsym.ss} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of \eqn{Z} test statistics for Pielou's first type of NN symmetry test
#' (in the upper-triangular form)}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative
#' (in the upper-triangular form)}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels (in the upper-triangular form) for the \eqn{N_{ij}-N_{ji}}
#' values for \eqn{i \ne j} at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits (i.e., conf.level) of the differences of the
#' off-diagonal entries.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the difference of the off-diagonal entries 
#' (in the upper-triangular form) of the \eqn{k \times k} NNCT, \eqn{N_{ij}-N_{ji}} for \eqn{i \ne j}.}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{ij})-E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT, which is 0 for this function.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Znnsym.ss.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Znnsym.ss} only}
#'  
#' @seealso \code{\link{Znnsym.dx.ct}}, \code{\link{Znnsym.dx}}, \code{\link{Znnsym2cl.ss.ct}} and
#' \code{\link{Znnsym2cl.ss}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnsym.ss
NULL
#'
#' @rdname funsZnnsym.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' Znnsym.ss(Y,cls)
#' Znnsym.ss.ct(ct)
#'
#' Znnsym.ss(Y,cls,method="max")
#'
#' Znnsym.ss(Y,cls,alt="g")
#' Znnsym.ss.ct(ct,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Znnsym.ss(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' Znnsym.ss(Y,cls)
#' Znnsym.ss.ct(ct)
#'
#' @export
Znnsym.ss.ct <- function(ct,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  k<-nrow(ct);
  if (k==1)
  {stop('contingency table must be kxk with k>=2')}
  
  diff.mat<-abs(ct[upper.tri(ct)]- t(ct)[upper.tri(ct)])
  estimate<-matrix(0,k,k)
  estimate[upper.tri(estimate, diag=FALSE)]<-diff.mat 
  
  estimate.name <-c("absolute differences of the off-diagonal entries of NNCT")
  estimate.name2 <-c("absolute differences of the off-diagonal entries of NNCT (in the upper-triangular form)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  nullij<-0
  names.null <-"differences of off-diagonal NNCT entries"
  
  ts<-stderr<- matrix(0,k,k);
  for (i in 1:(k-1))
    for (j in (i+1):k)
    {stderr[i,j] <-sqrt(ct[i,j]+ct[j,i])
    if (stderr[i,j]!=0)
    { ts[i,j] <-ts[i,j] +(ct[i,j]-ct[j,i])/stderr[i,j]}
    }
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderr
           ucl <-estimate+crit.val*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("Pielou's Pairwise NN Symmetry Tests (for Sparse Sampling)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  
  if (!is.null(lcl)) { lcl[!upper.tri(lcl)]<-0;
  rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) { ucl[!upper.tri(ucl)]<-0;
  rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"Pielou's pairwise NN tests of symmetry (in the upper-triangular form)"
  
  dname <-deparse(substitute(ct))
  
  pval[!upper.tri(pval)]<-0
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = nullij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZnnsym.ss
#'
#' @export
Znnsym.ss <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-Znnsym.ss.ct(ct,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

#' @title Index Matrix for Computing the Covariance of Dixon's Overall NN Symmetry Test
#'
#' @description Returns the index matrix for choosing the entries in the covariance matrix for NNCT 
#' used for computing the covariance for Dixon's NN symmetry test.
#' The matrix is \eqn{k(k-1)/2 \times 2} with each row is the \eqn{i,j} 
#' corresponding to \eqn{N_{ij}} in the NNCT.
#'
#' @param k An integer specifying the number of classes in the data set
#' 
#' @return The \eqn{k(k-1)/2 \times 2} index matrix with each row is the \eqn{i,j} 
#' corresponding to \eqn{N_{ij}} in the NNCT
#' 
#' @seealso \code{\link{cov.nnsym}} and \code{\link{ind.seg.coeff}}
#'
#' @author Elvan Ceyhan
#'
ind.nnsym <- function(k)
{
  ind1<-vector()
  if (k==1)
  {stop('k must be > 1')}
  for (i in 1:(k-1))
  {
    for (j in (i+1):k)
    {
      ind1<-c(ind1,c(i,j))
    }
  }
  ind.mat1<-matrix(ind1,ncol=2,byrow=T) #indices for the \code{T} vector
  ind.mat1
} #end for the function
#'

################################################

#' @title Variances of Differences of Off-Diagonal Entries in an NNCT
#'
#' @description Returns the variances of differences of off-diagonal cell counts \eqn{N_{ij}-N_{ji}} for \eqn{i,j=1,\ldots,k} and \eqn{i \ne j}
#' in the NNCT, \code{ct} in a vector of length \eqn{k(k-1)/2}, the order of \eqn{i,j} for \eqn{N_{ij}-N_{ji}}
#' is as in the output of \code{\link{ind.nnsym}(k)}.
#' These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat}).
#'
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#'
#' @return A \code{vector} of length \eqn{k(k-1)/2}, whose entries are the variances of differences of off-diagonal 
#' cell counts \eqn{N_{ij}-N_{ji}} for \eqn{i,j=1,\ldots,k} and \eqn{i \ne j} in the NNCT.
#'
#' @seealso \code{\link{var.nnct}}, \code{\link{var.tct}} and \code{\link{cov.nnct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' var.nnsym(covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' var.nnsym(covN)
#'
#' @export
var.nnsym <- function(covN)
{
  k<-sqrt(nrow(covN))
  ent2<-cbind(rep(1:k,rep(k,k)),rep(1:k,k))
  ind.mat1<-ind.nnsym(k)
  
  kl<-k*(k-1)/2
  varT<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the N vector
    r1<-rc1[1]
    c1<-rc1[2]
    ri1<-which(ent2[,1]==r1 & ent2[,2]==c1) #row indices for the N vector
    ri2<-which(ent2[,1]==c1 & ent2[,2]==r1)
    
    varT<-c(varT, covN[ri1,ri1]-covN[ri1,ri2]-covN[ri2,ri1]+covN[ri2,ri2])
  }
  varT
} #end for the function
#'

#################################################################

# funsZnnsym.dx
#'
#' @title Dixon's Pairwise NN Symmetry Test with Normal Approximation
#'
#' @description
#' Two functions: \code{Znnsym.dx.ct} and \code{Znnsym.dx}.
#' 
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the off-diagonal 
#' cell counts (i.e., entries) for each pair \eqn{i,j} of classes under RL or CSR in the NNCT for \eqn{k \ge 2} classes.
#' That is, each performs Dixon's NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the normal approximation of the difference of the off-diagonal entries
#' in the NNCT and are due to \insertCite{dixon:1994;textual}{nnspat}.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e. null value(s)), lower and upper confidence levels and sample estimates (i.e. observed values)
#' for the \eqn{N_{ij}-N_{ji}} values for \eqn{i \ne j} (all in the upper-triangular form except for the null value, which is 0
#' for all pairs) and also names of the test statistics, estimates, null values and the method and the data
#' set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT (i.e., symmetry in the 
#' mixed NN structure) for \eqn{k \ge 2}.
#' In the output, the test statistic, \eqn{p}-value and the lower and upper confidence limits are valid  
#' for completely mapped data.
#' 
#' See also
#' (\insertCite{dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Znnsym.dx.ct} only 
#' @param varS The variance vector of differences of off-diagonal cell counts in NNCT, \code{ct} , usually output 
#' of var.nnsym function.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnsym.dx} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnsym.dx} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{ij}-N_{ji}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnsym.dx} only
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of \eqn{Z} test statistics for Dixon's NN symmetry test
#' (in the upper-triangular form)}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative
#' (in the upper-triangular form)}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels (in the upper-triangular form) for the \eqn{N_{ij}-N_{ji}}
#' values for \eqn{i \ne j} at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits (i.e., conf.level) of the differences of the
#' off-diagonal entries.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the difference of the off-diagonal entries
#' (in the upper-triangular form) of the \eqn{k \times k} NNCT, \eqn{N_{ij}-N_{ji}} for \eqn{i \ne j}.}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{ij})-E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT, which is 0 for this function.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Znnsym.dx.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Znnsym.dx} only}
#'  
#' @seealso \code{\link{Znnsym2cl.dx.ct}}, \code{\link{Znnsym2cl.dx}}, \code{\link{Znnsym.ss.ct}},
#' \code{\link{Znnsym.ss}}, \code{\link{Xsq.nnsym.dx.ct}} and \code{\link{Xsq.nnsym.dx}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnsym.dx
NULL
#'
#' @rdname funsZnnsym.dx
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' varS<-var.nnsym(covN)
#'
#' Znnsym.dx(Y,cls)
#' Znnsym.dx.ct(ct,varS)
#'
#' Znnsym.dx(Y,cls,method="max")
#'
#' Znnsym.dx(Y,cls,alt="g")
#' Znnsym.dx.ct(ct,varS,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Znnsym.dx(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' varS<-var.nnsym(covN)
#'
#' Znnsym.dx(Y,cls)
#' Znnsym.dx.ct(ct,varS)
#'
#' @export
Znnsym.dx.ct <- function(ct,varS,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  k<-nrow(ct)
  
  diff.mat<-abs(ct[upper.tri(ct)]- t(ct)[upper.tri(ct)])
  estimate<-matrix(0,k,k)
  estimate[upper.tri(estimate, diag=FALSE)]<-diff.mat 
  
  estimate.name <-c("absolute differences of the off-diagonal entries of NNCT")
  estimate.name2 <-c("absolute differences of the off-diagonal entries of NNCT (in the upper-triangular form)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  nullij<-0
  names.null <-"differences of off-diagonal NNCT entries"
  
  ind.mat1<-ind.nnsym(k)
  
  kl<-k*(k-1)/2
  Tvec<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]
    c1<-rc1[2]
    Tvec<-c(Tvec,ct[r1,c1]-ct[c1,r1]) #entry in ct corresponding to index i of T
  }
  ZT<- Tvec/sqrt(varS)
  
  if (all(is.na(ZT)))
  {stop('The test statistics are all NaN, so the NN symmetry tests are not defined')}
  
  ts<-stderr<- matrix(0,k,k);
  stderr[lower.tri(stderr, diag=FALSE)] <- sqrt(varS)
  stderr <- t(stderr)
  ts[lower.tri(ts, diag=FALSE)] <- ZT
  ts <- t(ts)
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderr
           ucl <-estimate+crit.val*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("Dixon's Pairwise NN Symmetry Tests")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  
  if (!is.null(lcl)) { lcl[!upper.tri(lcl)]<-0;
  rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) { ucl[!upper.tri(ucl)]<-0;
  rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"Dixon's pairwise NN tests of symmetry (in the upper-triangular form)"
  
  dname <-deparse(substitute(ct))
  
  pval[!upper.tri(pval)]<-0
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = nullij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZnnsym.dx
#'
#' @export
Znnsym.dx <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
  
  varS<-var.nnsym(covN)
  
  rval<-Znnsym.dx.ct(ct,varS,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

#' @title NN Symmetry Test with Normal Approximation
#'
#' @description
#' An object of class \code{"cellhtest"} performing hypothesis test of equality of the expected values of the
#' off-diagonal cell counts (i.e., entries) for each pair \eqn{i,j} of classes under RL or CSR in the NNCT
#' for \eqn{k \ge 2} classes.
#' That is, the test performs Dixon's or Pielou's (first type of) NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data or for sparsely sample data, respectively.
#' (See \insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' The \code{type="dixon"} refers to Dixon's NN symmetry test and 
#' \code{type="pielou"} refers to Pielou's first type of NN symmetry test.
#' The symmetry test is based on the normal approximation of the difference of the off-diagonal entries
#' in the NNCT and are due to \insertCite{pielou:1961,dixon:1994;textual}{nnspat}.
#' 
#' The function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e. null value(s)), lower and upper confidence levels and sample estimate
#' for the \eqn{N_{ij}-N_{ji}} values for \eqn{i \ne j} (all in the upper-triangular form except for the null value, which is 0
#' for all pairs) and also names of the test statistics, estimates, null values and the method and the data
#' set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT (i.e., symmetry in the 
#' mixed NN structure) for \eqn{k \ge 2}.
#' In the output, if if \code{type="pielou"}, 
#' the test statistic, \eqn{p}-value and the lower and upper confidence limits are valid only 
#' for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param type The type of the NN symmetry test with default=\code{"dixon"}.
#' Takes on values \code{"dixon"} and \code{"pielou"} for Dixon's and Pielou's (first type) NN symmetry test
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function
#'   
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of \eqn{Z} test statistics for the NN symmetry test (in the upper-triangular form)}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative
#' (in the upper-triangular form)}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels (in the upper-triangular form) for the \eqn{N_{ij}-N_{ji}}
#' values for \eqn{i \ne j} at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits (i.e., conf.level) of the differences of the
#' off-diagonal entries.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the difference of the off-diagonal entries
#' (in the upper-triangular form) of the \eqn{k \times k} NNCT, \eqn{N_{ij}-N_{ji}} for \eqn{i \ne j}.}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null value for the expected difference between the off-diagonal entries, 
#' \eqn{E(N_{ij})-E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT, which is 0 for this function.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link{Znnsym.ss.ct}}, \code{\link{Znnsym.ss}}, \code{\link{Znnsym.dx.ct}},
#' \code{\link{Znnsym.dx}} and \code{\link{Znnsym2cl}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Znnsym(Y,cls)
#' Znnsym(Y,cls,method="max")
#' Znnsym(Y,cls,type="pielou")
#' Znnsym(Y,cls,type="pielou",method="max")
#'
#' Znnsym(Y,cls,alt="g")
#' Znnsym(Y,cls,type="pielou",alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Znnsym(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Znnsym(Y,cls)
#' Znnsym(Y,cls,type="pielou")
#'
#' @export
Znnsym <- function(dat,lab,type="dixon",alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  if ((type %in% c("dixon","pielou"))==FALSE)
  {stop("type is misspecified, should be one of dixon or pielou (in quotes)")}
  
  if (type == "dixon")
  {
    res<- Znnsym.dx(dat,lab,alternative=alternative,conf.level = conf.level,...) 
  } else
  {
    res<- Znnsym.ss(dat,lab,alternative=alternative,conf.level = conf.level,...) 
  }
  return(res)
} #end for the function
#'

#################################################################

#' @title Covariance Matrix of the Differences of the Off-Diagonal Cell Counts in an NNCT
#'
#' @description Returns the covariance matrix of the differences of the cell counts, \eqn{N_{ij}-N_{ji}} 
#' for \eqn{i,j=1,\ldots,k} and \eqn{i \ne j}, in the NNCT, \code{ct}.
#' The covariance matrix is of dimension \eqn{k(k-1)/2 \times k(k-1)/2} and its entries are
#' \eqn{cov(N_{ij}-N_{ji}, N_{kl}-N_{lk})} where the order of \eqn{i,j} for \eqn{N_{ij}-N_{ji}} is as
#' in the output of \code{\link{ind.nnsym}(k)}. 
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' The argument \code{covN} is the covariance matrix of \eqn{N_{ij}} (concatenated rowwise).
#'
#' See also (\insertCite{dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat}).
#'
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#'
#' @return The \eqn{k(k-1)/2 \times k(k-1)/2} covariance matrix of the differences of the off-diagonal cell counts \eqn{N_{ij}-N_{ji}} 
#' for \eqn{i,j=1,\ldots,k} and \eqn{i \ne j} in the NNCT, \code{ct} 
#'
#' @seealso \code{\link{var.nnsym}}, \code{\link{cov.tct}}, \code{\link{cov.nnct}} and \code{\link{cov.seg.coeff}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#'
#' cov.nnsym(covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.nnsym(covN)
#'
#' @export
cov.nnsym <- function(covN)
{
  k<-sqrt(nrow(covN))
  ent2<-cbind(rep(1:k,rep(k,k)),rep(1:k,k))
  ind.mat1<-ind.nnsym(k)
  
  kl<-k*(k-1)/2
  varcovT<-matrix(0,nrow=kl,ncol=kl)
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]
    c1<-rc1[2]
    ri1<-which(ent2[,1]==r1 & ent2[,2]==c1) #row indices for the N vector
    ri2<-which(ent2[,1]==c1 & ent2[,2]==r1)
    for (j in i:kl)
    {
      rc1<-ind.mat1[j,]
      r1n<-rc1[1]
      c1n<-rc1[2]
      ci1<-which(ent2[,1]==r1n & ent2[,2]==c1n) #column indices for the N vector
      ci2<-which(ent2[,1]==c1n & ent2[,2]==r1n) 
      
      varcovT[i,j]<-covN[ri1,ci1]-covN[ri1,ci2]-covN[ri2,ci1]+covN[ri2,ci2]
      varcovT[j,i]<-varcovT[i,j]
    }
  }
  
  varcovT
} #end for the function
#'

#################################################################

# funsXsq.nnsym.dx
#'
#' @title Dixon's NN Symmetry Test with Chi-square Approximation for multiple classes
#'
#' @description
#' Two functions: \code{Xsq.nnsym.dx.ct} and \code{Xsq.nnsym.dx}.
#' 
#' Both functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs the hypothesis test of equality of the expected value of the off-diagonal 
#' cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k \ge 2} classes.
#' That is, each performs Dixon's overall NN symmetry test.
#' The test is appropriate (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' Each symmetry test is based on the chi-squared approximation of the corresponding quadratic form
#' and is an extension of Dixon's NN symmetry test, which is extended by
#' \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat}.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)/2}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of differences of the off-diagonal entries,(which is
#' 0 for this function) and also the sample estimates (i.e. observed values) of absolute differences of the off-diagonal entries of 
#' NNCT (in the upper-triangular form).
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=E(N_{ji})} entries for all \eqn{i \ne j} (i.e., symmetry in the 
#' mixed NN structure).
#' 
#' See also
#' (\insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Xsq.nnsym.dx.ct} only
#' @param covS The \eqn{k(k-1)/2 \times k(k-1)/2} covariance matrix of the differences of the off-diagonal entries in the NNCT,
#' \code{ct}, usually the output of the function \code{\link{cov.nnsym}}.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Xsq.nnsym.dx} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.nnsym.dx} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Xsq.nnsym.dx} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for Dixon's overall NN symmetry test}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)/2} for this function.}
#' \item{estimate}{Estimates, i.e., absolute differences of the off-diagonal entries of 
#' NNCT (in the upper-triangular form).}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null values for the differences between the expected values of the off-diagonal 
#' entries, which is 0 for this function.}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Xsq.nnsym.dx.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.nnsym.dx} only}
#'  
#' @seealso \code{\link{Znnsym.dx.ct}}, \code{\link{Znnsym.dx}}, \code{\link{Znnsym}},
#' \code{\link{Xsq.nnsym}}, \code{\link{Xsq.nnsym.ss.ct}}, \code{\link{Xsq.nnsym.ss}}
#' and \code{\link{Qsym.test}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsXsq.nnsym.dx
NULL
#'
#' @rdname funsXsq.nnsym.dx
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
#' covS<-cov.nnsym(covN)
#'
#' Xsq.nnsym.dx(Y,cls)
#' Xsq.nnsym.dx.ct(ct,covS)
#'
#' Xsq.nnsym.dx(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Xsq.nnsym.dx(Y,fcls)
#' Xsq.nnsym.dx.ct(ct,covS)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#' covS<-cov.nnsym(covN)
#'
#' Xsq.nnsym.dx(Y,cls)
#' Xsq.nnsym.dx.ct(ct,covS)
#'
#' @export
Xsq.nnsym.dx.ct <- function(ct,covS)
{
  k<-nrow(ct)
  ind.mat1<-ind.nnsym(k)
  
  kl<-nu<-k*(k-1)/2
  Tvec<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]
    c1<-rc1[2]
    Tvec<-c(Tvec,ct[r1,c1]-ct[c1,r1]) #entry in ct corresponding to index i of T
  }
  
  ts<- t(Tvec) %*% ginv(covS,tol=1.490116e-20) %*% (Tvec)
  pval<- 1-pchisq(ts,df=nu)
  
  method <-"Dixon's Test of NN symmetry"
  
  dname <-deparse(substitute(ct))
  
  diff.mat<-abs(ct[upper.tri(ct)]- t(ct)[upper.tri(ct)])
  estimate<-matrix(0,k,k)
  estimate[upper.tri(estimate, diag=FALSE)]<-diff.mat 
  
  estimate.name <-c("absolute differences of the off-diagonal entries of NNCT")
  estimate.name2 <-c("absolute differences of the off-diagonal entries of NNCT (in the upper-triangular form)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = 0,
    method = method,
    data.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsXsq.nnsym.dx
#'
#' @export
Xsq.nnsym.dx <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
  covS<-cov.nnsym(covN)
  
  rval<-Xsq.nnsym.dx.ct(ct,covS)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

#' @title Overall NN Symmetry Test with Chi-square Approximation
#'
#' @description
#' An object of class \code{"Chisqtest"} performing the hypothesis test of equality of the expected
#' values of the off-diagonal cell counts (i.e., entries) under RL or CSR in the NNCT for \eqn{k \ge 2} classes.
#' That is, the test performs Dixon's or Pielou's (first type of) overall NN symmetry test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data or for sparsely sample data, respectively.
#' (See \insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' The \code{type="dixon"} refers to Dixon's overall NN symmetry test and 
#' \code{type="pielou"} refers to Pielou's first type of overall NN symmetry test.
#' The symmetry test is based on the chi-squared approximation of the corresponding quadratic form
#' and \code{type="dixon"} yields an extension of Dixon's NN symmetry test, which is extended by
#' \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} and \code{type="pielou"} yields
#' Pielou's overall NN symmetry test.
#'  
#' The function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k-1)/2}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of differences of the off-diagonal entries,(which is
#' 0 for this function) and also the sample estimates (i.e. observed values) of absolute differences of the off-diagonal entries of 
#' NNCT (in the upper-triangular form).
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=E(N_{ji})} for \eqn{i \ne j} in the \eqn{k \times k} NNCT (i.e., symmetry in the 
#' mixed NN structure) for \eqn{k \ge 2}.
#' In the output, if if \code{type="pielou"}, 
#' the test statistic, \eqn{p}-value and the df are valid only for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,dixon:1994,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param type The type of the overall NN symmetry test with default=\code{"dixon"}.
#' Takes on values \code{"dixon"} and \code{"pielou"} for Dixon's and Pielou's (first type) overall NN symmetry test
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function
#'   
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for Dixon's or Pielou's (first type of)
#' overall NN symmetry test}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k-1)/2} for this function.}
#' \item{estimate}{Estimates, i.e., absolute differences of the off-diagonal entries of 
#' NNCT (in the upper-triangular form).}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null values for the differences between the expected values of the off-diagonal 
#' entries, which is 0 for this function.}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}, or name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link{Znnsym.ss}}, \code{\link{Znnsym.dx}} and \code{\link{Znnsym2cl}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Xsq.nnsym(Y,cls)
#' Xsq.nnsym(Y,cls,method="max")
#' Xsq.nnsym(Y,cls,type="pielou")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#'
#' Xsq.nnsym(Y,fcls)
#' Xsq.nnsym(Y,fcls,type="pielou")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Xsq.nnsym(Y,cls)
#' Xsq.nnsym(Y,cls,type="pielou")
#'
#' @export
Xsq.nnsym <- function(dat,lab,type="dixon",...)
{
  if ((type %in% c("dixon","pielou"))==FALSE)
  {stop("type is misspecified, should be one of dixon or pielou (in quotes)")}
  
  if (type == "dixon")
  {
    res<- Xsq.nnsym.dx(dat,lab,...)
  } else
  {
    res<- Xsq.nnsym.ss(dat,lab,...)
  }
  
  return(res)
}
#'

#################################################################

#' @title The Shared NN Vectors for Multiple Classes
#'
#' @description 
#' Returns a \code{matrix} with \eqn{k} rows where each row is the vector of number of points with shared NNs,
#' \eqn{Q_i=(Q_{i0},Q_{i1},\ldots)} where \eqn{Q_{ij}} is the number of class \eqn{i} points that are NN to class \eqn{j} points.
#' The function also returns the indices of columns with nonzero sums as a vector.
#' 
#' The output matrix of shared NNs is used in testing symmetry in shared NN structure (i.e. \eqn{Q}-symmetry
#' or Pielou's second type of symmetry), 
#' e.g., in functions \code{\link{Qsym.ct}} and \code{\link{Qsym.test}}.
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix (IPD matrix), otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return \code{Qval} returns the \eqn{Q} value
#' \code{Qvec} returns a \code{list} with two elements
#'  \item{q}{the \eqn{Q} value, the number of shared NNs}
#'  \item{qvec}{the \code{vector} of \eqn{Q_j} values} 
#' \code{sharedNN} returns a \code{matrix} with 2 rows, where first row is the \eqn{j} values and second row is
#' the corresponding vector of \eqn{Q_j} values
#' \code{Rval}{the \eqn{R} value, the number of reflexive NNs}
#' 
#' @return Returns a \code{list} with two elements
#'  \item{Nv}{A \eqn{k}-row matrix of shared NNs by class where each row of the matrix is the vector of number of
#'  points with shared NNs \eqn{Q_i=(Q_{i0},Q_{i1},\ldots)} where \eqn{Q_{ij}} is the number of class \eqn{i} points that are NN
#'  to \eqn{j} points.}
#'  \item{col.ind}{The \code{vector} of indices of columns with nonzero sums} 
#' 
#' @seealso \code{\link{Qval}}, \code{\link{Qvec}} and \code{\link{sharedNN}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' #'
#' sharedNNmc(ipd,cls)
#' sharedNNmc(Y,cls,is.ipd = FALSE)
#' sharedNNmc(Y,cls,is.ipd = FALSE,method="max")
#' #'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' sharedNNmc(ipd,fcls)
#' #'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' #'
#' sharedNNmc(ipd,cls) 
#'
#' @export
sharedNNmc <- function(x,lab,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  n<-nrow(ipd)
  
  flab<-as.factor(lab) #converting class labels to factors
  clnames<-levels(flab)
  k<-length(clnames)
  
  lab<-as.numeric(flab)
  
  invnn <- matrix(0,k,n); 
  ind <- rep(0,n)
  
  for (i in 1:n)
    ind[i] <- NN(ipd,i);
  
  for (j in 1:n)
  {
    invnn[lab[ind[j]],ind[j]] <- invnn[lab[ind[j]],ind[j]]+1
  }
  
  N<-matrix(0,k,100) #should have been 7
  for (i in 1:k)
    for (j in 1:n)
    {
      if (lab[j]==i)
        N[i,invnn[i,j]+1] <- N[i,invnn[i,j]+1]+1
    }
  col.ind<-which(col.sum(N)!=0)
  Nv<-matrix(N[,col.ind],nrow=k)
  q<-length(col.ind) 
  colnm<-vector()
  for (i in 1:q)
  {
    colnm<-c(colnm,paste("Q",col.ind[i]-1,sep=""))
  }
  colnames(Nv)<-colnm
  rownames(Nv)<-clnames
  list(Nv=Nv,col.ind=col.ind)
} #end for the function
#'

#################################################################

#' @title \eqn{Q}-symmetry Contingency Table (QCT)
#'
#' @description
#' Returns the \eqn{k \times 3} contingency table for \eqn{Q}-symmetry (i.e. \eqn{Q}-symmetry contingency table (QCT)) given the 
#' IPD matrix or data set \code{x} where \eqn{k} is the number of classes in the data set.
#' Each row in the QCT is the vector of number of points with shared NNs,
#' \eqn{Q_i=(Q_{i0},Q_{i1},Q_{i2})} where \eqn{Q_{ij}} is the number of class \eqn{i} points that are NN to class \eqn{j} points
#' for \eqn{j=0,1} and \eqn{Q_{i2}} is the number of class \eqn{i} points that are NN to class \eqn{j} or more points.
#' That is, this function pools the cells 3 or larger together for \eqn{k} classes, so \eqn{Q_2}, \eqn{Q_3} etc. are pooled,
#' so the column labels are \eqn{Q_0}, \eqn{Q_1} and \eqn{Q_2} with the last one is actually sum of \eqn{Q_j} for \eqn{j \ge 2}.
#' Rows the QCT are labeled with the corresponding class labels.
#' 
#' \eqn{Q}-symmetry is also equivalent to Pielou's second type of NN symmetry
#' or the symmetry in the shared NN structure for all classes.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#'
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#'
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns the \eqn{k \times 3} QCT where \eqn{k} is the number of classes in the data set.
#'
#' @seealso \code{\link{sharedNNmc}}, \code{\link{Qsym.test}} and \code{\link{scct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(n*3),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' Qsym.ct(ipd,cls)
#' Qsym.ct(Y,cls,is.ipd = FALSE)
#' Qsym.ct(Y,cls,is.ipd = FALSE,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Qsym.ct(ipd,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' Qsym.ct(ipd,cls)
#'
#' @export
Qsym.ct <- function(x,lab,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  n<-nrow(ipd)
  
  flab<-as.factor(lab) #converting class labels to factors
  clnames<-levels(flab)
  k<-length(clnames)
  
  Nsh<-sharedNNmc(ipd,lab,...)
  CT0<-Nsh$Nv
  colind<-Nsh$col.ind
  
  nc<-length(colind)
  if (nc>3)
  {
    CT<-matrix(0,k,3)
    CT[1:k,1:2]<-CT0[1:k,1:2]
    
    for (i in 1:k)
      CT[i,3]<-sum(CT0[i,-(1:2)]);
    
    rownames(CT)<-clnames 
    colnm<-c(paste("Q",colind[1]-1,sep=""),paste("Q",colind[2]-1,sep=""),paste("Q",colind[3]-1,"+",sep=""))
    
    colnames(CT)<-colnm
  } else
  {CT<-CT0}
  CT
} #end for the function
#'

################################################################# 

#' @title Pielou's Second Type of NN Symmetry Test with Chi-square Approximation
#'
#' @description
#' An object of class \code{"Chisqtest"} performing the hypothesis test of equality of the probabilities for the rows
#' in the \eqn{Q}-symmetry contingency table (QCT).
#' Each row of the QCT is the vector of Qi\eqn{j} values where \eqn{Q_{ij}} is the number of class \eqn{i} points that are NN
#' to \eqn{j} points.
#' That is, the test performs Pielou's second type of NN symmetry test which is also equivalent to Pearson's
#' test on the QCT (\insertCite{pielou:1961;textual}{nnspat}).
#' Pielou's second type of NN symmetry is the symmetry in the shared NN structure for all classes, which is also 
#' called \eqn{Q}-symmetry.
#' The test is appropriate (i.e. have the appropriate asymptotic sampling distribution)
#' provided that data is obtained by sparse sampling, although simulations suggest it seems to work for
#' completely mapped data as well.
#' (See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more detail).
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#' 
#' The argument combine is a logical argument (default=\code{TRUE}) to determine whether to combine the 3rd column 
#' and the columns to the left.
#' If \code{TRUE}, this function pools the cells 3 or larger together for \eqn{k} classes in the QCT, 
#' so \eqn{Q_2}, \eqn{Q_3} etc. are pooled, so the column
#' labels are \eqn{Q_0}, \eqn{Q_1} and \eqn{Q_2} with the last one is actually sum of \eqn{Q_j} for \eqn{j \ge 2} in the QCT.
#' If \code{FALSE}, the function does not perform the pooling of the cells.
#' 
#' The function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{(k-1)(n_c-1)} where \eqn{n_c} is the number of
#' columns in QCT (which reduces to \eqn{2(k-1)}, if \code{combine=TRUE}). It also provides the description of
#' the alternative with the corresponding null values (i.e. expected values) of the entries of the QCT and also the sample estimates 
#' of the entries of QCT (i.e., the observed QCT).
#' The function also provides names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is the symmetry in the shared NN structure for each class, that is,
#' all \eqn{E(Q_{ij})=n_i Q_j/n} where \eqn{n_i} the size of class \eqn{i} and \eqn{Q_j} is the sum of column \eqn{j} 
#' in the QCT (i.e., the total number of points serving as NN to class \eqn{j} other points). (i.e., symmetry in the 
#' mixed NN structure).
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:SWJ-spat-sym2014;textual}{nnspat})
#' and the references therein.
#' 
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix (IPD matrix), otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param combine A logical parameter (default=\code{TRUE}). If \code{TRUE}, 
#' the cells in column 3 or columns to the left are merged in the QCT, so \eqn{Q_2}, \eqn{Q_3} etc. are pooled, so the column
#' labels are \eqn{Q_0}, \eqn{Q_1} and \eqn{Q_2} with the last one is actually sum of \eqn{Q_j} for \eqn{j \ge 2} in the QCT. 
#' If \code{FALSE}, the function does not perform the pooling of the cells.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for Pielou's second type of NN symmetry test (i.e., \eqn{Q}-symmetry 
#' which is equivalent to symmetry in the shared NN structure)}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{(k-1)(n_c-1)} where \eqn{n_c} is the number of
#' columns in QCT (which reduces to \eqn{2(k-1)} if \code{combine=TRUE}).}
#' \item{estimate}{Estimates, i.e., the observed QCT.}
#' \item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
#' \item{null.value}{Hypothesized null values for the entries of the QCT, i.e., the matrix with entries 
#' \eqn{E(Q_{ij})=n_i Q_j/n} where \eqn{n_i} the size of class \eqn{i} and \eqn{Q_j} is the sum of column \eqn{j} in the QCT (i.e., the total
#' number of points serving as NN to class \eqn{j} other points).}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{x}}
#'  
#' @seealso \code{\link{Znnsym}} and \code{\link{Xsq.nnsym}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' Qsym.ct(ipd,cls)
#'
#' Qsym.test(ipd,cls)
#' Qsym.test(Y,cls,is.ipd = FALSE)
#' Qsym.test(Y,cls,is.ipd = FALSE,method="max")
#'
#' Qsym.test(ipd,cls,combine = FALSE)
#'
#' #cls as a faqctor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Qsym.test(ipd,fcls)
#' Qsym.test(Y,fcls,is.ipd = FALSE)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' Qsym.test(ipd,cls)
#' Qsym.test(Y,cls,is.ipd = FALSE)
#'
#' @export
Qsym.test <- function(x,lab,is.ipd=TRUE,combine=TRUE,...)
{
  ifelse(combine,qct<-Qsym.ct(x,lab,is.ipd,...) ,qct<-sharedNNmc(x,lab,is.ipd,...)$Nv)
  
  TS<-chisq.test(qct)
  ts<-TS$stat
  nu<-TS$par
  pval<-TS$p.val
  
  ifelse(combine,method <-"Pielou's Test of Q Symmetry (i.e., Symmetry in Shared NN structure)\n 
        with the entries in columns Qk for k>=2 (if exist) are combined in the Q symmetry contingency table",
         method <-"Pielou's Test of Q Symmetry (i.e., Symmetry in Shared NN structure)\n 
        with the entire Q symmetry contingency table")
  
  dname <-TS$data.name
  
  estimate<-qct
  estimate.name <-c("Q symmetry Contingency Table (QsymCT) entries")
  
  EN<-TS$exp
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value = EN,
    method = method,
    data.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'

#################################################################

#' @title Reflexivity Contingency Table (RCT)
#'
#' @description
#' Returns the RCT given the IPD matrix or data set \code{x}, the RCT is \eqn{2 \times 2} regardless of the 
#' number of classes in the data set.
#' 
#' RCT is constructed by categorizing the NN pairs according to pair type as self or mixed and whether
#' the pair is reflexive or non-reflexive. 
#' A base-NN pair is called a reflexive pair, if the elements of the pair are NN to each other;
#' a non-reflexive pair, if the elements of the pair are not NN to each other;
#' a self pair, if the elements of the pair are from the same class; a mixed pair, if the
#' elements of the pair are from different classes.
#' Row labels in the RCT are \code{"ref"} for reflexive and \code{"non-ref"} for non-reflexive and 
#' column labels are \code{"self"} and \code{"mixed"}.
#' 
#' The argument \code{is.ipd} is a logical argument (default=\code{TRUE}) to determine the structure of the argument \code{x}.
#' If \code{TRUE}, \code{x} is taken to be the inter-point distance (IPD) matrix, and if \code{FALSE}, \code{x} is taken to be the data set
#' with rows representing the data points.
#'
#' See also (\insertCite{ceyhan:NNreflexivity2017,ceyhan:NNreflex1D2018;textual}{nnspat})
#' and the references therein.
#'
#' @param x The IPD matrix (if \code{is.ipd=TRUE}) or a data set of points in matrix or data frame form where points
#' correspond to the rows (if \code{is.ipd = FALSEALSE}).
#' @param lab The \code{vector} of class labels (numerical or categorical)
#' @param is.ipd A logical parameter (default=\code{TRUE}). If \code{TRUE}, \code{x} is taken as the inter-point distance
#' matrix, otherwise, \code{x} is taken as the data set with rows representing the data points. 
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Returns the \eqn{2 \times 2} RCT, see the description above for more detail.
#'
#' @seealso \code{\link{nnct}}, \code{\link{tct}} and \code{\link{scct}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' rct(ipd,cls)
#' rct(Y,cls,is.ipd = FALSE)
#' rct(Y,cls,is.ipd = FALSE,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' rct(ipd,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' rct(ipd,cls)
#'
#' @export
rct <- function(x,lab,is.ipd=TRUE,...)
{
  ifelse(is.ipd,ipd<-x,ipd<-ipd.mat(x,...))
  n<-nrow(ipd)
  if (n<=1)
  {stop('n must be >=2 for the reflexivity contingency table, rct, to be defined')}
  
  Rf <- 0; SRf<-0; SNR<-0;
  nn.ind <- list()
  Lind <- vector()
  
  for (i in 1:n)
  {
    nn.pti<-NN(ipd,i)#,...)
    nn.ind[[i]] <- nn.pti;
    Lind<-c(Lind,length(nn.pti))
  }
  
  Nnn<-sum(Lind)
  
  for (j in 1:n)
  {
    for (k in 1:Lind[j])
    {
      if (sum(j==nn.ind[[ nn.ind[[j]][k] ]]) == 1)
      {Rf<- Rf +(1/Lind[j])*(1/Lind[nn.ind[[j]][k]]);
      if (lab[j]==lab[nn.ind[[j]][k]])
        SRf<-SRf +(1/Lind[j])*(1/Lind[nn.ind[[j]][k]])
      }
      else
      {if (lab[j]==lab[nn.ind[[j]][k]])
        SNR<-SNR+1/Lind[j]
      }
    }
  }
  rfct<-matrix(c(SRf,Rf-SRf,SNR,n-Rf-SNR),nrow=2,byrow="T")
  
  rownm<-c("ref","non-ref")
  rownames(rfct)<-rownm
  colnm<-c("self","mixed")
  colnames(rfct)<-colnm
  
  rfct
} #end for the function
#'

#################################################################

#' @title Expected Values of the Cell Counts in RCT
#'
#' @description Returns a \code{matrix} of same dimension as the RCT, \code{rfct}, 
#' whose entries are the expected cell counts of
#' the RCT under RL or CSR.
#' 
#' See also (\insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}).
#'
#' @param rfct An RCT
#' @param nvec The \code{vector} of class sizes
#'
#' @return A \code{matrix} of the expected values of cell counts in the RCT.
#'
#' @seealso \code{\link{rct}}, \code{\link{EV.nnct}} and \code{\link{EV.tct}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#' EV.rct(rfct,nvec)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' nvec<-as.numeric(table(fcls))
#' rfct<-rct(ipd,fcls)
#' EV.rct(rfct,nvec)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' rfct<-rct(ipd,cls)
#' EV.rct(rfct,nvec)
#'
#' @export
EV.rct <- function(rfct,nvec)
{
  n<-sum(nvec)
  k<-length(nvec)
  
  if (k<4) 
  {nvec<-c(nvec,rep(0,4-k))} 
  
  R<-sum(rfct[1,]) #number of reflexive NNs
  Paa<-P11(nvec)
  ENsr<-R*Paa
  
  Pab<-P12(nvec)
  ENmn<-(n-R)*Pab
  
  erfct<-matrix(c(ENsr,R-ENsr,n-R-ENmn,ENmn),nrow=2,byrow="T")
  
  rownm<-c("ref","non-ref")
  rownames(erfct)<-rownm
  colnm<-c("self","mixed")
  colnames(erfct)<-colnm
  
  erfct
} #end for the function
#'

#################################################################

# funsZnnref
#'
#' @title Z Tests for NN Reflexivity
#'
#' @description
#' Two functions: \code{Znnref.ct} and \code{Znnref}.
#'
#' Both functions are objects of class \code{"refhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the
#' diagonal cell counts (i.e., entries) under RL or CSR in the RCT for \eqn{k \ge 2} classes.
#' That is, each test performs NN reflexivity test (i.e., a test of self reflexivity and a test of
#' mixed non-reflexivity, corresponding to entries \eqn{(1,1)} and \eqn{(2,2)}, respectively, in the RCT) which is
#' appropriate (i.e. have the appropriate asymptotic sampling distribution) for completely mapped data.
#' (See \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat} for more detail).
#' 
#' The reflexivity test is based on the normal approximation of the diagonal entries
#' in the RCT and are due to \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}.
#' 
#' Each function yields the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values (i.e. null value(s)), confidence intervals and sample estimates (i.e. observed values)for the
#' self reflexivity and mixed non-reflexivity values (i.e., entries \eqn{(1,1)} and \eqn{(2,2)} values, respectively)
#' in the RCT. Each function also gives names of the test statistics, null values and the method and the data
#' set used.
#' 
#' The null hypothesis is that \eqn{E(N_{11})=R P_{aa}} and \eqn{E(N_{22})=R P_{ab}} in the RCT, where \eqn{R} is the number of reflexive
#' NNs and \eqn{P_{aa}} is the probability of any two points selected are being from the same class
#' and \eqn{P_{ab}} is the probability of any two points selected are being from two different classes.
#' 
#' The \code{Znnref} functions (i.e. \code{Znnref.ct} and \code{Znnref}) are different from 
#' the \code{Znnself} functions (i.e. \code{\link{Znnself.ct}} and \code{\link{Znnself}}) and 
#' from \code{Zself.ref} functions (i.e. \code{\link{Zself.ref.ct}} and \code{\link{Zself.ref}}), and also
#' from \code{Znnself.sum} functions (i.e. \code{Znnself.sum.ct} and \code{Znnself.sum}).
#' \code{Znnref} functions are for testing the self reflexivity and mixed non-reflexivity
#' using the diagonal entries in the RCT while \code{Znnself} functions are testing the self reflexivity at a
#' class-specific level (i.e. for each class) using the first column in the SCCT, and
#' \code{Zself.ref} functions are for testing the self reflexivity for the entire data set
#' using entry \eqn{(1,1)} in RCT, and \code{Znnself.sum} functions are testing the cumulative species correspondence using
#' the sum of the self column (i.e., the first column) in the SCCT.
#' 
#' @param rfct An RCT, used in \code{Znnref.ct} only
#' @param nvec The \code{vector} of class sizes, used in \code{Znnref.ct} only
#' @param Qv The number of shared NNs, used in \code{Znnref.ct} only
#' @param Tv \eqn{T} value, which is the number of triplets \eqn{(z_i, z_j, z_k)} with 
#' "\eqn{NN(z_i) = NN(z_j) = z_k} and \eqn{NN(z_k) = z_j} where \eqn{NN(\cdot)} is the nearest neighbor function, used in \code{Znnref.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnref} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnref} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Znnref} only
#'
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistics for self reflexivity and mixed non-reflexivity, corresponding to entries
#' \eqn{(1,1)} and \eqn{(2,2)} in the RCT}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \eqn{p}-values for self reflexivity and mixed non-reflexivity tests}
#' \item{conf.int}{Confidence intervals for the self reflexivity and mixed non-reflexivity values
#' (i.e., diagonal entries \eqn{(1,1)} and \eqn{(2,2)} values, respectively) in the RCT at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{cnf.lvl}{Level of the onfidence intervals of the diagonal entries, provided in \code{conf.level}.}
#' \item{estimate}{Estimates of the parameters, i.e., the observed diagonal entries \eqn{(1,1)} and \eqn{(2,2)}
#' in the RCT, \code{rfct}.}
#' \item{null.value}{Hypothesized null values for the self reflexivity and mixed non-reflexivity values
#' (i.e., expected values of the diagonal entries \eqn{(1,1)} and \eqn{(2,2)} values, 
#' which are \eqn{E(N_{11})=R P_{aa}} and \eqn{E(N_{22})=R P_{ab}}, respectively) in the RCT.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{rfct}, returned by \code{Znnref.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Znnref} only}
#' 
#' @seealso \code{\link{Znnself.ct}}, \code{\link{Znnself}}, \code{\link{Zmixed.nonref.ct}},
#' \code{\link{Zmixed.nonref}}, \code{\link{Xsq.nnref.ct}} and \code{\link{Xsq.nnref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnref
NULL
#'
#' @rdname funsZnnref
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' Tv<-Tval(W,Rv)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Znnref(Y,cls)
#' Znnref(Y,cls,method="max")
#'
#' Znnref.ct(rfct,nvec,Qv,Tv)
#' Znnref.ct(rfct,nvec,Qv,Tv,alt="g")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' R<-Rval(W)
#' Tv<-Tval(W,R)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Znnref(Y,cls,alt="g")
#'
#' Znnref.ct(rfct,nvec,Qv,Tv)
#' Znnref.ct(rfct,nvec,Qv,Tv,alt="l")
#'
#' @export
Znnref.ct <- function(rfct,nvec,Qv,Tv,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-diag(rfct)
  
  clnames<-c("self-ref","mixed-non-ref") #names of the tests for the diagonal entries of rfct
  names(estimate)<-clnames
  
  n<-sum(nvec)
  k<-length(nvec)
  
  if (k<4) 
  {nvec<-c(nvec,rep(0,4-k))} 
  
  Nsr<-rfct[1,1]
  Nmn<-rfct[2,2]
  
  R<-sum(rfct[1,]) #row sum for rfct
  Paa<-P11(nvec)
  ENsr<-R*Paa
  Paaaa<-P1111(nvec)
  Paabb<-P1122(nvec)
  VNsr<-R^2*(Paaaa+Paabb-Paa^2)+2*R*(Paa-Paaaa-Paabb)
  se.VNsr<-sqrt(VNsr)
  Zsr<-(Nsr-ENsr)/se.VNsr
  
  Pab<-P12(nvec)
  ENmn<-(n-R)*Pab
  Paabc<-P1123(nvec)
  Pabcd<-P1234(nvec)
  Paab<-P112(nvec)
  Pabc<-P123(nvec)
  VNmn<-(n-R)^2*(2*Paabb+4*Paabc+Pabcd-Pab^2)+(n-R)*Pab+(2*n-2*R+Qv-4*Tv)*(Paab+Pabc)+(-3*n+3*R-Qv+4*Tv)*(2*Paabb+4*Paabc+Pabcd)
  se.VNmn<-sqrt(VNmn)
  Zmn<-(Nmn-ENmn)/se.VNmn
  ts<-c(Zsr,Zmn)
  stderr<-c(se.VNsr,se.VNmn)
  
  if (all(is.na(ts)))
  {stop('Both test statistics are NaN, so the Reflexivity Z-tests are not defined')}
  
  null.val<-c(ENsr,ENmn)
  names(null.val)<-clnames
  names.null <-"diagonal rct entries"
  
  cint<-matrix(0,ncol=2,nrow=2)

  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint[1,] <-estimate[1]+c(-Inf, qnorm(conf.level))*stderr[1]
           cint[2,] <-estimate[2]+c(-Inf, qnorm(conf.level))*stderr[2]
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint[1,] <-estimate[1]+c(-qnorm(conf.level),Inf)*stderr[1]
           cint[2,] <-estimate[2]+c(-qnorm(conf.level),Inf)*stderr[2]
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           Cint <-qnorm(1 - alpha/2)
           cint[1,] <-estimate[1]+c(-Cint, Cint)*stderr[1]
           cint[2,] <-estimate[2]+c(-Cint, Cint)*stderr[2]
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-"Z-tests for NN Reflexivity Patterns"
  ts.names <-"Reflexivity Z-tests (for Diagonals of Reflexivity Contingency Table (rct)) for"
  
  dname <-deparse(substitute(rfct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    conf.int = cint,
    cnf.lvl=conf.level,
    estimate = estimate,
    null.value = null.val,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  class(rval) <- "refhtest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZnnref
#'
#' @export
Znnref <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  R<-Rval(W)
  Tv<-Tval(W,R)
  
  nvec<-as.numeric(table(lab))
  rfct<-rct(ipd,lab,...)
  
  rval<-Znnref.ct(rfct,nvec,Qv,Tv,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZself.ref
#'
#' @title Self-Reflexivity Test with Normal Approximation
#'
#' @description
#' Two functions: \code{Zself.ref.ct} and \code{Zself.ref}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of self reflexivity in the NN structure using the
#' number of self-reflexive NN pairs (i.e. the first diagonal entry, \eqn{(1,1)}) in the RCT for \eqn{k \ge 2} classes.
#' That is, each test performs a test of self reflexivity corresponding to entry \eqn{(1,1)} in the RCT)
#' which is appropriate (i.e. have the appropriate asymptotic sampling distribution) for completely mapped data.
#' (See \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat} for more detail).
#' 
#' The self reflexivity test is based on the normal approximation of the diagonal entry \eqn{(1,1)}
#' in the RCT and are due to \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, sample estimate (i.e. observed value) and null (i.e., expected) value for the
#' self reflexivity value (i.e., diagonal entry \eqn{(1,1)} value, respectively) in the RCT, 
#' and method and name of the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{11})=R P_{aa}} in the RCT, where \eqn{R} is the number of reflexive
#' NNs and \eqn{P_{aa}} is the probability of any two points selected are being from the same class.
#' 
#' The \code{Zself.ref} functions (i.e. \code{Zself.ref.ct} and \code{Zself.ref}) are different from the \code{Znnref}
#' functions (i.e. \code{\link{Znnref.ct}} and \code{\link{Znnref}}) and from \code{Znnself} functions (i.e. \code{\link{Znnself.ct}} and \code{\link{Znnself}}), and also
#' from \code{Znnself.sum} functions (i.e. \code{Znnself.sum.ct} and \code{Znnself.sum}).
#' \code{Zself.ref} functions are for testing the self reflexivity for the entire data set
#' using entry \eqn{(1,1)} in RCT while \code{Znnself} functions are testing the self reflexivity at a class-specific level
#' (i.e. for each class) using the first column in the SCCT, \code{Znnref} functions are for testing the self
#' reflexivity and mixed non-reflexivity using the diagonal entries in the RCT, and
#' \code{Znnself.sum} functions are testing the cumulative species correspondence using the sum of the self column (i.e.,
#' the first column) in the SCCT.
#' 
#' @param rfct An RCT, used in \code{Zself.ref.ct} only
#' @param nvec The \code{vector} of class sizes, used in \code{Zself.ref.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zself.ref} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zself.ref} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Zself.ref} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for self reflexivity corresponding to entry \eqn{(1,1)} in the RCT}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the self reflexivity value (i.e., diagonal entry \eqn{(1,1)} value)
#' in the RCT at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the observed diagonal entry \eqn{(1,1)} in the RCT, \code{rfct}.}
#' \item{null.value}{Hypothesized null value for the self reflexivity value (i.e., expected value of the 
#' diagonal entry \eqn{(1,1)} which is \eqn{E(N_{11})=R P_{aa}}) in the RCT.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{rfct}, returned by \code{Zself.ref.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zself.ref} only}
#'  
#' @seealso \code{\link{Znnref.ct}}, \code{\link{Znnref}}, \code{\link{Zmixed.nonref.ct}} and
#' \code{\link{Zmixed.nonref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZself.ref
NULL
#'
#' @rdname funsZself.ref
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Zself.ref(Y,cls)
#' Zself.ref(Y,cls,method="max")
#'
#' Zself.ref.ct(rfct,nvec)
#' Zself.ref.ct(rfct,nvec,alt="g")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Zself.ref(Y,cls,alt="g")
#'
#' Zself.ref.ct(rfct,nvec)
#' Zself.ref.ct(rfct,nvec,alt="l")
#'
#' @export
Zself.ref.ct <- function(rfct,nvec,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-Nsr<-rfct[1,1]
  names(estimate) <-"number of self-reflexive NN pairs"
  
  n<-sum(nvec)
  k<-length(nvec)
  
  if (k<4) 
  {nvec<-c(nvec,rep(0,4-k))} 
  
  R<-sum(rfct[1,]) #row sum for rfct
  Paa<-P11(nvec)
  ENsr<-R*Paa
  Paaaa<-P1111(nvec)
  Paabb<-P1122(nvec)
  VNsr<-R^2*(Paaaa+Paabb-Paa^2)+2*R*(Paa-Paaaa-Paabb)
  stderr<-sqrt(VNsr)
  ts<-(Nsr-ENsr)/stderr
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the self reflexivity test is not defined')}
  
  null.val<-ENsr
  names(null.val)<- "(expected) number of self-reflexive NN pairs"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )

  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  method <-"Z-test for Self-Reflexive NN Pattern"
  names(ts) <-"Self-Reflexivity Z-test"
  
  dname <-c("contingency table = ",deparse(substitute(rfct)))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZself.ref
#'
#' @export
Zself.ref <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  nvec<-as.numeric(table(lab))
  rfct<-rct(dat,lab,is.ipd = FALSE,...)
  
  rval<-Zself.ref.ct(rfct,nvec,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZmixed.nonref
#'
#' @title Mixed-Non-Reflexivity Test with Normal Approximation
#'
#' @description
#' Two functions: \code{Zmixed.nonref.ct} and \code{Zmixed.nonref}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of mixed non-reflexivity in the NN structure using the
#' number of mixed-non-reflexive NN pairs (i.e. the second diagonal entry, \eqn{(2,2)}) in the RCT for \eqn{k \ge 2} classes.
#' That is, each test performs a test of mixed non-reflexivity corresponding to entry \eqn{(2,2)} in the RCT)
#' which is appropriate (i.e. have the appropriate asymptotic sampling distribution) for completely mapped data.
#' (See \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat} for more detail).
#' 
#' The mixed non-reflexivity test is based on the normal approximation of the diagonal entry \eqn{(2,2)}
#' in the RCT and are due to \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, sample estimate (i.e. observed value) and null (i.e., expected) value for the
#' mixed non-reflexivity value (i.e., diagonal entry \eqn{(2,2)} value, respectively) in the RCT, 
#' and method and name of the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{22})=R P_{ab}} in the RCT, where \eqn{R} is the number of reflexive
#' NNs and \eqn{P_{ab}} is the probability of any two points selected are being from two different classes.
#' 
#' @param rfct An RCT, used in \code{Zmixed.nonref.ct} only
#' @param nvec The \code{vector} of class sizes, used in \code{Zmixed.nonref.ct} only
#' @param Qv The number of shared NNs, used in \code{Zmixed.nonref.ct} only
#' @param Tv \eqn{T} value, which is the number of triplets \eqn{(z_i, z_j, z_k)} with 
#' "\eqn{NN(z_i) = NN(z_j) = z_k} and \eqn{NN(z_k) = z_j} where \eqn{NN(\cdot)} is the nearest neighbor function, used in \code{Zmixed.nonref.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zmixed.nonref} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zmixed.nonref} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference of the off-diagonal entries, \eqn{N_{12}-N_{21}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Zmixed.nonref} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for mixed non-reflexivity corresponding to entry \eqn{(2,2)} in the RCT}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the mixed non-reflexivity value (i.e., diagonal entry \eqn{(2,2)} value)
#' in the RCT at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the observed diagonal entry \eqn{(2,2)} in the RCT, \code{rfct}.}
#' \item{null.value}{Hypothesized null value for the mixed non-reflexivity value (i.e., expected value of the 
#' diagonal entry \eqn{(2,2)} which is \eqn{E(N_{22})=R P_{ab}}) in the RCT.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{rfct}, returned by \code{Zmixed.nonref.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zmixed.nonref} only}
#'  
#' @seealso \code{\link{Zself.ref.ct}}, \code{\link{Zself.ref}}, \code{\link{Znnref.ct}} and
#' \code{\link{Znnref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZmixed.nonref
NULL
#'
#' @rdname funsZmixed.nonref
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' R<-Rval(W)
#' Tv<-Tval(W,R)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Zmixed.nonref(Y,cls)
#' Zmixed.nonref.ct(rfct,nvec,Qv,Tv)
#' Zmixed.nonref(Y,cls,alt="g")
#'
#' Zmixed.nonref(Y,cls,method="max")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' R<-Rval(W)
#' Tv<-Tval(W,R)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Zmixed.nonref(Y,cls,alt="g")
#'
#' Zmixed.nonref.ct(rfct,nvec,Qv,Tv)
#' Zmixed.nonref.ct(rfct,nvec,Qv,Tv,alt="l")
#'
#' @export
Zmixed.nonref.ct <- function(rfct,nvec,Qv,Tv,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-Nmn<-rfct[2,2]
  names(estimate) <-"number of mixed non-reflexive NN pairs"
  
  n<-sum(nvec)
  k<-length(nvec)
  
  if (k<4) 
  {nvec<-c(nvec,rep(0,4-k))} 
  
  R<-sum(rfct[1,]) #row sum for rfct
  Pab<-P12(nvec)
  ENmn<-(n-R)*Pab
  Paabb<-P1122(nvec)
  Paabc<-P1123(nvec)
  Pabcd<-P1234(nvec)
  Paab<-P112(nvec)
  Pabc<-P123(nvec)
  VNmn<-(n-R)^2*(2*Paabb+4*Paabc+Pabcd-Pab^2)+(n-R)*Pab+(2*n-2*R+Qv-4*Tv)*(Paab+Pabc)+
   (-3*n+3*R-Qv+4*Tv)*(2*Paabb+4*Paabc+Pabcd)
  stderr<-sqrt(VNmn)
  ts<-(Nmn-ENmn)/stderr
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the mixed non-reflexivity test is not defined')}
  
  null.val<-ENmn
  names(null.val)<- "(expected) number of mixed non-reflexive NN pairs"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  method <-"Z-test for Mixed Non-Reflexive NN Pattern"
  names(ts) <-"Mixed Non-Reflexivity Z-test"
  
  dname <-c("contingency table = ",deparse(substitute(rfct)))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZmixed.nonref
#'
#' @export
Zmixed.nonref <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  R<-Rval(W)
  Tv<-Tval(W,R)
  
  nvec<-as.numeric(table(lab))
  rfct<-rct(ipd,lab,...)
  
  rval<-Zmixed.nonref.ct(rfct,nvec,Qv,Tv,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsXsq.nnref
#'
#' @title Reflexivity Test with Chi-square Approximation
#'
#' @description
#' Two functions: \code{Xsq.nnref.ct} and \code{Xsq.nnref}.
#' 
#' Both functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the
#' diagonal cell counts (i.e., entries) under RL or CSR in the RCT for \eqn{k \ge 2} classes.
#' That is, each test performs an overall NN reflexivity test (for the vector of entries \eqn{(1,1)} and \eqn{(2,2)},
#' respectively, in the RCT) which is
#' appropriate (i.e. have the appropriate asymptotic sampling distribution) for completely mapped data.
#' (See \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat} for more detail).
#' 
#' Each reflexivity test is based on the chi-squared approximation of the corresponding quadratic form
#' for the vector of diagonal entries 
#' in the RCT and are due to \insertCite{ceyhan:NNreflexivity2017;textual}{nnspat}.
#'
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is 2, description of the 
#' alternative with the corresponding null values (i.e. expected values) of the diagonal entries
#' and also the sample estimates (i.e. observed values) of the diagonal entries of RCT (as a vector).
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{11},N_{22})=(R P_{aa},R P_{ab})} in the RCT, where \eqn{R} is the number of reflexive
#' NNs and \eqn{P_{aa}} is the probability of any two points selected are being from the same class
#' and \eqn{P_{ab}} is the probability of any two points selected are being from two different classes.
#' 
#' @param rfct An RCT, used in \code{Xsq.nnref.ct} only
#' @param nvec The \code{vector} of class sizes, used in \code{Xsq.nnref.ct} only
#' @param Qv The number of shared NNs, used in \code{Xsq.nnref.ct} only
#' @param Tv \eqn{T} value, which is the number of triplets \eqn{(z_i, z_j, z_k)} with 
#' \eqn{NN(z_i) = NN(z_j) = z_k} and \eqn{NN(z_k) = z_j} where \eqn{NN(\cdot)} is the nearest neighbor function, used in \code{Xsq.nnref.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Xsq.nnref} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.nnref} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{Xsq.nnref} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for overall NN reflexivity test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is 2 for this function.}
#' \item{estimate}{Estimates of the parameters, i.e., the observed diagonal entries \eqn{(1,1)} and \eqn{(2,2)}
#' in the RCT, \code{rfct}.}
#' \item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
#' \item{null.value}{Hypothesized null values for the diagonal entries \eqn{(1,1)} and \eqn{(2,2)} in the RCT, 
#' which are \eqn{E(N_{11})=R P_{aa}} and \eqn{E(N_{22})=R P_{ab}}, respectively).}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{rfct}, returned by \code{Xsq.nnref.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.nnref} only}
#'  
#' @seealso \code{\link{Znnref.ct}}, \code{\link{Znnref}}, \code{\link{Zself.ref.ct}},
#' \code{\link{Zself.ref}}, \code{\link{Zmixed.nonref.ct}} and \code{\link{Zmixed.nonref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsXsq.nnref
NULL
#'
#' @rdname funsXsq.nnref
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' R<-Rval(W)
#' Tv<-Tval(W,R)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Xsq.nnref(Y,cls)
#' Xsq.nnref.ct(rfct,nvec,Qv,Tv)
#'
#' Xsq.nnref(Y,cls,method="max")
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#'
#' ipd<-ipd.mat(Y)
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' R<-Rval(W)
#' Tv<-Tval(W,R)
#'
#' nvec<-as.numeric(table(cls))
#' rfct<-rct(ipd,cls)
#'
#' Xsq.nnref(Y,cls)
#' Xsq.nnref.ct(rfct,nvec,Qv,Tv)
#'
#' @export
Xsq.nnref.ct <- function(rfct,nvec,Qv,Tv)
{
  n<-sum(nvec)
  k<-length(nvec)
  
  if (k<4) 
  {nvec<-c(nvec,rep(0,4-k))} 
  
  Nsr<-rfct[1,1]
  Nmn<-rfct[2,2]
  
  R<-sum(rfct[1,]) #row sum for rfct
  Paa<-P11(nvec)
  ENsr<-R*Paa
  Paaaa<-P1111(nvec)
  Paabb<-P1122(nvec)
  VNsr<-R^2*(Paaaa+Paabb-Paa^2)+2*R*(Paa-Paaaa-Paabb)
  
  Pab<-P12(nvec)
  ENmn<-(n-R)*Pab
  Paabc<-P1123(nvec)
  Pabcd<-P1234(nvec)
  Paab<-P112(nvec)
  Pabc<-P123(nvec)
  VNmn<-(n-R)^2*(2*Paabb+4*Paabc+Pabcd-Pab^2)+(n-R)*Pab+(2*n-2*R+Qv-4*Tv)*(Paab+Pabc)+
   (-3*n+3*R-Qv+4*Tv)*(2*Paabb+4*Paabc+Pabcd)
  
  Paaab<-P1112(nvec)
  CovNsrNmn<-R*(n-R)*(2*Paaab+Paabc-Paa*Pab)+2*T*(Paab-2*Paaab-Paabc)
  ref.vec<-diag(rfct)
  covar<- matrix(c(VNsr,CovNsrNmn,CovNsrNmn,VNmn),nrow=2,byrow="T")
  ts<- t(ref.vec-c(ENsr,ENmn)) %*% solve(covar) %*% (ref.vec-c(ENsr,ENmn))
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the chi-squared test of NN self reflexivity is not defined')}
  
  nu<-2
  pval<- 1-pchisq(ts,df=nu)
  
  method <-"chi-squared Test of NN Self-Reflexivity and Mixed Non-Reflexivity"
  
  dname <-deparse(substitute(rfct))
  
  estimate<-ref.vec
  null.val<-c(ENsr,ENmn)
  names(estimate)<-names(null.val)<-c("self-reflexive","mixed non-reflexive")
  
  estimate.name <-c("diagonal entries of rct")
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value =null.val,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsXsq.nnref
#'
#' @export
Xsq.nnref <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  R<-Rval(W)
  Tv<-Tval(W,R)
  
  nvec<-as.numeric(table(lab))
  rfct<-rct(ipd,lab,...)
  
  rval<-Xsq.nnref.ct(rfct,nvec,Qv,Tv)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funs.scct
#'
#' @title Species Correspondence Contingency Table (SCCT)
#'
#' @description
#' Two functions: \code{scct.ct} and \code{scct}.
#' 
#' Both functions return the \eqn{k \times 2} species correspondence contingency table (SCCT) 
#' but have different arguments (see the parameter list below).
#' 
#' SCCT is constructed by categorizing the NN pairs according to pair type as self or mixed. 
#' A base-NN pair is called a self pair, if the elements of the pair are from the same class;
#' a mixed pair, if the elements of the pair are from different classes.
#' Row labels in the RCT are the class labels and the column labels are \code{"self"} and \code{"mixed"}.
#' The \eqn{k \times 2} SCCT (whose first column is self column with entries \eqn{S_i} and second column is mixed with entries \eqn{M_i})
#' is closely related to the \eqn{k \times k} nearest neighbor contingency table (NNCT) whose entries are \eqn{N_{ij}},
#' where \eqn{S_i=N_{ii}} and \eqn{M_i=n_i-N_{ii}} with \eqn{n_i} is the size of class \eqn{i}.
#' 
#' The function \code{scct.ct} returns the SCCT given the inter-point distance (IPD) matrix or data set \code{x},
#' and the  function \code{scct} returns the SCCT given the IPD matrix. SCCT is a \eqn{k \times 2} matrix where \eqn{k} is 
#' number of classes in the data set.
#' (See \insertCite{ceyhan:NNCorrespond2018;textual}{nnspat} for more detail,
#' where SCCT is labeled as CCT for correspondence contingency table).
#' 
#' The argument \code{ties} is a logical argument (default=\code{FALSE} for both functions) to take ties into account or not.
#' If \code{TRUE} a NN contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' 
#' The argument nnct is a logical argument for \code{scct.ct} only (default=\code{FALSE}) to determine the structure of the
#' argument \code{x}. If \code{TRUE}, \code{x} is taken to be the \eqn{k \times k} NNCT, and if \code{FALSE}, \code{x} is taken to be the IPD matrix.
#' 
#' The argument lab is the \code{vector} of class labels (default=\code{NULL} when \code{nnct=TRUE} in the function \code{scct.ct} and no default
#' specified for scct).
#' 
#' @param x The IPD matrix (if \code{nnct=FALSE}) or the NNCT (if \code{nnct=TRUE}), used in \code{scct.ct} only
#' @param lab The \code{vector} of class labels (numerical or categorical), default=\code{NULL} when \code{nnct=FALSE} in the function \code{scct.ct} and no default
#' specified for scct.
#' @param ties A logical argument (default=\code{FALSE}) to take ties into account or not. If \code{TRUE} a NN 
#' contributes \eqn{1/m} to the NN count if it is one of the \eqn{m} tied NNs of a subject.
#' @param nnct A logical parameter (default=\code{FALSE}). If \code{TRUE}, \code{x} is taken to be the \eqn{k \times k} NNCT, 
#' and if \code{FALSE}, \code{x} is taken to be the IPD matrix, used in \code{scct.ct} only.
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{scct} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{scct} only
#' 
#' @return Returns the \eqn{k \times 2} SCCT where \eqn{k} is the number of classes in the data set.
#'  
#' @seealso \code{\link{nnct}}, \code{\link{tct}}, \code{\link{rct}} and \code{\link{Qsym.ct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.scct
NULL
#'
#' @rdname funs.scct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' NNCT<-nnct(ipd,cls)
#' NNCT
#'
#' scct(Y,cls)
#' scct(Y,cls,method="max")
#'
#' scct.ct(ipd,cls)
#' scct.ct(ipd,cls,ties = TRUE)
#' scct.ct(NNCT,nnct=TRUE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' scct.ct(ipd,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' NNCT<-nnct(ipd,cls)
#' NNCT
#'
#' scct(Y,cls)
#'
#' scct.ct(ipd,cls)
#' scct.ct(NNCT,nnct=TRUE)
#'
#' @export
scct.ct <- function(x,lab=NULL,ties=FALSE,nnct=FALSE)
{
  if (nnct==FALSE)
  {
    if (is.null(lab))
    {stop('lab must be provided as the list of class labels')}
    
    ipd<-x
    n<-nrow(ipd)
    
    ord<-order(lab)#ordering the class label lab first 
    #(to be consistent with row and column labeling in the NNCT)
    lab<-sort(lab) 
    ipd<-ipd[ord,ord]
    flab<-as.factor(lab) #converting class labels to factors
    clnames<-levels(flab)
    k<-length(clnames)
    
    ct<-matrix(0,k,2)  
    if (n<=1)
    { colnames(ct)<-c("self","mixed")
    return(ct)}
    
    nlab<-as.numeric(flab)  #converting class labels to numbers
    for(i in 1:n)
    {
      ind <- NN(ipd,i);
      lind<-length(ind)
      for (j in 1:lind)
      {
        addend<-ifelse(ties==FALSE,1,1/lind)
        ifelse(nlab[i]==nlab[ind[j]],ct[nlab[i],1]<- ct[nlab[i],1]  + addend,ct[nlab[i],2]<- ct[nlab[i],2]  + addend)
      }
    }
    rownames(ct)<-clnames #row names for the SCCT
  } else
  { k<-nrow(x)
  ct<-matrix(0,k,2)  
  ct[,1]<-diag(x)
  rs<-row.sum(x)
  ct[,2]<-rs-ct[,1]
  rownames(ct)<-rownames(x)
  }
  colnames(ct)<-c("self","mixed")
  ct
} #end for the function
#'
#' @rdname funs.scct
#'
#' @export
scct <- function(dat,lab,ties=FALSE,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-scct.ct(ipd,lab,ties,nnct=FALSE)
  ct
} #end for the function
#'

#################################################################

# funs.varNii
#'
#' @title Variances of the Self Entries in a Species Correspondence Contingency Table (SCCT)
#'
#' @description
#' Two functions: \code{varNii.ct} and \code{varNii}.
#' 
#' Both functions return a \code{vector} of length \eqn{k} of variances of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the variances of the diagonal entries \eqn{N_{ii}} in an NNCT,
#' but have different arguments (see the parameter list below).
#' These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' The argument \code{ct} which is used in \code{varNii.ct} only, can be either the NNCT or SCCT.
#' 
#' See also (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat}).
#' 
#' @param ct The NNCT or SCCT, used in \code{varNii.ct} only
#' @param Q The number of shared NNs, used in \code{varNii.ct} only
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs), used in \code{varNii.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{varNii} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{varNii} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{varNii} only
#' 
#' @return A \code{vector} of length \eqn{k} whose entries are the variances of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or of the diagonal entries in an NNCT.
#'  
#' @seealso \code{\link{scct}}, \code{\link{var.nnct}}, \code{\link{var.tct}}, \code{\link{var.nnsym}}
#' and \code{\link{covNii}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.varNii
NULL
#'
#' @rdname funs.varNii
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' varNii(Y,cls)
#' varNii.ct(ct,Qv,Rv)
#'
#' varNii(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' varNii(Y,fcls)
#' varNii.ct(ct,Qv,Rv)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' varNii(Y,cls)
#' varNii.ct(ct,Qv,Rv)
#'
#' @export
varNii.ct <- function(ct,Q,R) 
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct) 
  
  var<-vector();
  for (i in 1:k)
  {
    Paa <- p11(rs[i],n); Paaa <- p111(rs[i],n); Paaaa <- p1111(rs[i],n)
    var<- c(var,(n+R)*Paa+(2*n-2*R+Q)*Paaa+(n^2-3*n-Q+R)*Paaaa-(n*Paa)^2)
  }
  var
} #end for the function
#'
#' @rdname funs.varNii
#'
#' @export
varNii <- function(dat,lab,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  var<-varNii.ct(ct,Qv,Rv) 
var
} #end for the function
#'

#################################################################

# funs.covNii
#'
#' @title Covariance Matrix of the Self Entries in a Species Correspondence Contingency Table (SCCT)
#'
#' @description
#' Two functions: \code{covNii.ct} and \code{covNii}.
#' 
#' Both functions return the covariance matrix of the self entries (i.e. first column entries) in a
#' species correspondence contingency table (SCCT)
#' but have different arguments (see the parameter list below).
#' The covariance matrix is of dimension \eqn{k \times k} and its entries are \eqn{cov(S_i,S_j)} where \eqn{S_i} values are
#' the entries in the first column of SCCT (recall that \eqn{S_i} equals diagonal entry \eqn{N_{ii}} in the NNCT).
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' The argument \code{ct} which is used in \code{covNii.ct} only, can be either the NNCT or SCCT.
#' And the argument \code{Vsq} is the vector of variances of the diagonal entries \eqn{N_{ii}} in the NNCT or the self entries
#' (i.e. the first column) in the SCCT.
#' 
#' See also (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat}).
#' 
#' @param ct The NNCT or SCCT, used in \code{covNii.ct} only
#' @param Vsq The \code{vector} of variances of the diagonal entries \eqn{N_{ii}} in the NNCT or the self entries
#' (i.e. the first column) in the SCCT, used in \code{covNii.ct} only
#' @param Q The number of shared NNs, used in \code{covNii.ct} only
#' @param R The number of reflexive NNs (i.e., twice the number of reflexive NN pairs), used in \code{covNii.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{covNii} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{covNii} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function, used in \code{covNii} only
#' 
#' @return A \code{vector} of length \eqn{k} whose entries are the variances of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT).
#' 
#' @return The \eqn{k \times k} covariance matrix of cell counts \eqn{S_i} in the self (i.e., first) column of the SCCT
#' or of the diagonal cell counts \eqn{N_{ii}} for \eqn{i=1,\ldots,k} in the NNCT.
#'  
#' @seealso \code{\link{scct}}, \code{\link{cov.nnct}}, \code{\link{cov.tct}} and \code{\link{cov.nnsym}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.covNii
NULL
#'
#' @rdname funs.covNii
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' covNii(Y,cls)
#' covNii.ct(ct,vsq,Qv,Rv)
#'
#' covNii(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' covNii(Y,fcls)
#' covNii.ct(ct,vsq,Qv,Rv)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' covNii(Y,cls)
#' covNii.ct(ct,vsq,Qv,Rv)
#'
#' @export
covNii.ct <- function(ct,Vsq,Q,R)
{
  rs <- row.sum(ct); 
  n<-sum(ct)
  m<-nrow(ct); 
  cov.diag<-matrix(0,m,m);
  
  for (i in 1:m)
  { 
    for (j in i:m)
      if (i==j)
        cov.diag[i,j] <- Vsq[i]
      else 
      {Paa<- p11(rs[i],n); Pbb<- p11(rs[j],n); Paabb<- p1122(rs[i],rs[j],n)
      cov.diag[i,j] <- (n^2-3*n-Q+R)*Paabb-n^2*Paa*Pbb
      cov.diag[j,i] <- cov.diag[i,j]
      }
  }
  cov.diag
} #end for the function 
#'
#' @rdname funs.covNii
#'
#' @export
covNii <- function(dat,lab,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Q<-Qvec(W)$q
  R<-Rval(W)
  
  Vsq<-varNii.ct(ct,Q,R) 
  
  rs <- row.sum(ct); 
  n<-sum(ct)
  m<-nrow(ct); 
  cov.diag<-matrix(0,m,m);
  
  for (i in 1:m)
  { 
    for (j in i:m)
      if (i==j)
        cov.diag[i,j] <- Vsq[i]
      else 
      {Paa<- p11(rs[i],n); Pbb<- p11(rs[j],n); Paabb<- p1122(rs[i],rs[j],n)
      cov.diag[i,j] <- (n^2-3*n-Q+R)*Paabb-n^2*Paa*Pbb
      cov.diag[j,i] <- cov.diag[i,j]
      }
  }
  cov.diag
} #end for the function 
#'

#################################################################

#' @title Expected Values of the Self Entries in a Species Correspondence Contingency Table (SCCT)
#'
#' @description Returns a \code{vector} of length \eqn{k} of expected values of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the expected values of the diagonal entries \eqn{N_{ii}} in an NNCT.
#' These expected values are valid under RL or CSR.
#' 
#' The argument \code{ct} can be either the NNCT or SCCT.
#' 
#' See also (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat}).
#' 
#' @param ct The NNCT or SCCT
#' 
#' @return A \code{vector} of length \eqn{k} whose entries are the expected values of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or of the diagonal entries in an NNCT.
#'  
#' @seealso \code{\link{scct}} and \code{\link{EV.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' EV.Nii(ct)
#' ct<-scct(ipd,cls)
#' EV.Nii(ct)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' EV.Nii(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' EV.Nii(ct)
#' ct<-scct(ipd,cls)
#' EV.Nii(ct)
#'
#' @export
EV.Nii <- function(ct)
{
  nvec<-row.sum(ct)
  k<-length(nvec);
  n<-sum(nvec) 
  
  EN<- vector();
  for (i in 1:k)
  {
    EN<- c(EN,nvec[i]*(nvec[i]-1)/(n-1))  
  }
  
  EN
} #end for the function
#'

#################################################################

# funsZnnself
#'
#' @title Self-Reflexivity Tests with Normal Approximation
#'
#' @description
#' Two functions: \code{Znnself.ct} and \code{Znnself}.
#' 
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the self entries (i.e. first column)
#' in a species correspondence contingency table (SCCT) or the expected values of the diagonal entries \eqn{N_{ii}} in
#' an NNCT to the ones under RL or CSR.
#' That is, each performs NN self reflexivity for each class test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' NN self reflexivity is for each class can be viewed as a decomposition of species correspondence for
#' each class.
#' (See \insertCite{ceyhan:NNCorrespond2018;textual}{nnspat} for more detail).
#' 
#' Each test is based on the normal approximation of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the diagonal entries \eqn{N_{ii}} in an NNCT and
#' are due to \insertCite{ceyhan:NNCorrespond2018}{nnspat}.
#' 
#' Each function yields a \code{vector} of length \eqn{k} of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, null values (i.e. expected values), sample estimates (i.e. observed values) of self entries 
#' in the SCCT or diagonal entries in the NNCT, a \eqn{k \times 2} matrix of confidence intervals (where each row is the
#' confidence interval for self entry \eqn{S_i} in the SCCT or diagonal entry \eqn{N_{ii}} in the NNCT) and
#' also names of the test statistics, estimates, null values and the method and the data
#' set used.
#' 
#' The null hypothesis is that all \eqn{E[S_i] = E[N_{ii}] = n_i(n_i - 1)/(n - 1)} where \eqn{n_i} is the size of class \eqn{i} and
#' \eqn{n} is the data size.
#' 
#' The \code{Znnself} functions (i.e. \code{Znnself.ct} and \code{Znnself}) are different from the \code{Znnref} functions 
#' (i.e. \code{\link{Znnref.ct}} and \code{\link{Znnref}}) and from \code{Zself.ref} functions (i.e. \code{\link{Zself.ref.ct}} and \code{\link{Zself.ref}}) and also from
#' \code{Znnself.sum} functions (i.e. \code{Znnself.sum.ct} and \code{Znnself.sum}).
#' \code{Znnself} functions are testing the self reflexivity at a class-specific level (i.e. for each class) using the
#' first column in the SCCT, while \code{Zself.ref} functions are for testing the self reflexivity for the entire data set
#' using entry \eqn{(1,1)} in RCT, and \code{Znnref} functions are for testing the self reflexivity and mixed non-reflexivity
#' using the diagonal entries in the RCT, and
#' \code{Znnself.sum} functions are testing the cumulative species correspondence using the sum of the self column (i.e.,
#' the first column) in the SCCT.
#' 
#' @param ct The NNCT or SCCT, used in \code{Znnself.ct} only 
#' @param VarNii The variance vector of differences of self entries in the SCCT or diagonal entries in the NNCT,
#' used in \code{Znnself.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnself} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnself} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the self entries in the SCCT or diagonal entries in the NNCT
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnself} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{vector} (of length k) of \eqn{Z} test statistics for NN self reflexivity test}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{vector} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Lower and Upper Confidence Levels, it is \code{NULL} here since we provide confidence intervals
#' as a \eqn{k \times 2} matrix.} 
#' \item{conf.int}{The \eqn{k \times 2} matrix of confidence intervals for the estimates, (where each row is the
#' confidence interval for self entry \eqn{S_i} in the SCCT or diagonal entry \eqn{N_{ii}} in the NNCT).}
#' \item{cnf.lvl}{Level of the confidence intervals (i.e., conf.level) for the self entries in the SCCT or
#' diagonal entries in the NNCT.}
#' \item{estimate}{The \code{vector} of estimates of the parameters, i.e., observed values of self entries 
#' in the SCCT or diagonal entries in the NNCT.}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function.}
#' \item{null.value}{The \code{vector} of null values of the parameters, i.e., expected values of self entries 
#' in the SCCT or diagonal entries in the NNCT under RL or CSR.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Znnself.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Znnself} only}
#'  
#' @seealso \code{\link{Zself.ref.ct}}, \code{\link{Zself.ref}}, \code{\link{Znnref.ct}},
#' \code{\link{Znnref}}, \code{\link{Xsq.spec.cor}} and \code{\link{Xsq.spec.cor.ct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnself
NULL
#'
#' @rdname funsZnnself
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' VarN.diag<-varNii.ct(ct,Qv,Rv)
#'
#' Znnself(Y,cls)
#' Znnself(Y,cls,alt="g")
#'
#' Znnself.ct(ct,VarN.diag)
#' Znnself.ct(ct,VarN.diag,alt="g")
#'
#' Znnself(Y,cls,method="max")
#'
#' ct<-scct(ipd,cls)
#' Znnself.ct(ct,VarN.diag)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Znnself(Y,fcls)
#' Znnself.ct(ct,VarN.diag)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' VarN.diag<-varNii.ct(ct,Qv,Rv)
#'
#' Znnself(Y,cls,alt="l")
#' Znnself.ct(ct,VarN.diag)
#' Znnself.ct(ct,VarN.diag,alt="l")
#'
#' @export
Znnself.ct <- function(ct,VarNii,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-diag(ct) #diagonal vector of NNCT
  estimate.name <-"diagonal entries of NNCT (i.e. number of self NNs)"
  
  # rs<-row.sum(ct)
  ENii<-EV.Nii(ct)
  null.val<-ENii
  names.null <-"diagonal entries of NNCT"
  
  k<-nrow(ct);
  Nii<-diag(ct)
  
  stderr<-sqrt(VarNii)
  ts<-(Nii-ENii)/stderr
  
  if (all(is.na(ts)))
  {stop('Both test statistics are NaN, so the self NN Z-tests are not defined')}
  
  cint<-matrix(0,ncol=2,nrow=k)
  alt<- switch(alternative,
         less = { 
           for (i in 1:k)
           {cint[i,] <-estimate[i]+c(-Inf, qnorm(conf.level))*stderr[i]}
           pval <-pnorm(ts)
        },
         greater = { 
           for (i in 1:k)
           {cint[i,] <-estimate[i]+c(-qnorm(conf.level),Inf)*stderr[i]}
           pval <-pnorm(ts, lower.tail = FALSE)
         },
         two.sided = { 
           alpha <-1 - conf.level
           Cint <-qnorm(1 - alpha/2)
           for (i in 1:k)
           {cint[i,] <-estimate[i]+c(-Cint, Cint)*stderr[i]}
           pval <-2 * pnorm(-abs(ts))
        }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("z-Tests for Diagonal NNCT Entries (for Self NN Pairs)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  names(null.val)<-clnames #names for the diagonal NNC\eqn{T} values
  rownames(cint)<-clnames
  colnames(cint)<-c("LCL","UCL")
  ts.names <-"Self NN Z-tests (for diagonal entries of NNCT)"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = NULL,UCL = NULL,
    conf.int = cint,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name, #this is for other functions to have a different description for the sample estimates
    null.value = null.val,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZnnself
#'
#' @export
Znnself <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  VarN.diag<-varNii.ct(ct,Qv,Rv) 
  
  rval<-Znnself.ct(ct,VarN.diag,alternative=alternative,conf.level=conf.level)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsXsq.spec.cor
#'
#' @title Overall Species Correspondence Test with Chi-square Approximation
#'
#' @description
#' Two functions: \code{Xsq.spec.cor.ct} and \code{Xsq.spec.cor}.
#' 
#' Each one performs hypothesis tests of (simultaneous) equality of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the expected values of the diagonal entries \eqn{N_{ii}} in an NNCT
#' to the ones under RL or CSR.
#' That is, each performs the overall species correspondence test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:NNCorrespond2018;textual}{nnspat} for more detail).
#' 
#' Each test is based on the Chi-square approximation of the corresponding quadratic form for the first column
#' in a species correspondence contingency table (SCCT) or the diagonal entries \eqn{N_{ii}} in an NNCT and
#' are due to \insertCite{ceyhan:NNCorrespond2018}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of the self entries (i.e. first column) in the SCCT
#' or the diagonal entries in the NNCT and also the sample estimates (i.e. observed values) of these entries.
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that all 
#' \eqn{E[S_1,S_2,\ldots,S_k] = E[N_{11},N_{22},\ldots,N_{kk}] = ((n_1(n_1 - 1)/(n - 1),(n_2(n_2 - 1)/(n - 1),\ldots,(n_k(n_k - 1)/(n - 1) )}
#' where \eqn{n_i} is the size of class \eqn{i} and \eqn{n} is the data size.
#' 
#' @param ct The NNCT or SCCT, used in \code{Xsq.spec.cor.ct} only
#' @param covSC The covariance matrix for the self entries (i.e. first column) in the SCCT
#' or the diagonal entries in the NNCT, used in \code{Xsq.spec.cor.ct} only. Usually output of the functions 
#' \code{\link{covNii.ct}} or \code{\link{covNii}}.
#' @param nnct A logical parameter (default=\code{FALSE}). If \code{TRUE}, \code{x} is taken to be the \eqn{k \times k} NNCT, 
#' and if \code{FALSE}, \code{x} is taken to be the IPD matrix, used in \code{Xsq.spec.cor.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Xsq.spec.cor} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.spec.cor} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Xsq.spec.cor} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for overall species correspondence test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k} for this function.}
#' \item{estimate}{The \code{vector} of estimates of the parameters, i.e., observed values of self entries 
#' in the SCCT or diagonal entries in the NNCT.}
#' \item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
#' \item{null.value}{The \code{vector} of null values of the parameters, i.e., expected values of self entries 
#' in the SCCT or diagonal entries in the NNCT under RL or CSR.}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Xsq.spec.cor.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.spec.cor} only}
#'  
#' @seealso \code{\link{Zself.ref.ct}}, \code{\link{Zself.ref}}, \code{\link{Xsq.nnref.ct}} and \code{\link{Xsq.nnref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsXsq.spec.cor
NULL
#'
#' @rdname funsXsq.spec.cor
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-scct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' cv<-covNii.ct(ct,vsq,Qv,Rv)
#' Xsq.spec.cor.ct(ct,cv)
#' Xsq.spec.cor(Y,cls)
#' Xsq.spec.cor(Y,cls,method="max")
#'
#' ct<-nnct(ipd,cls)
#' Xsq.spec.cor.ct(ct,cv,nnct = TRUE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-scct(ipd,fcls)
#' Xsq.spec.cor.ct(ct,cv)
#' Xsq.spec.cor(Y,fcls)
#'
#' ct<-nnct(ipd,fcls)
#' Xsq.spec.cor.ct(ct,cv,nnct=TRUE)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-scct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' cv<-covNii.ct(ct,vsq,Qv,Rv)
#' Xsq.spec.cor.ct(ct,cv)
#'
#' ct<-nnct(ipd,cls)
#' Xsq.spec.cor.ct(ct,cv,nnct = TRUE)
#' Xsq.spec.cor(Y,cls)
#'
#' @export
Xsq.spec.cor.ct <- function(ct,covSC,nnct=FALSE)
{
  if (nnct==FALSE)
  {
    nsq<-ct[,1]
    estimate.name <-c("self NN pairs (i.e. first column) in SCCT")
  } else
  {
    nsq<-diag(ct)
    estimate.name <-c("self NN pairs (i.e. diagonal entries) in NNCT")
  }
  Ensq<-EV.Nii(ct)
  
  k<-length(nsq)
  ts<- t(nsq-Ensq) %*% ginv(covSC,tol=1.490116e-20) %*% (nsq-Ensq)
  nu<-k
  pval<- 1-pchisq(ts,df=nu)
  
  method <-"chi-squared Test of (Overall) Species Correspondence in NN Structure"
  
  dname <-deparse(substitute(ct))
  
  estimate<-nsq
  null.val<-Ensq
  names(estimate)<-names(null.val)<-row.names(ct)
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value =null.val,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsXsq.spec.cor
#'
#' @export
Xsq.spec.cor <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-scct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  
  vsq<-varNii.ct(ct,Qv,Rv) 
  cv<-covNii.ct(ct,vsq,Qv,Rv) 
  rval<-Xsq.spec.cor.ct(ct,cv)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZnnself.sum
#'
#' @title Cumulative Species Correspondence Test with Normal Approximation
#'
#' @description
#' Two functions: \code{Znnself.sum.ct} and \code{Znnself.sum}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected value of the sum of the self entries (i.e.
#' first column) in a species correspondence contingency table (SCCT) or the expected values of the sum of the 
#' diagonal entries \eqn{N_{ii}} in an NNCT to the one under RL or CSR.
#' That is, each performs a cumulative species correspondence test which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:NNCorrespond2018;textual}{nnspat} for more detail).
#' 
#' Each test is based on the normal approximation of the sum of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the sum of the diagonal entries \eqn{N_{ii}} in an NNCT and
#' are due to \insertCite{ceyhan:NNCorrespond2018}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, sample estimate (i.e. observed value) and null (i.e., expected) value for the
#' sum of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the sum of the diagonal entries \eqn{N_{ii}} in an NNCT, 
#' and method and name of the data set used.
#' 
#' The null hypothesis is that all 
#' \eqn{E[S] = \sum_{i=1}^k n_i(n_i - 1)/(n - 1)} where \eqn{S} is the sum of the self column
#' in the SCCT, \eqn{n_i} is the size of class \eqn{i} and \eqn{n} is the data size. 
#' 
#' The \code{Znnself.sum} functions (i.e. \code{Znnself.sum.ct} and \code{Znnself.sum}) are different from the Znnself
#' functions (i.e. \code{Znnself.ct} and \code{Znnself}), and from the \code{Znnref} functions 
#' (i.e. \code{\link{Znnref.ct}} and \code{\link{Znnref}}) and also from \code{Zself.ref} functions (i.e. \code{\link{Zself.ref.ct}} and \code{\link{Zself.ref}}).
#' \code{Znnself.sum} functions are testing the cumulative species correspondence using the sum of the self column (i.e.,
#' the first column) in the SCCT, while \code{Znnself} functions are testing the self reflexivity at a class-specific level (i.e. for each class) using the
#' first column in the SCCT, while \code{Zself.ref} functions are for testing the self reflexivity for the entire data set
#' using entry \eqn{(1,1)} in RCT, and \code{Znnref} functions are for testing the self reflexivity and mixed non-reflexivity
#' using the diagonal entries in the RCT.
#' 
#' @param ct The NNCT or SCCT, used in \code{Znnself.sum.ct} only 
#' @param covSC The covariance matrix for the self entries (i.e. first column) in the SCCT
#' or the diagonal entries in the NNCT, used in \code{Znnself.sum.ct} only. Usually output of the functions 
#' \code{\link{covNii.ct}} or \code{\link{covNii}}.
#' @param nnct A logical parameter (default=\code{FALSE}). If \code{TRUE}, \code{x} is taken to be the \eqn{k \times k} NNCT, 
#' and if \code{FALSE}, \code{x} is taken to be the IPD matrix, used in \code{Znnself.sum.ct} only
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Znnself.sum} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Znnself.sum} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the self entries in the SCCT or diagonal entries in the NNCT
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Znnself.sum} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the overall species correspondence test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the sum of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the sum of the diagonal entries \eqn{N_{ii}} in an NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the observed sum of the self entries (i.e. first column)
#' in a species correspondence contingency table (SCCT) or the sum of the diagonal entries \eqn{N_{ii}} in an NNCT.}
#' \item{null.value}{Hypothesized null value for the sum of the self entries (i.e. first column) in a
#' species correspondence contingency table (SCCT) or the sum of the diagonal entries \eqn{N_{ii}} in an NNCT
#' which is \eqn{E[S] = \sum_{i=1}^k n_i(n_i - 1)/(n - 1)} where \eqn{S} is the sum of the self column
#' in the SCCT, \eqn{n_i} is the size of class \eqn{i} and \eqn{n} is the data size.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Znnself.sum.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Znnself.sum} only}
#'  
#' @seealso \code{\link{Znnself.ct}}, \code{\link{Znnself}}, \code{\link{Znnref.ct}}, \code{\link{Znnref}},
#' \code{\link{Zself.ref.ct}} and \code{\link{Zself.ref}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZnnself.sum
NULL
#'
#' @rdname funsZnnself.sum
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-scct(ipd,cls)
#' ct
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' cv<-covNii.ct(ct,vsq,Qv,Rv)
#'
#' Znnself.sum(Y,cls)
#'
#' Znnself.sum.ct(ct,cv)
#' Znnself.sum.ct(ct,cv,alt="g")
#'
#' Znnself.sum(Y,cls,method="max")
#'
#' ct<-nnct(ipd,cls)
#' Znnself.sum.ct(ct,cv,nnct = TRUE)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-scct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#'
#' vsq<-varNii.ct(ct,Qv,Rv)
#' cv<-covNii.ct(ct,vsq,Qv,Rv)
#'
#' Znnself.sum(Y,cls)
#'
#' Znnself.sum.ct(ct,cv)
#' Znnself.sum.ct(ct,cv,alt="g")
#'
#' ct<-nnct(ipd,cls)
#' Znnself.sum.ct(ct,cv,nnct = TRUE)
#'
#' Znnself.sum(Y,cls,alt="g")
#'
#' @export
Znnself.sum.ct <- function(ct,covSC,nnct=FALSE,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  ifelse(nnct==FALSE, nsq<-ct[,1], nsq<-diag(ct))
  
  estimate<-sum.dg<-sum(nsq)
  
  ifelse(nnct==FALSE, names(estimate) <-"sum of self NN pairs (i.e. sum of first column in SCCT)",
         names(estimate) <-"sum of self NN pairs (i.e. diagonal entries in NNCT)")
  
  Ensq<-EV.Nii(ct)
  Esum.dg<-sum(Ensq)
  var.sum.dg<-sum(covSC)
  stderr<-sqrt(var.sum.dg)
  
  ts<-(sum.dg-Esum.dg)/stderr
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the Z-test for sum of self NN pairs is not defined')}
  
  null.val<-Esum.dg
  names(null.val)<- "(expected) sum of self NN pairs"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  method <-"Z-test for Sum of Self NN Pairs"
  names(ts) <-"Cumulative Species-Correspondence Z-test"
  
  dname <-c("contingency table = ",deparse(substitute(ct)))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZnnself.sum
#'
#' @export
Znnself.sum <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-scct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  
  vsq<-varNii.ct(ct,Qv,Rv) 
  cv<-covNii.ct(ct,vsq,Qv,Rv) 
  
  rval<- Znnself.sum.ct(ct,cv,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

################################################

# funs.seg.coeff
#'
#' @title Pielou's Segregation Coefficients for NNCTs
#'
#' @description
#' Two functions: \code{Pseg.coeff} and \code{seg.coeff}.
#' 
#' Each function computes segregation coefficients based on NNCTs.
#'The function \code{Pseg.coeff} computes Pielou's segregation coefficient (\insertCite{pielou:1961;textual}{nnspat})
#' for the two-class case (i.e., based on \eqn{2 \times 2} NNCTs)
#' and \code{seg.coeff} is the extension of \code{Pseg.coeff} to the multi-class case (i.e. for \eqn{k \times k} NNCTs with \eqn{k \ge 2})
#' and provides a \eqn{k \times k} matrix of segregation coefficients
#' (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' Both functions use the same argument, \code{ct}, for NNCT.
#' 
#' Pielou's segregation coefficient (for two classes) is \eqn{S_P = 1-(N_{12} + N_{21})/(E[N_{12}] + E[N_{21}])}
#' and the extended segregation coefficents (for \eqn{k \ge 2} classes) are 
#' \eqn{S_c = 1 -(N_{ii})/(E[N_{ii}])} for the diagonal cells in the NNCT
#' and
#' \eqn{S_c = 1 -(N_{ij} + N_{ji})/(E[N_{ij}] + E[N_{ji}])} for the off-diagonal cells in the NNCT.
#'   
#' @param ct A nearest neighbor contingency table, used in both functions
#' 
#' @return
#' \code{Pseg.coeff} returns Pielou's segregation coefficient for \eqn{2 \times 2} NNCT
#' \code{seg.coeff} returns a \eqn{k \times k} matrix of segregation coefficients (which are extended versions
#' of Pielou's segregation coefficient)
#'   
#' @seealso \code{\link{seg.ind}}, \code{\link{Zseg.coeff.ct}} and \code{\link{Zseg.coeff}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.seg.coeff
NULL
#'
#' @rdname funs.seg.coeff
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for Pseg.coeff
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#' Pseg.coeff(ct)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Pseg.coeff(ct)
#'
#' #############
#' ct<-matrix(sample(1:25,9),ncol=3)
#' #Pseg.coeff(ct)
#'
#' @export
Pseg.coeff <- function(ct)
{
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  rs <- row.sum(ct)
  n<-sum(rs)
  
  SP<-1-(ct[1,2]+ct[2,1])/(2*rs[1]*rs[2]/(n-1))
  names(SP)<-NULL
  SP
} #end for the function
#'
#' @rdname funs.seg.coeff
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for seg.coeff
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#' seg.coeff(ct)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' seg.coeff(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-nnct(ipd,cls)
#'
#' seg.coeff(ct)
#'
#' @export
seg.coeff <- function(ct)
{
  rs <- row.sum(ct); 
  k<-nrow(ct);
  n<-sum(ct)
  
  if (n<=1)
  {stop('sample size, n, must be >=2 for the coefficient of segregation to be defined')}
  Scf<- matrix(0,k,k); 
  for (i in 1:k)
    for (j in i:k)
    {
      if (i == j)
        Scf[i,i]<- 1-(n-1)*ct[i,i]/(rs[i]*(rs[i]-1)) 
      else 
        Scf[i,j]<- 1-(n-1)*(ct[i,j]+ct[j,i])/(2*rs[i]*(rs[j])) 
      Scf[j,i]<-Scf[i,j]
    }
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(Scf)<-colnames(Scf)<-clnames #row and column names for the segregation coefficient matrix
  
  Scf
} #end for the function
#'

#################################################################

#' @title Variance of Pielou's Segregation Coefficient for 2 Classes
#'
#' @description Returns the variance of Pielou's coefficient of segregation for the two-class case
#' (i.e., based on \eqn{2 \times 2} NNCTs)in a \eqn{2 \times 2} NNCT. 
#' This variance is valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}) for more detail.
#'
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#'
#' @return The variance of Pielou's coefficient of segregation for the two-class case.
#'
#' @seealso \code{\link{Pseg.coeff}}, \code{\link{seg.coeff}} and \code{\link{var.seg.coeff}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' varPseg.coeff(ct,covN)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' varPseg.coeff(ct,covN)
#'
#' #############
#' ct<-matrix(sample(1:25,9),ncol=3)
#' #varPseg.coeff(ct,covN)
#'
#' @export
varPseg.coeff <- function(ct,covN) 
{
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  rs <- row.sum(ct); 
  n<-sum(ct) 
  
  Var.SP<- ((n-1)/(2*rs[1]*rs[2]))^2*(covN[2,2]+covN[3,3]+2*covN[2,3])
  
  names(Var.SP)<-NULL
  Var.SP
} #end for the function
#'

#################################################################

#' @title Index Matrix for Computing the Covariance of Segregation Coefficients
#'
#' @description Returns the index matrix for choosing the entries in the covariance matrix for NNCT 
#' used for computing the covariance for the extension of Pielou's segregation coefficient to the multi-class
#' case. The matrix is \eqn{k(k+1)/2 \times 2} with each row is the \eqn{i,j} 
#' corresponding to \eqn{N_{ij}} in the NNCT.
#'
#' @param k An integer specifying the number of classes in the data set
#' 
#' @return The \eqn{k(k+1)/2 \times 2} index matrix with each row is the \eqn{i,j} 
#' corresponding to \eqn{N_{ij}} in the NNCT
#' 
#' @seealso \code{\link{cov.seg.coeff}}, \code{\link{seg.coeff}} and \code{\link{ind.nnsym}} 
#'
#' @author Elvan Ceyhan
#'
ind.seg.coeff <- function(k)
{
  ind1<-vector()
  for (i in 1:k)
  {
    for (j in i:k)
    {
      ind1<-c(ind1,c(i,j))
    }
  }
  ind.mat1<-matrix(ind1,ncol=2,byrow=T) #indices for the \code{T} vector
  ind.mat1
} #end for the function
#'

################################################

#' @title Variances of Segregation Coefficients in a Multi-class Case
#'
#' @description Returns the variances of segregation coefficients in a multi-class case based on the NNCT, \code{ct} 
#' in a \code{vector} of length \eqn{k(k+1)/2}, the order of the variances are as in the order of rows output of 
#' \code{\link{ind.seg.coeff}(k)}. These variances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' 
#' The argument \code{covN} is the covariance matrix of \eqn{N_{ij}} (concatenated rowwise).
#'
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#'
#' @return A \code{vector} of length \eqn{k(k+1)/2}, whose entries are the variances of segregation coefficients for the
#' entry \eqn{i,j} in the NNCT, where the order of the variances are as in the order of rows output of 
#' \code{\link{ind.seg.coeff}(k)}.
#'
#' @seealso \code{\link{seg.coeff}}, \code{\link{cov.seg.coeff}}, \code{\link{var.nnsym}}
#' and \code{\link{var.nnct}} and 
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' var.seg.coeff(ct,covN)
#' varPseg.coeff(ct,covN)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' var.seg.coeff(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' var.seg.coeff(ct,covN)
#'
#' @export
var.seg.coeff <- function(ct,covN)
{
  k<-nrow(ct)
  ent2<-cbind(rep(1:k,rep(k,k)),rep(1:k,k))
  ind.mat1<-ind.seg.coeff(k)
  
  EN<-EV.nnct(ct)
  
  kl<-k*(k+1)/2
  varT<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]
    c1<-rc1[2]
    ri1<-which(ent2[,1]==r1 & ent2[,2]==c1) #row indices for the N vector
    ri2<-which(ent2[,1]==c1 & ent2[,2]==r1)
    
    varT<-c(varT, (covN[ri1,ri1]+covN[ri1,ri2]+covN[ri2,ri1]+covN[ri2,ri2])/((2*EN[r1,c1])^2))
  }
  
  varT
} #end for the function
#'

#################################################################

#' @title Covariance Matrix of Segregation Coefficients in a Multi-class Case
#'
#' @description Returns the covariance matrix of the segregation coefficients in a multi-class case based on
#' the NNCT, \code{ct}. The covariance matrix is of dimension \eqn{k(k+1)/2 \times k(k+1)/2} and its entry \eqn{i,j} correspond to the
#' entries in the rows \eqn{i} and \eqn{j} of the output of \code{\link{ind.seg.coeff}(k)}. 
#' The segregation coefficients in the multi-class case are the extension of Pielou's segregation coefficient
#' for the two-class case.
#' These covariances are valid under RL or conditional on \eqn{Q} and \eqn{R} under CSR.
#' 
#' The argument \code{covN} is the covariance matrix of \eqn{N_{ij}} (concatenated rowwise).
#'
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#'
#' @param ct A nearest neighbor contingency table
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#'
#' @return The \eqn{k(k+1)/2} x \eqn{k(k+1)/2} covariance matrix of the segregation coefficients for the multi-class case
#' based on the NNCT, \code{ct} 
#'
#' @seealso \code{\link{seg.coeff}}, \code{\link{var.seg.coeff}}, \code{\link{cov.nnct}}
#' and \code{\link{cov.nnsym}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.seg.coeff(ct,covN)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' cov.seg.coeff(ct,covN)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' cov.seg.coeff(ct,covN)
#'
#' @export
cov.seg.coeff <- function(ct,covN)
{
  k<-sqrt(nrow(covN))
  ent2<-cbind(rep(1:k,rep(k,k)),rep(1:k,k))
  ind.mat1<-ind.seg.coeff(k)
  
  EN<-EV.nnct(ct)
  
  kl<-k*(k+1)/2
  varcovT<-matrix(0,nrow=kl,ncol=kl)
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]
    c1<-rc1[2]
    ri1<-which(ent2[,1]==r1 & ent2[,2]==c1) #row indices for the N vector
    ri2<-which(ent2[,1]==c1 & ent2[,2]==r1)
    for (j in i:kl)
    {
      rc1<-ind.mat1[j,]
      r1n<-rc1[1]
      c1n<-rc1[2]
      ci1<-which(ent2[,1]==r1n & ent2[,2]==c1n) #column indices for the N vector
      ci2<-which(ent2[,1]==c1n & ent2[,2]==r1n) 
      
      varcovT[i,j]<-(covN[ri1,ci1]+covN[ri1,ci2]+covN[ri2,ci1]+covN[ri2,ci2])/(2*(EN[r1,c1]+EN[r1n,c1n]))
      varcovT[j,i]<-varcovT[i,j]
    }
  }
  
  varcovT
} #end for the function
#'

#################################################################

# funsZseg.coeff
#'
#' @title Z Tests for Segregation Coefficients
#'
#' @description
#' Two functions: \code{Zseg.coeff.ct} and \code{Zseg.coeff}.
#' 
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' segregation coefficients from their expected values under RL or CSR for each segregation coefficient
#' in the NNCT.
#' 
#' The test for each cell \eqn{i,j} is based on the normal approximation of the corresponding segregation coefficient.
#' That is, each performs the segregation coefficient tests which are appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' The segregation coefficients in the multi-class case are the extension of Pielou's segregation coefficient
#' for the two-class case.
#' (See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).
#' 
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, lower and upper confidence levels, sample estimates (i.e. observed values) and null value
#' (i.e. expected value, which is 0) for the segregation coefficients
#' and also names of the test statistics, estimates, null value and the method and the data set used.
#'
#' The null hypothesis for each cell \eqn{i,j} is that the corresponding segregation coefficient equal to the expected value
#' (which is 0) under RL or CSR.
#'
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Zseg.coeff.ct} only 
#' @param VarSC The variance matrix for the segregation coefficients in the NNCT, \code{ct} ; used in \code{Zseg.coeff.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zseg.coeff} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zseg.coeff} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, for the segregation 
#' coefficients
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Zseg.coeff} only
#'
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of test statistics for the segregation coefficients}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the segregation coefficients at the given
#' confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{conf.int}{Confidence interval for segregation coefficients, it is \code{NULL} here since we provide the upper 
#' and lower confidence limits as \eqn{k \times k} matrices.} 
#' \item{cnf.lvl}{Level of the upper and lower confidence limits of the segregation coefficients,
#' provided in \code{conf.level}.}
#' \item{estimate}{Estimate of the parameter, i.e., matrix of the observed segregation coefficients}
#' \item{est.name,est.name2}{Names of the estimates, both are same in this function}
#' \item{null.value}{Hypothesized null values for the parameters, i.e. expected values of the segregation 
#' coefficients, which are all 0 under RL or CSR.}
#' \item{null.name}{Name of the null value}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zseg.coeff.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zseg.coeff} only}
#' 
#' @seealso \code{\link{seg.coeff}} and \code{\link{Zseg.ind}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZseg.coeff
NULL
#'
#' @rdname funsZseg.coeff
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' varT<-var.seg.coeff(ct,covN)
#'
#' Zseg.coeff(Y,cls)
#' Zseg.coeff.ct(ct,varT)
#'
#' Zseg.coeff(Y,cls,method="max")
#'
#' Zseg.coeff(Y,cls,alt="g")
#' Zseg.coeff.ct(ct,varT,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Zseg.coeff.ct(ct,varT)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' varT<-var.seg.coeff(ct,covN)
#'
#' Zseg.coeff(Y,cls)
#' Zseg.coeff.ct(ct,varT)
#'
#' Zseg.coeff(Y,cls,alt="g")
#' Zseg.coeff.ct(ct,varT,alt="g")
#'
#' @export
Zseg.coeff.ct <- function(ct,VarSC,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  k<-nrow(ct)
  
  estimate<-scf<-seg.coeff(ct)
  estimate.name <-c("segregation coefficients for pairs of classes")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  nullij<-0
  names.null <-"segregation coefficients"
  
  ind.mat1<-ind.seg.coeff(k)
  
  kl<-k*(k+1)/2
  Tvec<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]; c1<-rc1[2]
    Tvec<-c(Tvec,scf[r1,c1]) #entry in ct corresponding to index i of T
  }
  ZT<- Tvec/sqrt(VarSC)
  
  if (all(is.na(ZT)))
  {stop('The test statistics are all NaN, so the segregation coefficient tests are not defined')}
  
  TS<-stderr0<- matrix(0,k,k);
  stderr0[lower.tri(stderr0, diag=TRUE)] <- sqrt(VarSC)
  stderr0 <- t(stderr0)
  TS[lower.tri(TS, diag=TRUE)] <- ZT
  tstat <- t(TS)
  
  ts<-tstat+t(tstat)-diag(diag(tstat),k,k)
  stderr<-stderr0+t(stderr0)-diag(diag(stderr0),k,k)
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderr
           ucl <-estimate+crit.val*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("Z-Tests for Segregation Coefficients between Pairs of Classes")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  
  if (!is.null(lcl)) {rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) {rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"z-tests for segregation coefficients"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value = nullij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funsZseg.coeff
#'
#' @export
Zseg.coeff <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  varT<-var.seg.coeff(ct,covN) 
  
  rval<- Zseg.coeff.ct(ct,varT,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsXsq.seg.coeff
#'
#' @title Chi-square Test for Segregation Coefficients
#' 
#' @description
#' Two functions: \code{Xsq.seg.coeff.ct} and \code{Xsq.seg.coeff}.
#' 
#' Each one performs hypothesis tests of (simultaneous) equality of the segregation coefficients in an NNCT
#' to the ones under RL or CSR.
#' That is, each performs the combined Chi-square test for segregation coefficients which is appropriate 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).
#' 
#' Each test is based on the Chi-square approximation of the corresponding quadratic form for the segregation
#' coefficients in an NNCT.
#' The segregation coefficients in the multi-class case are the extension of Pielou's segregation coefficient
#' for the two-class case.
#' (See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail).
#' 
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{k(k+1)/2-1}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of the segregation coefficients in the NNCT
#' (which are 0 for this function) and also the sample estimates (i.e. observed values) of the segregation
#' coefficients. The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis for all cells \eqn{(i,j)} is that the corresponding segregation coefficients are all 
#' equal to the expected value (which is 0) under RL or CSR.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Xsq.seg.coeff.ct} only 
#' @param covSC The covariance matrix for the segregation coefficients in the NNCT, used in \code{Xsq.seg.coeff.ct} only.
#' Usually output of the function \code{\link{cov.seg.coeff}}
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Xsq.seg.coeff} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Xsq.seg.coeff} only
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Xsq.seg.coeff} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for the combined segregation coefficients}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{k(k+1)/2-1} for this function.}
#' \item{estimate}{The \code{vector} of estimates of the parameters, i.e., observed values of segregation coefficients 
#' in the NNCT.}
#' \item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
#' \item{null.value}{The null value of the parameters, i.e., expected values of segregation coefficients
#' in the NNCT under RL or CSR (which is 0).}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Xsq.seg.coeff.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Xsq.seg.coeff} only}
#'  
#' @seealso \code{\link{seg.coeff}}, \code{\link{Zseg.coeff.ct}} and \code{\link{Zseg.coeff}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsXsq.seg.coeff
NULL
#'
#' @rdname funsXsq.seg.coeff
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' covSC<-cov.seg.coeff(ct,covN)
#'
#' Xsq.seg.coeff(Y,cls)
#' Xsq.seg.coeff.ct(ct,covSC)
#'
#' Xsq.seg.coeff(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Xsq.seg.coeff.ct(ct,covSC)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ipd<-ipd.mat(Y)
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' covSC<-cov.seg.coeff(ct,covN)
#'
#' Xsq.seg.coeff(Y,cls)
#' Xsq.seg.coeff.ct(ct,covSC)
#'
#' @export
Xsq.seg.coeff.ct <- function(ct,covSC)
{
  k<-nrow(ct)
  ind.mat1<-ind.seg.coeff(k)
  
  estimate<-scf<-seg.coeff(ct)
  
  kl<-k*(k+1)/2
  Tvec<-vector()
  for (i in 1:kl)
  {
    rc1<-ind.mat1[i,] #row and column indices for the \code{T} vector
    r1<-rc1[1]; c1<-rc1[2]
    Tvec<-c(Tvec,scf[r1,c1]) #entry in ct corresponding to index i of T
  }
  
  ts<- t(Tvec) %*% ginv(covSC,tol=1.490116e-20) %*% (Tvec)
  nu<-kl-1
  pval<- 1-pchisq(ts,df=nu)
  
  method <-"Combined Test for Segregation Coefficients"
  
  dname <-deparse(substitute(ct))
  
  estimate.name <-c("segregation coefficients for pairs of classes")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value = 0,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsXsq.seg.coeff
#'
#' @export
Xsq.seg.coeff <- function(dat,lab,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv) #default is byrow
  
  covSC<-cov.seg.coeff(ct,covN) 
  
  rval<-Xsq.seg.coeff.ct(ct,covSC)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funs.cell.spec.ss
#'
#' @title Pielou's Cell-specific Segregation Test with Normal Approximation (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{cell.spec.ss.ct} and \code{cell.spec.ss}.
#' 
#' Both functions are objects of class \code{"cellhtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected values of the 
#' cell counts (i.e., entries) in the NNCT for \eqn{k \ge 2} classes.
#' Each test is appropriate (i.e. have the appropriate asymptotic sampling distribution)
#' when that data is obtained by sparse sampling.
#' 
#' Each cell-specific segregation test is based on the normal approximation of the entries
#' in the NNCT and are due to \insertCite{pielou:1961;textual}{nnspat}.
#'
#' Each function yields a contingency table of the test statistics, \eqn{p}-values for the corresponding 
#' alternative, expected values, lower and upper confidence levels, sample estimates (i.e. observed values)
#' and null value(s) (i.e. expected values) for the \eqn{N_{ij}} values for \eqn{i,j=1,2,\ldots,k} and also names of the test
#' statistics, estimates, null values and the method and the data set used.
#' 
#' The null hypothesis is that all \eqn{E(N_{ij})=n_i c_j /n} where \eqn{n_i} is the sum of row \eqn{i} (i.e. size of class \eqn{i})
#' \eqn{c_j} is the sum of column \eqn{j} in the \eqn{k \times k} NNCT for \eqn{k \ge 2}.
#' In the output, the test statistic, \eqn{p}-value and the lower and upper confidence limits are valid only 
#' for (properly) sparsely sampled data.
#' 
#' See also
#' (\insertCite{pielou:1961,ceyhan:eest-2010;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{cell.spec.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{cell.spec.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{cell.spec.ss} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the entries, \eqn{N_{ij}} in the NNCT
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{cell.spec.ss} only
#'    
#' @return A \code{list} with the elements
#' \item{statistic}{The \code{matrix} of \eqn{Z} test statistics for cell-specific tests}
#' \item{stat.names}{Name of the test statistics}
#' \item{p.value}{The \code{matrix} of \eqn{p}-values for the hypothesis test for the corresponding alternative}
#' \item{LCL,UCL}{Matrix of Lower and Upper Confidence Levels for the entries \eqn{N_{ij}} in the NNCT 
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{conf.int}{The confidence interval for the estimates, it is \code{NULL} here, since we provide the \code{UCL} and \code{LCL}
#' in \code{matrix} form.}
#' \item{cnf.lvl}{Level of the upper and lower confidence limits (i.e., conf.level) of the NNCT entries.}
#' \item{estimate}{Estimates of the parameters, i.e., matrix of the NNCT entries of the \eqn{k \times k} NNCT, Nij
#' for i,j=1,2,\ldots,k.}
#' \item{est.name,est.name2}{Names of the estimates, former is a shorter description of the estimates
#' than the latter.}
#' \item{null.value}{Hypothesized null value for the expected values of the NNCT entries, 
#' E(Nij) for i,j=1,2,\ldots,k.}
#' \item{null.name}{Name of the null values}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{cell.spec.ss.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{cell.spec.ss} only}
#'  
#' @seealso \code{\link{cell.spec.ct}} and \code{\link{cell.spec}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funs.cell.spec.ss
NULL
#'
#' @rdname funs.cell.spec.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' cell.spec.ss(Y,cls)
#' cell.spec.ss.ct(ct)
#' cell.spec.ss.ct(ct,alt="g")
#'
#' cell.spec.ss(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' cell.spec.ss(Y,fcls)
#' cell.spec.ss.ct(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' cell.spec.ss(Y,cls,alt="l")
#' cell.spec.ss.ct(ct)
#' cell.spec.ss.ct(ct,alt="l")
#'
#' @export
cell.spec.ss.ct <- function(ct,alternative=c("two.sided", "less", "greater"),conf.level = 0.95) 
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  estimate<-ct
  estimate.name <-c("NNCT entries")
  
  k<-nrow(ct)
  rs <- row.sum(ct); cs <- col.sum(ct)
  n<-sum(rs)
  
  ts<- stderrN<- Exp<- matrix(0,k,k);
  for (i in 1:k)
    for (j in 1:k)
    { stderrN[i,j]<-sqrt(ct[i,j]*(n-ct[i,j])/n)
    Exp[i,j]<-rs[i]*cs[j]/n
    ts[i,j] <- (ct[i,j]-Exp[i,j])/stderrN[i,j]
    }
  
  if (all(is.na(ts)))
  {stop('All of the test stat statistics are NaN, these cell-specific tests are not well defined')}
  
  nullNij<-Exp
  names.null <-"NNCT entries"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           lcl <-NULL
           ucl <-estimate+qnorm(conf.level)*stderrN
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ucl <-NULL
           lcl <-estimate-qnorm(conf.level)*stderrN
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           crit.val <-qnorm(1-alpha/2)
           lcl <-estimate-crit.val*stderrN
           ucl <-estimate+crit.val*stderrN
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  cnf.lvl<-conf.level
  
  method <-c("NNCT Cell-Specific Tests (for Sparse Sampling)")
  
  clnames<-rownames(ct) #row and column names for the NNCT, \code{ct} 
  rownames(ts)<-colnames(ts)<-clnames #row and column names for the segregation indices test stat matrix
  rownames(pval)<-colnames(pval)<-clnames
  rownames(nullNij)<-colnames(nullNij)<-clnames
  if (!is.null(lcl)) {rownames(lcl)<-colnames(lcl)<-clnames}
  if (!is.null(ucl)) {rownames(ucl)<-colnames(ucl)<-clnames}
  ts.names <-"(NNCT) cell-specific segregation test statistics (for sparse sampling)"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    stat.names=ts.names,
    p.value=pval,
    LCL = lcl,UCL = ucl,
    conf.int = NULL,
    cnf.lvl=conf.level,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name, #this is for other functions to have a different description for the sample estimates
    null.value = nullNij,
    null.name=names.null,
    alternative = alternative,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"cellhtest"
  return(rval)
} #end for the function
#'
#' @rdname funs.cell.spec.ss
#'
#' @export
cell.spec.ss <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...) 
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-cell.spec.ss.ct(ct,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsPseg.ss
#'
#' @title Pielou's Overall Test of Segregation for NNCT (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{Pseg.ss.ct} and \code{Pseg.ss}.
#'
#' Both functions are objects of class \code{"Chisqtest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of deviations of 
#' cell counts from the expected values under independence for all cells (i.e., entries) combined in the NNCT.
#' That is, each test is Pielou's overall test of segregation based on NNCTs for \eqn{k \ge 2} classes.
#' This overall test is based on the chi-squared approximation,
#' is equivalent to Pearson's chi-squared test on NNCT and
#' is due to \insertCite{pielou:1961;textual}{nnspat}.
#' Each test is appropriate (i.e. have the appropriate asymptotic sampling distribution)
#' when that data is obtained by sparse sampling.
#' 
#' Each function yields the test statistic, \eqn{p}-value and \code{df} which is \eqn{(k-1)^2}, description of the 
#' alternative with the corresponding null values (i.e. expected values) of NNCT entries,
#' sample estimates (i.e. observed values) of the entries in NNCT.
#' The functions also provide names of the test statistics, the method and the data set used.
#' 
#' The null hypothesis is that \eqn{E(N_{ij})=n_i c_j /n} for all entries in the NNCT
#' where \eqn{n_i} is the sum of row \eqn{i} (i.e. size of class \eqn{i}), \eqn{c_j} is the sum of column \eqn{j} in the \eqn{k \times k} NNCT for \eqn{k \ge 2}.
#' In the output, the test statistic and the \eqn{p}-value are valid only 
#' for (properly) sparsely sampled data.
#'
#' See also (\insertCite{pielou:1961,ceyhan:eest-2010;textual}{nnspat})
#' and the references therein.
#' 
#' @param ct A nearest neighbor contingency table, used in \code{Pseg.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Pseg.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Pseg.ss} only
#' @param yates A logical parameter (default=\code{TRUE}). If \code{TRUE}, Yates continuity correction is applied,
#' and if \code{FALSE} the continuity correction is not applied. 
#' Equivalent to the \code{correct} argument in the base function \code{\link[stats]{chisq.test}}
#' @param sim A logical parameter (default=\code{FALSE}). If \code{TRUE}, \eqn{p}-values are computed by Monte Carlo simulation
#' and if \code{FALSE} the \eqn{p}-value is based on the chi-squared approximation.
#' Equivalent to the \code{simulate.p.value} argument in the base function \code{\link[stats]{chisq.test}}
#' @param Nsim A positive integer specifying the number of replicates used in the Monte Carlo test.
#' Equivalent to the \code{B} argument in the base function \code{\link[stats]{chisq.test}}
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. 
#' used in \code{Pseg.ss} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The overall chi-squared statistic}
#' \item{stat.names}{Name of the test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is (k-1)^2 for this function. 
#' Yields \code{NA} if \code{sim=TRUE} and \code{NSim} is provided.}
#' \item{estimate}{Estimates of the parameters, NNCT, i.e., matrix of the observed \eqn{N_{ij}} values
#' which is the NNCT.}
#' \item{est.name,est.name2}{Names of the estimates, they are identical for this function.}
#' \item{null.value}{Matrix of hypothesized null values for the parameters which are expected values of the
#' the \eqn{N_{ij}} values in the NNCT.}
#' \item{null.name}{Name of the null values}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Pseg.ss.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Pseg.ss} only}
#' 
#' @seealso \code{\link{overall.nnct.ct}}, \code{\link{overall.nnct}}, \code{\link{overall.seg.ct}},
#' \code{\link{overall.seg}} and \code{\link[stats]{chisq.test}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsPseg.ss
NULL
#'
#' @rdname funsPseg.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' Pseg.ss(Y,cls)
#'
#' Pseg.ss.ct(ct)
#' Pseg.ss.ct(ct,yates=FALSE)
#'
#' Pseg.ss.ct(ct,yates=FALSE,sim=TRUE)
#' Pseg.ss.ct(ct,yates=FALSE,sim=TRUE,Nsim=10000)
#'
#' Pseg.ss(Y,cls,method="max")
#' Pseg.ss(Y,cls,yates=FALSE,sim=TRUE,Nsim=10000,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Pseg.ss(Y,fcls)
#' Pseg.ss.ct(ct)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:4,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' Pseg.ss(Y,cls)
#' Pseg.ss.ct(ct,yates=FALSE)
#'
#' Pseg.ss(Y,cls, sim = TRUE, Nsim = 2000)
#' Pseg.ss.ct(ct,yates=FALSE) 
#'
#' @export
Pseg.ss.ct <- function(ct,yates=TRUE,sim = FALSE, Nsim = 2000)
{
  TS<-chisq.test(ct,correct=yates,simulate.p.value=sim,B=Nsim)
  ts<-TS$stat
  nu<-TS$par
  pval<-TS$p.val
  
  ifelse(yates==TRUE,method <-"Pielou's Test of Segregation with Yates correction (for Sparse Sampling)",
         method <-"Pielou's Test of Segregation without Yates correction (for Sparse Sampling)")
  
  dname <-TS$data.name
  
  estimate<-ct
  estimate.name <-c("NNCT entries")
  
  EN<-TS$exp
  
  clnames<-rownames(ct) #row and column names for the nnct, ct
  rownames(estimate)<-colnames(estimate)<-clnames #row and column names for the difference matrix
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name,
    null.value = EN,
    method = method,
    ct.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'
#' @rdname funsPseg.ss
#'
#' @export
Pseg.ss <- function(dat,lab,yates=TRUE,sim = FALSE, Nsim = 2000,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-Pseg.ss.ct(ct,yates,sim,Nsim)
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZdir.nnct.ss
#'
#' @title Directional Segregation Test for Two Classes with Normal Approximation (for Sparse Sampling)
#'
#' @description
#' Two functions: \code{Zdir.nnct.ss.ct} and \code{Zdir.nnct.ss}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of independence in the \eqn{2 \times 2} NNCT which implies \eqn{Z_P=0}
#' or equivalently \eqn{N_{11}/n_1=N_{21}/n_2}.
#' \eqn{Z_P=(N_{11}/n_1-N_{21}/n_2)\sqrt{n_1 n_2 n/(C_1 C_2)}} 
#' where \eqn{N_{ij}} is the cell count in entry \eqn{i,j}, \eqn{n_i} is the sum of row \eqn{i} (i.e. size of class \eqn{i}),
#' \eqn{c_j} is the sum of column \eqn{j} in the \eqn{2 \times 2} NNCT;
#' \eqn{N_{11}/n_1} and \eqn{N_{21}/n_2} are also referred to as the phat estimates in row-wise binomial framework
#' for \eqn{2 \times 2} NNCT (see \insertCite{ceyhan:jnps-NNCT-2010;textual}{nnspat}).
#' 
#' That is, each performs directional (i.e. one-sided) tests based on the \eqn{2 \times 2} NNCT and is appropriate
#' (i.e. have the appropriate asymptotic sampling distribution)
#' when that data is obtained by sparse sampling.
#' (See \insertCite{ceyhan:jnps-NNCT-2010;textual}{nnspat} for more detail).
#' 
#' Each test is based on the normal approximation of \eqn{Z_P} which is the directional \eqn{Z}-tests for the chi-squared
#' tests of independence for the contingency tables \insertCite{bickel:1977}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, sample estimate (i.e. observed value) and
#' null (i.e., expected) value for the difference in the phat values (which is 0 for this test) in an NNCT, 
#' and method and name of the data set used.
#' 
#' The null hypothesis is that \eqn{E[Z_P] = 0} or equivalently \eqn{N_{11}/n_1 = N_{21}/n_2}.
#' 
#' @param ct The NNCT, used in \code{Zdir.nnct.ss.ct} only 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zdir.nnct.ss} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zdir.nnct.ss} only
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the confidence limits, default is \code{0.95}, 
#' for the difference in phat values in the NNCT
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Zdir.nnct.ss} only
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the directional (i.e. one-sided) test of segregation based on
#' the NNCT}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the difference in phat values in the NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the observed difference in phat values in the NNCT.}
#' \item{null.value}{Hypothesized null value for the difference in phat values in the NNCT
#' which is 0 for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zdir.nnct.ss.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zdir.nnct.ss} only}
#'  
#' @seealso \code{\link{Zdir.nnct.ct}}, \code{\link{Zdir.nnct}}, \code{\link{Pseg.ss.ct}} and \code{\link{Pseg.ss}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZdir.nnct.ss
NULL
#'
#' @rdname funsZdir.nnct.ss
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' Zdir.nnct.ss(Y,cls)
#' Zdir.nnct.ss.ct(ct)
#' Zdir.nnct.ss(Y,cls,alt="g")
#'
#' Zdir.nnct.ss(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Zdir.nnct.ss(Y,fcls)
#' Zdir.nnct.ss.ct(ct)
#'
#' #############
#' ct<-matrix(1:4,ncol=2)
#' Zdir.nnct.ss.ct(ct) #gives an error message if ct<-matrix(1:9,ncol=3)
#'
#' @export
Zdir.nnct.ss.ct <- function(ct,alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  rs <- row.sum(ct); cs <- col.sum(ct)
  n<-sum(rs)
  stderr <-sqrt((cs[1]*cs[2])/(rs[1]*rs[2]*n))
  estimate<-(ct[1,1]/rs[1]-ct[2,1]/rs[2])
  ts<-estimate/stderr
  
  names(ts) <-"directional Z-test statistic for 2 x 2 NNCT"
  
  names(estimate) <-c("difference between phat estimates")
  method <-c("Directional Z-Test for 2 x 2 NNCT (for Sparse Sampling)")
  
  null.val<-0
  names(null.val) <-"(expected) difference between phat estimates in row-wise binomial framework for 2 x 2 NNCT"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )

  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-c("contingency table = ",deparse(substitute(ct)))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZdir.nnct.ss
#'
#' @export
Zdir.nnct.ss <- function(dat,lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  rval<-  Zdir.nnct.ss.ct(ct,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

# funsZdir.nnct
#'
#' @title Directional Segregation Test for Two Classes with Normal Approximation
#'
#' @description
#' Two functions: \code{Zdir.nnct.ct} and \code{Zdir.nnct}.
#' 
#' Both functions are objects of class \code{"htest"} but with different arguments (see the parameter list below).
#' Each one performs hypothesis tests of equality of the expected value of the the difference between the
#' phat estimates in a \eqn{2 \times 2} NNCT to the one under RL or CSR (which is \eqn{-1/(n-1)}) where
#' phat estimates are \eqn{N_{11}/n_1} and \eqn{N_{21}/n_2}.
#' That is, each performs directional (i.e. one-sided) tests based on the \eqn{2 \times 2} NNCT 
#' (i.e. have the appropriate asymptotic sampling distribution)
#' for completely mapped data.
#' (See \insertCite{ceyhan:jnps-NNCT-2010;textual}{nnspat} for more detail).
#' 
#' The one-sided (or directional) test has two types, specified with the type argument, with default
#' \code{type="II"}. The second type is 
#' \eqn{Z_{II}=(T_n-E T_n)/\sqrt{Var(T_n)}} where \eqn{T_n=N_{11}/n_1 - N_{21}/n_2}
#' (which is the difference between
#' phat values) and the first type is \eqn{Z_I=U_n T_n} where \eqn{U_n=\sqrt{n_1 n_2/(C_1 C_2)}}.
#' Each test is based on the normal approximation of the \eqn{Z_I} and \eqn{Z_{II}} based on the \eqn{2 \times 2} NNCT and
#' are due to \insertCite{ceyhan:jnps-NNCT-2010}{nnspat}.
#' 
#' Each function yields the test statistic, \eqn{p}-value for the
#' corresponding alternative, the confidence interval, sample estimate (i.e. observed value) and null
#' (i.e., expected) value for the difference in phat values which is \eqn{-1/(n-1)} for this function 
#' and method and name of the data set used.
#' 
#' The null hypothesis is that all \eqn{E[Z_{II}] = 0} and \eqn{E[Z_I]} converges to 0 as class sizes go to infinity (or
#' \eqn{T_n} has mean equal to \eqn{-1/(n-1)} where \eqn{n} is the data size. 
#' 
#' @param ct The NNCT, used in \code{Zdir.nnct.ct} only 
#' @param covN The \eqn{k^2 \times k^2} covariance matrix of row-wise vectorized entries of NNCT
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point,
#' used in \code{Zdir.nnct} only
#' @param lab The \code{vector} of class labels (numerical or categorical), used in \code{Zdir.nnct} only
#' @param type The type of the directional (i.e. one-sided) test with default=\code{"II"}.
#' Takes on values \code{"I"} and \code{"II"} for types I and II directional tests (see the description above).
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the difference in phat estimates in the NNCT
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' used in \code{Zdir.nnct} only
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the directional (i.e. one-sided) test of segregation based on
#' the NNCT}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the difference in phat values in an NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the observed difference in phat values in an NNCT.}
#' \item{null.value}{Hypothesized null value for the difference in phat values in an NNCT
#' which is \eqn{-1/(n-1)} for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{ct.name}{Name of the contingency table, \code{ct}, returned by \code{Zdir.nnct.ct} only}
#' \item{data.name}{Name of the data set, \code{dat}, returned by \code{Zdir.nnct} only}
#'  
#' @seealso \code{\link{Zdir.nnct.ss.ct}}, \code{\link{Zdir.nnct.ss}}, \code{\link{overall.nnct.ct}}
#' and \code{\link{overall.nnct}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZdir.nnct
NULL
#'
#' @rdname funsZdir.nnct
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#'
#' W<-Wmat(ipd)
#' Qv<-Qvec(W)$q
#' Rv<-Rval(W)
#' varN<-var.nnct(ct,Qv,Rv)
#' covN<-cov.nnct(ct,varN,Qv,Rv)
#'
#' Zdir.nnct(Y,cls)
#' Zdir.nnct.ct(ct,covN)
#'
#' Zdir.nnct(Y,cls,alt="g")
#' Zdir.nnct.ct(ct,covN,type="I",alt="l")
#'
#' Zdir.nnct(Y,cls,method="max")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ct<-nnct(ipd,fcls)
#'
#' Zdir.nnct(Y,fcls)
#' Zdir.nnct.ct(ct,covN)
#'
#' #############
#' ct<-matrix(1:4,ncol=2) 
#' Zdir.nnct.ct(ct,covN) #gives an error message if ct is defined as ct<-matrix(1:9,ncol=3)
#'
#' @export
Zdir.nnct.ct <- function(ct,covN,type="II",alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  rs <- row.sum(ct); cs <- col.sum(ct); n<-sum(rs)
  varN11<-covN[1,1]
  varN21<-covN[3,3]
  covN11N21<-covN[1,3]
  VarTn<-varN11/rs[1]^2+varN21/rs[2]^2-2*covN11N21/(rs[1]*rs[2])
  stderr <-sqrt(VarTn)
  estimate<-ct[1,1]/rs[1]-ct[2,1]/rs[2] #Tn value
  Un<-sqrt(rs[1]*rs[2]/(cs[1]*cs[2]))
  ifelse(type=="I",ts<-Un*estimate/stderr,ts<-estimate/stderr)
  
  if (is.na(ts))
  {stop('The test statistic is NaN, so the directional Z-test is not defined')}
  
  names(ts) <-paste("Type ",type," directional Z-test statistic for 2 x 2 NNCT",sep="")
  
  names(estimate) <-c("difference between phat estimates")
  method <-paste("Type ",type," Directional Z-Test for 2 x 2 NNCT",sep="")
  
  null.val<- -1/(n-1) #ETn value
  names(null.val) <-"(expected) difference between phat estimates under RL for 2 x 2 NNCT"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           ifelse(type=="I",cint <-estimate+c(-Inf, qnorm(conf.level))*stderr/Un,
                  cint <-estimate+c(-Inf, qnorm(conf.level))*stderr)
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           ifelse(type=="I",cint <-estimate+c(-qnorm(conf.level),Inf)*stderr/Un,
                  cint <-estimate+c(-qnorm(conf.level),Inf)*stderr)
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           ifelse(type=="I",cint <-estimate+c(-cint, cint)*stderr/Un,cint <-estimate+c(-cint, cint)*stderr)
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'
#' @rdname funsZdir.nnct
#'
#' @export
Zdir.nnct <- function(dat,lab,type="II",alternative=c("two.sided", "less", "greater"),conf.level = 0.95,...)
{
  ipd<-ipd.mat(dat,...)
  ct<-nnct(ipd,lab)
  
  W<-Wmat(ipd)
  Qv<-Qvec(W)$q
  Rv<-Rval(W)
  varN<-var.nnct(ct,Qv,Rv) 
  covN<-cov.nnct(ct,varN,Qv,Rv)
  
  rval<-  Zdir.nnct.ct(ct,covN,type=type,alternative=alternative,conf.level=conf.level) 
  
  dname <-deparse(substitute(dat))
  
  rval$data.name<-dname
  return(rval)
} #end for the function
#'

#################################################################

#' @title Generation of Points Associated in the Type I Sense with a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_2} 2D points associated with the given set of points (i.e. reference points) \eqn{X_1} in the
#' type I fashion with circular (or radial) between class attraction parameter \code{p}, which is a probability value between 0 and 1.
#' The generated points are intended to be from a different class, say class 2 (or \eqn{X_2} points) than the reference
#' (i.e. \eqn{X_1} points, say class 1 points, denoted as \code{X1} as an argument of the function). 
#' To generate \eqn{n_2} (denoted as \code{n2} as an argument of the function) \eqn{X_2} points, \eqn{n_2} of \eqn{X_1} points are randomly selected (possibly with replacement) and
#' for a selected \code{X1} point, say \eqn{x_{1ref}}, a \eqn{Uniform(0,1)} number, \eqn{U}, is generated.
#' If \eqn{U \le p}, a new point from the class 2, say \eqn{x_{2new}}, is generated within a
#' circle with radius equal to the distance to the closest \eqn{X_1} point (uniform in the polar coordinates),
#' else the new point is generated uniformly
#' within the smallest bounding box containing \eqn{X_1} points.
#' That is, if \eqn{U \le p}, \eqn{x_{2new} = x_{1ref}+r_u c(\cos(t_u),\sin(t_u))}
#' where \eqn{r_u \sim U(0,rad)} and \eqn{t_u \sim U(0, 2\pi)} with \eqn{rad=\min(d(x_{1ref},X_1\setminus \{x_{1ref}\}))},
#' else \eqn{x_{2new} \sim rect(X_1)} where \eqn{rect(X_1)} is the smallest bounding box containing \eqn{X_1} points.
#' Note that, the level of association increases as \code{p} increases, and the association vanishes
#' when \code{p} approaches to 0.
#'
#' Type I association is closely related to Type C association in
#' \insertCite{ceyhan:serra-2014;textual}{nnspat}, see the function \code{\link{rassocC}}
#' and also other association types.
#' In the type C association pattern
#' the new point from the class 2, \eqn{x_{2new}}, is generated (uniform in the polar coordinates) within a circle
#' centered at \eqn{x_{1ref}} with radius equal to \eqn{r_0},
#' in type U association pattern \eqn{x_{2new}} is generated similarly except it is uniform in the circle.
#' In type G association, \eqn{x_{2new}} is generated from the bivariate normal distribution centered at \eqn{x_{1ref}} with covariance
#' \eqn{\sigma I_2} where \eqn{I_2} is \eqn{2 \times 2} identity matrix. 
#' 
#' @param X1 A set of 2D points representing the reference points, also referred as class 1 points.
#' The generated points are associated in a type I sense (in a circular/radial fashion) with these points.
#' @param n2 A positive integer representing the number of class 2 (i.e. \eqn{X_2}) points to be generated.
#' @param p A real number between 0 and 1 representing the attraction probability of class 2 points associated
#' with a randomly selected class 1 point (see the description below).  
#'
#' @return A \code{list} with the elements 
#' \item{pat.type}{equals \code{"ref.gen"} for the bivariate pattern of association of class 2 (i.e. \eqn{X_2}) points with the reference
#' points (i.e. \eqn{X_1}), indicates reference points are required to be entered as an argument in the function}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{Radial (i.e. circular) between class attraction parameter controlling the level of association}
#' \item{gen.points}{The output set of generated points (i.e. class 2 points) associated with reference (i.e.
#' \eqn{X_1} points)}
#' \item{ref.points}{The input set of reference points \eqn{X_1}, i.e., points with which generated class 2
#' points are associated.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{lab}{The class labels of the generated points, it is \code{NULL} for this function, since only class 2
#' points are generated in this pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 2 points and
#' the number of reference (i.e. \eqn{X_1}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and
#' the reference points}
#'
#' @seealso \code{\link{rassocC}}, \code{\link{rassocG}}, \code{\link{rassocU}}, and \code{\link{rassoc}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20; n2<-1000;  #try also n1<-10; n2<-1000;
#'
#' p<- .75 #try also .25, .5, .9, runif(1)
#' #with default bounding box (i.e., unit square)
#' X1<-cbind(runif(n1),runif(n1))  #try also X1<-1+cbind(runif(n1),runif(n1))
#'
#' Xdat<-rassocI(X1,n2,p)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' @export 
rassocI <- function(X1,n2,p)
{
  ipdx1<-ipd.mat(X1)
  n1<-nrow(X1)
  
  rx<-range(X1[,1])
  ry<-range(X1[,2])
  
  X2<-matrix(0,n2,2)
  for (i in 1:n2)
  {ind<-sample(1:n1,1)
  rmax<- min(ipdx1[ind,]) #this is equivalent to case with r0=rmax for Type C association in SERRA article
  U<-runif(1)
  if ( U <= p )
  {
    ru<-runif(1,0,rmax) 
    tu<-runif(1,0,2*pi)
    X2[i,]<-X1[ind,]+ru*c(cos(tu),sin(tu))
  }
  else
  { X2[i,]<-c(runif(1,rx[1],rx[2]),runif(1,ry[1],ry[2])) }
  }
  
  Xlim<-range(X1[,1],X2[,1])
  Ylim<-range(X1[,2],X2[,2])
  
  pname <-deparse(substitute(p))
  param<-p
  
  rparam<-round(param,2)
  
  names(param)<-"attraction probability"
  typ<-paste("Type I Association of ",n2, " Class 2 points with ", n1," Class 1 points with circular (or radial) between class attraction parameter = ",p,sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  txt<-"Type I Association with the Reference Points"
  main.txt<-paste("Type I Association of Two Classes \n with Attraction Parameter p = ",rparam,sep="")
  
  res<-list(
    pat.type="ref.gen", #ref.gen for bivariate pattern of association of X2 wrt reference points X1
    type=typ,
    parameters=param,
    gen.points=X2, #generated points associated with Y points
    ref.points=X1, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    lab=NULL,
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Points Associated in the Type C Sense with a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_2} 2D points associated with the given set of points (i.e. reference points) \eqn{X_1} in the
#' type C fashion with a radius of association \eqn{r_0} (denoted as \code{r0} as
#' an argument of the function) which is a positive real number.
#' The generated points are intended to be from a different class, say class 2 (or \eqn{X_2} points) than the reference
#' (i.e. \eqn{X_1} points, say class 1 points, denoted as \code{X1} as an argument of the function), say class 1 points). 
#' To generate \eqn{n_2} \eqn{X_2} points, \eqn{n_2} of \eqn{X_1} points are randomly selected (possibly with replacement) and
#' for a selected \code{X1} point, say \eqn{x_{1ref}},
#' a new point from the class 2, say \eqn{x_{2new}}, is generated within a
#' circle with radius equal to \eqn{r_0} (uniform in the polar coordinates).
#' That is, \eqn{x_{2new} = x_{1ref}+r_u c(\cos(t_u),\sin(t_u))}
#' where \eqn{r_u \sim U(0,r_0)} and \eqn{t_u \sim U(0, 2\pi)}.
#' Note that, the level of association increases as \eqn{r_0} decreases, and the association vanishes when \eqn{r_0} is 
#' sufficiently large.
#' 
#' For type C association, it is recommended to take \eqn{r_0 \le 0.25} times length of the shorter
#' edge of a rectangular study region, or take \eqn{r_0 = 1/(k \sqrt{\hat \rho})} with the appropriate choice of \eqn{k} 
#' to get an association pattern more robust to differences in relative abundances
#' (i.e. the choice of \eqn{k} implies \eqn{r_0 \le 0.25} times length of the shorter edge to have alternative patterns more 
#' robust to differences in sample sizes).
#' Here \eqn{\hat \rho} is the 
#' estimated intensity of points in the study region (i.e., # of points divided by the area of the region). 
#'
#' Type C association is closely related to Type U association, see the function \code{\link{rassocC}}
#' and the other association types.
#' In the type U association pattern
#' the new point from the class 2, \eqn{x_{2new}}, is generated uniformly within a circle
#' centered at \eqn{x_{1ref}} with radius equal to \eqn{r_0}.
#' In type G association, \eqn{x_{2new}} is generated from the bivariate normal distribution centered at \eqn{x_{1ref}} with covariance
#' \eqn{\sigma I_2} where \eqn{I_2} is \eqn{2 \times 2} identity matrix. 
#' In type I association, first a \eqn{Uniform(0,1)} number, \eqn{U}, is generated.
#' If \eqn{U \le p}, \eqn{x_{2new}} is generated (uniform in the polar coordinates) within a
#' circle with radius equal to the distance to the closest \eqn{X_1} point,
#' else it is generated uniformly within the smallest bounding box containing \eqn{X_1} points.
#' 
#' See \insertCite{ceyhan:serra-2014;textual}{nnspat} for more detail.
#'
#' @param X1 A set of 2D points representing the reference points, also referred as class 1 points.
#' The generated points are associated in a type C sense (in a circular/radial fashion) with these points.
#' @param n2 A positive integer representing the number of class 2 points to be generated.
#' @param r0 A positive real number representing the radius of association of class 2 points associated with a
#' randomly selected class 1 point (see the description below). 
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{=\code{"ref.gen"} for the bivariate pattern of association of class 2 points with the reference points
#' (i.e. \eqn{X_1}), indicates reference points are required to be entered as an argument in the function}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{Radius of association controlling the level of association}
#' \item{gen.points}{The output set of generated points (i.e. class 2 points) associated with reference (i.e.
#' \eqn{X_1} points)}
#' \item{ref.points}{The input set of reference points \eqn{X_1}, i.e., points with which generated class 2 points
#' are associated.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 2 points and
#' the number of reference (i.e. \eqn{X_1}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rassocI}}, \code{\link{rassocG}}, \code{\link{rassocU}}, and \code{\link{rassoc}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20; n2<-1000;  #try also n1<-10; n2<-1000;
#'
#' r0<-.15 #try also .10 and .20, runif(1)
#' #with default bounding box (i.e., unit square)
#' X1<-cbind(runif(n1),runif(n1))  #try also X1<-1+cbind(runif(n1),runif(n1))
#'
#' Xdat<-rassocC(X1,n2,r0)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #radius adjusted with the expected NN distance
#' x<-range(X1[,1]); y<-range(X1[,2])
#' ar<-(y[2]-y[1])*(x[2]-x[1]) #area of the smallest rectangular window containing X1 points
#' rho<-n1/ar
#' r0<-1/(2*sqrt(rho)) #r0=1/(2rho) where \code{rho} is the intensity of X1 points
#' Xdat<-rassocC(X1,n2,r0)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' @export 
rassocC <- function(X1,n2,r0)
{
  n1<-nrow(X1)
  
  ind<-sample(1:n1,max(n1,n2),replace=TRUE)
  
  X2<-matrix(0,n2,2)
  for (i in 1:n2)
  {
    ru<-runif(1,0,r0) 
    tu<-runif(1,0,2*pi)
    X2[i,]<-X1[ind[i],]+ru*c(cos(tu),sin(tu))
  }
  
  Xlim<-range(X1[,1],X2[,1])
  Ylim<-range(X1[,2],X2[,2])
  
  pname <-deparse(substitute(r0))
  param<-r0
  names(param)<-"radius of association"
  typ<-paste("Type C Association of ",n2, " Class 2 points with ", n1," Class 1 points with radius = ",r0,sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  txt<-"Type C Association with the Reference Points"
  main.txt<-paste("Type C Association of Two Classes \n with radius = ",round(r0,2),sep="")
  res<-list(
    pat.type="ref.gen", #ref.gen for bivariate pattern of association of X2 wrt reference points X1
    type=typ,
    parameters=param,
    gen.points=X2, #generated points associated with Y points
    ref.points=X1, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Points Associated in the Type U Sense with a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_2} 2D points associated with the given set of points (i.e. reference points) \eqn{X_1} in the
#' type U fashion with a radius of association \eqn{r_0} (denoted as \code{r0} as an argument of the function) which is a positive real number.
#' The generated points are intended to be from a different class, say class 2 (or \eqn{X_2} points) than the reference
#' (i.e. \eqn{X_1} points, say class 1 points, denoted as \code{X1} as an argument of the function), say class 1 points). 
#' To generate \eqn{n_2} (denoted as \code{n2} as an argument of the function)\eqn{X_2} points, \eqn{n_2} of \eqn{X_1} points are randomly selected (possibly with replacement) and
#' for a selected \code{X1} point, say \eqn{x_{1ref}},
#' a new point from the class 2, say \eqn{x_{2new}}, is generated uniformly within a
#' circle with radius equal to \eqn{r_0}.
#' That is, \eqn{x_{2new} = x_{1ref}+r_u c(\cos(t_u),\sin(t_u))}
#' where\eqn{r_u=sqrt(U)*r_0} with \eqn{U \sim U(0,1)} and \eqn{t_u \sim U(0, 2\pi)}.
#' Note that, the level of association increases as \eqn{r_0} decreases, and the association vanishes when \eqn{r_0} is 
#' sufficiently large.
#' 
#' For type U association, it is recommended to take \eqn{r_0 \le 0.10} times length of the shorter
#' edge of a rectangular study region, or take \eqn{r_0 = 1/(k \sqrt{\hat \rho})} with the appropriate choice of \eqn{k} 
#' to get an association pattern more robust to differences in relative abundances
#' (i.e. the choice of \eqn{k} implies \eqn{r_0 \le 0.10} times length of the shorter edge to have alternative patterns more 
#' robust to differences in sample sizes).
#' Here \eqn{\hat \rho} is the 
#' estimated intensity of points in the study region (i.e., # of points divided by the area of the region). 
#'
#' Type U association is closely related to Type C association, see the function \code{\link{rassocC}}
#' and the other association types.
#' In the type C association pattern
#' the new point from the class 2, \eqn{x_{2new}}, is generated (uniform in the polar coordinates) within a circle
#' centered at \eqn{x_{1ref}} with radius equal to \eqn{r_0}.
#' In type G association, \eqn{x_{2new}} is generated from the bivariate normal distribution centered at \eqn{x_{1ref}} with covariance
#' \eqn{\sigma I_2} where \eqn{I_2} is \eqn{2 \times 2} identity matrix.
#' In type I association, first a \eqn{Uniform(0,1)} number, \eqn{U}, is generated.
#' If \eqn{U \le p}, \eqn{x_{2new}} is generated (uniform in the polar coordinates) within a
#' circle with radius equal to the distance to the closest \eqn{X_1} point,
#' else it is generated uniformly within the smallest bounding box containing \eqn{X_1} points.  
#' 
#' See \insertCite{ceyhan:serra-2014;textual}{nnspat} for more detail.
#'
#' @param X1 A set of 2D points representing the reference points, also referred as class 1 points.
#' The generated points are associated in a type U sense (in a circular/radial fashion) with these points.
#' @param n2 A positive integer representing the number of class 2 points to be generated.
#' @param r0 A positive real number representing the radius of association of class 2 points associated with a
#' randomly selected class 1 point (see the description below).
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{=\code{"ref.gen"} for the bivariate pattern of association of class 2 points with the reference points
#' (i.e. \eqn{X_1}), indicates reference points are required to be entered as an argument in the function}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{Radius of association controlling the level of association}
#' \item{gen.points}{The output set of generated points (i.e. class 2 points) associated with reference (i.e.
#' \eqn{X_1} points)}
#' \item{ref.points}{The input set of reference points \eqn{X_1}, i.e., points with which generated class 2 points
#' are associated.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 2 points and
#' the number of reference (i.e. \eqn{X_1}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rassocI}}, \code{\link{rassocG}}, \code{\link{rassocC}}, and \code{\link{rassoc}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20; n2<-1000;  #try also n1<-10; n2<-1000;
#'
#' r0<-.15 #try also .10 and .20
#' #with default bounding box (i.e., unit square)
#' X1<-cbind(runif(n1),runif(n1))  #try also X1<-1+cbind(runif(n1),runif(n1))
#'
#' Xdat<-rassocU(X1,n2,r0)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #radius adjusted with the expected NN distance
#' x<-range(X1[,1]); y<-range(X1[,2])
#' ar<-(y[2]-y[1])*(x[2]-x[1]) #area of the smallest rectangular window containing X1 points
#' rho<-n1/ar
#' r0<-1/(2*sqrt(rho)) #r0=1/(2rho) where \code{rho} is the intensity of X1 points
#' Xdat<-rassocU(X1,n2,r0)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rassocU <- function(X1,n2,r0)
{
  n1<-nrow(X1)
  
  ind<-sample(1:n1,max(n1,n2),replace=TRUE)
  
  X2<-matrix(0,n2,2)
  for (i in 1:n2)
  {
    U<-runif(1) 
    ru<-sqrt(U)*r0 
    tu<-runif(1,0,2*pi)
    X2[i,]<-X1[ind[i],]+ru*c(cos(tu),sin(tu))
  }
  
  Xlim<-range(X1[,1],X2[,1])
  Ylim<-range(X1[,2],X2[,2])
  
  pname <-deparse(substitute(r0))
  param<-r0
  names(param)<-"radius of association"
  typ<-paste("Type U Association of ",n2, " Class 2 points with ", n1," Class 1 points with radius = ",r0,sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  txt<-"Type U Association with the Reference Points"
  main.txt<-paste("Type U Association of Two Classes \n with radius = ",round(r0,2),sep="")
  res<-list(
    pat.type="ref.gen", #ref.gen for bivariate pattern of association of X2 wrt reference points X1
    type=typ,
    parameters=param,
    gen.points=X2, #generated points associated with Y points
    ref.points=X1, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Points Associated in the Type G Sense with a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_2} 2D points associated with the given set of points (i.e. reference points) \eqn{X_1} in the
#' type G fashion with the parameter sigma which is a positive real number representing the variance of the
#' Gaussian marginals.
#' The generated points are intended to be from a different class, say class 2 (or \eqn{X_2} points) than the reference
#' (i.e. \eqn{X_1} points, say class 1 points, denoted as \code{X1} as an argument of the function), say class 1 points). 
#' To generate \eqn{n_2} (denoted as \code{n2} as an argument of the function)\eqn{X_2} points, \eqn{n_2} of \eqn{X_1} points are randomly selected (possibly with replacement) and
#' for a selected \code{X1} point, say \eqn{x_{1ref}},
#' a new point from the class 2, say \eqn{x_{2new}}, is generated from a bivariate normal distribution centered at \eqn{x_{1ref}}
#' where the covariance matrix of the bivariate normal is a diagonal matrix with sigma in the diagonals.
#' That is, \eqn{x_{2new} = x_{1ref}+V} where \eqn{V \sim BVN((0,0),\sigma I_2)} with \eqn{I_2} being the \eqn{2 \times 2} identity matrix.
#' Note that, the level of association increases as \code{sigma} decreases, and the association vanishes when \code{sigma} 
#' goes to infinity.
#' 
#' For type G association, it is recommended to take \eqn{\sigma \le 0.10}  times length of the shorter
#' edge of a rectangular study region, or take \eqn{r_0 = 1/(k \sqrt{\hat \rho})} with the appropriate choice of \eqn{k} 
#' to get an association pattern more robust to differences in relative abundances
#' (i.e. the choice of \eqn{k} implies \eqn{\sigma \le 0.10} times length of the shorter edge to have alternative patterns more 
#' robust to differences in sample sizes).
#' Here \eqn{\hat \rho} is the 
#' estimated intensity of points in the study region (i.e., # of points divided by the area of the region). 
#' 
#' Type G association is closely related to Types C and U association,
#' see the functions \code{\link{rassocC}} and \code{\link{rassocU}} and
#' the other association types.
#' In the type C association pattern
#' the new point from the class 2, \eqn{x_{2new}}, is generated (uniform in the polar coordinates) within a circle
#' centered at \eqn{x_{1ref}} with radius equal to \eqn{r_0},
#' in type U association pattern \eqn{x_{2new}} is generated similarly except it is uniform in the circle.
#' In type I association, first a \eqn{Uniform(0,1)} number, \eqn{U}, is generated.
#' If \eqn{U \le p}, \eqn{x_{2new}} is generated (uniform in the polar coordinates) within a
#' circle with radius equal to the distance to the closest \eqn{X_1} point,
#' else it is generated uniformly within the smallest bounding box containing \eqn{X_1} points.
#' 
#' See \insertCite{ceyhan:serra-2014;textual}{nnspat} for more detail.
#'
#' @param X1 A set of 2D points representing the reference points, also referred as class 1 points.
#' The generated points are associated in a type G sense with these points.
#' @param n2 A positive integer representing the number of class 2 points to be generated.
#' @param sigma A positive real number representing the variance of the Gaussian marginals, where
#' the bivariate normal distribution has covariance \code{BVN((0,0),sigma*I_2)} with \eqn{I_2} being the \eqn{2 \times 2} identity matrix.
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{=\code{"ref.gen"} for the bivariate pattern of association of class 2 points with the reference points
#' (i.e. \eqn{X_1}), indicates reference points are required to be entered as an argument in the function}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{The variance of the Gaussian marginals controlling the level of association, where
#' the bivariate normal distribution has covariance \eqn{\sigma I_2} with \eqn{I_2} being the \eqn{2 \times 2} identity matrix.}
#' \item{gen.points}{The output set of generated points (i.e. class 2 points) associated with reference (i.e.
#' \eqn{X_1} points)}
#' \item{ref.points}{The input set of reference points \eqn{X_1}, i.e., points with which generated class 2 points
#' are associated.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 2 points and
#' the number of reference (i.e. \eqn{X_1}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rassocI}}, \code{\link{rassocG}}, \code{\link{rassocC}}, and \code{\link{rassoc}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20; n2<-1000;  #try also n1<-10; n2<-1000;
#'
#' stdev<-.05  #try also .075 and .15
#' #with default bounding box (i.e., unit square)
#' X1<-cbind(runif(n1),runif(n1))   #try also X1<-1+cbind(runif(n1),runif(n1))
#'
#' Xdat<-rassocG(X1,n2,stdev)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #sigma adjusted with the expected NN distance
#' x<-range(X1[,1]); y<-range(X1[,2])
#' ar<-(y[2]-y[1])*(x[2]-x[1]) #area of the smallest rectangular window containing X1 points
#' rho<-n1/ar
#' stdev<-1/(4*sqrt(rho)) #r0=1/(2rho) where \code{rho} is the intensity of X1 points
#' Xdat<-rassocG(X1,n2,stdev)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rassocG <- function(X1,n2,sigma)
{
  n1<-nrow(X1)
  
  ind<-sample(1:n1,max(n1,n2),replace=TRUE)
  
  X2<-matrix(0,n2,2)
  for (i in 1:n2)
  {
    x<-rnorm(1, mean = 0, sd = sigma)
    y<-rnorm(1, mean = 0, sd = sigma)
    X2[i,]<-X1[ind[i],]+c(x,y)
  }
  
  Xlim<-range(X1[,1],X2[,1])
  Ylim<-range(X1[,2],X2[,2])
  
  pname <-deparse(substitute(sigma))
  param<-sigma
  names(param)<-"standard deviation"
  typ<-paste("Type G Association of ",n2, " Class 2 points with ", n1," Class 1 points with standard deviation = ",sigma,sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  txt<-"Type G Association with the Reference Points"
  main.txt<-paste("Type G Association of Two Classes \n with Standard Deviation = ",round(sigma,2),sep="")
  res<-list(
    pat.type="ref.gen", #ref.gen for bivariate pattern of association of X2 wrt reference points X1
    type=typ,
    parameters=param,
    gen.points=X2, #generated points associated with Y points
    ref.points=X1, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Points Associated with a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_2} 2D points associated with the given set of points (i.e. reference points) \eqn{X_1} in the
#' \code{type=type} fashion with the parameter=asc.par which specifies the level of association.
#' The generated points are intended to be from a different class, say class 2 (or \eqn{X_2} points) 
#' than the reference (i.e. \eqn{X_1} points, say class 1 points, denoted as \code{X1} as an argument 
#' of the function), say class 1 points). 
#' 
#' To generate \eqn{n_2} (denoted as \code{n2} as an argument of the function)\eqn{X_2} points, 
#' \eqn{n_2} of \eqn{X_1} points are randomly selected (possibly with replacement) and
#' for a selected \code{X1} point, say \eqn{x_{1ref}},
#' a new point from the class 2, say \eqn{x_{2new}}, is generated from a distribution specified
#' by the type argument.
#' 
#' In type I association, i.e., if \code{type="I"}, first a \eqn{Uniform(0,1)} number, \eqn{U}, is generated.
#' If \eqn{U \le p}, \eqn{x_{2new}} is generated (uniform in the polar coordinates) within a
#' circle with radius equal to the distance to the closest \eqn{X_1} point,
#' else it is generated uniformly within the smallest bounding box containing \eqn{X_1} points.
#' 
#' In the type C association pattern
#' the new point from the class 2, \eqn{x_{2new}}, is generated (uniform in the polar coordinates) within a circle
#' centered at \eqn{x_{1ref}} with radius equal to \eqn{r_0},
#' in type U association pattern \eqn{x_{2new}} is generated similarly except it is uniform in the circle.
#' 
#' In type G association, \eqn{x_{2new}} is generated from the bivariate normal distribution centered at \eqn{x_{1ref}} with covariance
#' \eqn{\sigma I_2} where \eqn{I_2} is \eqn{2 \times 2} identity matrix.
#' 
#' See \insertCite{ceyhan:serra-2014;textual}{nnspat} for more detail.
#'
#' @param X1 A set of 2D points representing the reference points, also referred as class 1 points.
#' The generated points are associated in a type=type sense with these points.
#' @param n2 A positive integer representing the number of class 2 points to be generated.
#' @param asc.par A positive real number representing the association parameter. For \code{type="I"},
#' it is attraction probability, \code{p}, of class 2 points associated with a randomly selected class 1 point; 
#' for \code{type="C"} or \code{"U"}, it is the radius of association, \code{r0}, of class 2 points associated with a
#' randomly selected class 1 point;
#' for \code{type="G"}, it is the variance of the Gaussian marginals, where
#' the bivariate normal distribution has covariance \eqn{\sigma I_2} with \eqn{I_2} being the \eqn{2 \times 2} identity matrix.
#' @param type The type of the association pattern. Takes on values \code{"I"}, \code{"C"}, \code{"U"} and \code{"G"} 
#' for types I, C, U and G association patterns (see the description above).
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{=\code{"ref.gen"} for the bivariate pattern of association of class 2 points with the reference points
#' (i.e. \eqn{X_1}), indicates reference points are required to be entered as an argument in the function}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{The \code{asc.par} value specifying the level of association}
#' \item{ref.points}{The input set of reference points \eqn{X_1}, i.e., points with which generated class 2 points
#' are associated.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 2 points and
#' the number of reference (i.e. \eqn{X_1}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rassocI}}, \code{\link{rassocC}}, \code{\link{rassocU}}, and \code{\link{rassocG}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20; n2<-1000;  #try also n1<-10; n2<-1000;
#'
#' #with default bounding box (i.e., unit square)
#' X1<-cbind(runif(n1),runif(n1))
#'
#' Xdat<-rassoc(X1,n2,asc.par=.05,type="G") #try other types as well
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #with type U association
#' Xdat<-rassoc(X1,n2,asc.par=.1,type="U")
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #with type C association
#' Xdat<-rassoc(X1,n2,asc.par=.1,type=2) #2 is for "C"
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rassoc <- function(X1,n2,asc.par,type)
{
  
 res<- switch(type,
         I = { res <- rassocI(X1,n2,asc.par) },
         C = { res <- rassocC(X1,n2,asc.par) },
         U = { res <- rassocU(X1,n2,asc.par) },
         G = { res <- rassocG(X1,n2,asc.par)  }
  )
 
 if (is.null(res)) stop("Enter numbers 1-4 or I, C, U, G in quotes for type")
 
  res
} #end for the function
#'

#################################################################

#' @title Type I Non-Random Labeling of a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Given the set of \eqn{n} points, \code{dat}, in a region, this function assigns \eqn{n_1=}\code{round(n*prop,0)} of them as cases,
#' and the rest as controls with first selecting a point, \eqn{Z_i}, as a case and assigning the
#' label case to the remaining points with infection probabilities 
#' \code{prob=c(prop+((1-prop)*rho)/(1:k))} where \code{rho} is a
#' parameter adjusting the NN dependence of infection probabilities.
#' The number of cases will be \eqn{n_1} on the average if the argument \code{poisson=TRUE}
#' (i.e., \eqn{n_1=}\code{rpois(1,round(n*prop,0))}), otherwise \eqn{n_1=}\code{round(n*prop,0)}.
#' We stop when we first exceed \eqn{n_1} cases. \code{rho} must be between \code{-prop/(1-prop)} and 1 for the infection
#' probabilities to be valid.
#' The \code{init.from.cases} is a logical argument (with default=\code{TRUE}) to determine the initial cases are from the
#' cases or controls (the first initial case is always from controls), so if \code{TRUE}, initial cases (other than
#' the first initial case) are selected randomly among the cases (as if they are contagious), otherwise,
#' they are selected from controls as new cases infecting their \code{k}NNs.
#' otherwise first entry is chosen as the case (or case is recorded as the first entry) in the data set, \code{dat}.  
#' 
#' Algorithmically, first all dat points are treated as non-cases (i.e. controls or healthy subjects).
#' Then the function follows the following steps for labeling of the points:  
#' 
#' step 0: \eqn{n_1} is generated randomly from a Poisson distribution with \code{mean = n*prop}, so that the 
#' average number of cases is \code{n*prop}. 
#' 
#' step 0: \eqn{n_1} is generated randomly from a Poisson distribution with \code{mean = round(n*prop,0)}, so that the 
#' average number of cases will be \code{round(n*prop,0)}
#' if the argument \code{poisson=TRUE}, else \eqn{n_1=}\code{round(n*prop,0)}.
#' 
#' step 1: Initially, one point from dat is selected randomly as a case. In the first round this point is selected
#' from the controls, and the subsequent rounds, it is selected from cases if the argument \code{init.from.cases=TRUE},
#' and from controls otherwise. Then it assigns the label case to the \code{k}NNs among controls of the initial case
#' selected in step 1 with infection probabilities \code{prob=c(prop+((1-prop)*rho)/(1:k))}, see the description for the details
#' of the parameters in the \code{prob}.
#'  
#' step 2: Then this initial case and cases among its \code{k}NNs (possibly all \eqn{k+1} points) in step 2 are removed from
#' the data, and for the remaining control points step 1 is applied where initial point is from cases or control
#' based on the argument init.from.cases.
#' 
#' step 3: The procedure ends when number of cases \eqn{n_c} exceeds \eqn{n_1}, and \eqn{n_c-n_1} of the cases (other than the
#' initial cases) are randomly selected and relabeled as controls, i.e. 0s,
#' so that the number of cases is exactly \eqn{n_1}.
#' 
#' In the output cases are labeled as 1 and controls as 0.
#' Note that the infection probabilities of the \code{k}NNs of each initial case increase
#' with increasing rho, and infection probability decreases for increasing k in the \code{k}NNs. 
#'  
#' See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail where type I non-RL pattern is the 
#' case 1 of non-RL pattern considered in Section 6 with \eqn{n_1} is
#' fixed as a parameter rather than being generated from a Poisson distribution and \code{init=FALSEALSE}.
#' 
#' Although the non-RL pattern is described for the case-control setting, it can be adapted for any two-class
#' setting when it is appropriate to treat one of the classes as cases or one of the classes behave like cases
#' and other class as controls.
#' 
#' @param dat A set of points the non-RL procedure is applied to obtain cases and controls randomly in the 
#' type I fashion (see the description).
#' @param prop A real number between 0 and 1 (inclusive) representing the proportion of new cases (on the average)
#' infected by the initial cases, i.e., number of newly infected cases (in addition to the initial cases) is
#' Poisson with \code{mean=round(n*prop)} where \eqn{n} is the number of points in \code{dat}, if the argument \code{poisson=TRUE},
#' else it is \code{round(n*prop)}.
#' @param k An integer representing the number of NNs considered for each initial case, i.e., \code{k}NNs of each
#' initial case are candidates to be infected to become cases.
#' @param rho A parameter for labeling the \code{k}NNs of each initial case as cases such that \code{k}NNs of each initial case
#' is infected with decreasing probabilities \code{prob=c(prop+((1-prop)*rho)/(1:k))} where
#' \code{rho} has to be between \code{-prop/(1-prop)} and 1 for \code{prob} to be a \code{vector} of probabilities. 
#' @param poisson A logical argument (default is \code{FALSE}) to determine whether the number of cases \eqn{n_1},
#' will be random or fixed. If \code{poisson=TRUE} then the \eqn{n_1} is from a Poisson distribution, 
#' \eqn{n_1=}\code{rpois(1,round(n*prop,0))} 
#' otherwise it is fixed, \eqn{n_1=}\code{round(n*prop,0)}.
#' @param init.from.cases A logical argument (default is \code{TRUE}) to determine whether the initial cases at each
#' round will be take from cases or controls. At first round, the initial cases are taken from controls.
#' And in the subsequent rounds, the initial cases are taken from cases if \code{init.from.cases=TRUE},
#' and from controls otherwise.
#' 
#' @return A \code{list} with the elements 
#' \item{pat.type}{\code{="cc"} for the case-control patterns for RL or non-RL of the given data points, \code{dat}}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{\code{prop}, \code{rho}, and \code{k} values for this non-RL pattern, see the description for these
#' parameters.}
#' \item{dat.points}{The set of points non-RL procedure is applied to obtain cases and controls randomly in the 
#' type I fashion}
#' \item{lab}{The labels of the points as 1 for cases and 0 for controls after the type I nonRL procedure is
#' applied to the data set, \code{dat}. Cases are denoted as red dots and controls as black circles in the plot.}
#' \item{init.cases}{The initial cases in the data set, \code{dat}. Marked with red crosses in the plot of the points.}
#' \item{gen.points,ref.points}{Both are \code{NULL} for this function, as initial set of points, \code{dat}, are provided
#' for the non-RL procedure.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of cases and controls.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rnonRLII}}, \code{\link{rnonRLIII}}, \code{\link{rnonRLIV}}, and \code{\link{rnonRL}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-40;  #try also n<-20; n<-100;
#' #data generation
#' dat<-cbind(runif(n,0,1),runif(n,0,1))
#'
#' prop<-.5; #try also .25, .75
#' rho<- .3
#' knn<-3 #try 2 or 5
#'
#' Xdat<-rnonRLI(dat,prop,knn,rho,poisson=FALSE,init=FALSE) 
#' #labeled data try also poisson=TRUE or init=FALSE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #normal original data
#' n<-40;  #try also n<-20; n<-100;
#' #data generation
#' dat<-cbind(rnorm(n,0,1),rnorm(n,0,1))
#'
#' prop<-.50; #try also .25, .75
#' rho<- .3
#' knn<-5 #try 2 or 3
#'
#' Xdat<-rnonRLI(dat,prop,knn,rho,poisson=FALSE) #labeled data try also poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rnonRLI <- function(dat,prop=.50,k,rho,poisson=FALSE,init.from.cases=TRUE)
{
  n<-nrow(dat)
  
  if (prop<0 | prop>1)
  {stop('prop must be in [0,1]')}
  
  lam <- max(round(n*prop-1,0),0) #to avoid -1 for lam
  
  ifelse(poisson==TRUE, N1<-rpois(1,lam), N1<-lam)
  
  if (k>=n)
  {stop('k must be less than the data size, n')}
  
  if (rho < -prop/(1-prop) | rho > 1)
  {stop(paste('rho must be between ',-prop/(1-prop),' and 1 for infection probabilities to be valid'))}
  
  lab<-rep(0,n)
  init.cases<-NULL  #initial cases
  
  cond<-which(c(N1==0, #1
                N1+1>=n, #2
                N1>0 & N1+1<n #3
  )==TRUE)
  
  ind.list<-1:n
  switch(cond,
         "1" = { lab<-lab },
         "2" = { 
           ind0<-sample(ind.list,1) #initial case
           ind.cases<-ind0
           init.cases<- dat[ind0,]
           
           lab<-rep(1,n) },
         "3" = { 
           ind0<-sample(ind.list,1) #initial case
           
           ind.cases<-ind.init.cases<-ind0
           ind.ctrl<-ind.list[-ind0]
           lab[ind0]<-1
           
           ipd<-ipd.mat(dat)
           
           pr<-c(prop+((1-prop)*rho)/(1:k)) #case probability for k NNs
           
           cnt<-1
           while ( cnt <= N1+1 )
           {
             ifelse(init.from.cases==TRUE,ind0<-sample(ind.cases,1),ind0<-sample(ind.ctrl,1)) #initial case at each round
             ind.init.cases<-unique(c(ind.init.cases,ind0)) # the set of initial cases
             ind.ctrl<-ind.list[lab==0]  #index of the controls
             
             lab.can<- rbinom(k,1,prob=pr) #to have cases as 1 and controls as 0
             ind.can<-ind.ctrl[order(ipd[ind0,ind.ctrl])[1:k]]
             lab[ind.can]<-lab.can
             cnt<-sum(lab==1)
             ind.cases<-ind.list[lab==1]
             
           }
           init.cases<-dat[ind.init.cases,]
           if (cnt>N1+1) 
           {relab<-sample(ind.cases[-ind.init.cases],cnt-1-N1)
           lab[relab]<-0
           }
         }
  )
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  n1<-sum(lab==1); n0<-n-n1
  npts<-c(n1,n0)
  names(npts)<-c("n1","n0")
  
  pname <-"parameters"
  prop<-n1/n
  param<-c(prop,rho,k)
  names(param)<-c("proportion of cases","adjusting parameter for NN infection probabilities","Number of NNs")
  typ<-paste("Type I non-RL pattern with ",n1, " cases and ", n0," controls with ",pname, 
             ": proportion of cases = ",param[1],", adj param for NN infection prob = ",param[2],
             " and number of NNs = ",param[3],sep="")
  
  txt<-"type I non-RL pattern (for disease clustering)"
  
  rparam<-sub('^(-)?0[.]', '\\1.', round(param,2)) # to write decimal 0.1 as .1
  main.txt<-paste("Type I Non-RL Pattern with Parameters\n ","infected prop=",rparam[1],", adj param for infect prob=",rparam[2],", # of NNs=",param[3],sep="")
  
  res<-list(
    pat.type="cc", #cc for case-control patterns for RL or non-RL of the given data points
    type=typ,
    parameters=param,
    dat.points=dat,
    lab=lab, #labels of the data points
    init.cases=init.cases, #initial cases
    gen.points=NULL,
    ref.points=NULL,
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Type II Non-Random Labeling of a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Given the set of \eqn{n} points, \code{dat}, in a region, this function assigns \eqn{n_1=}\code{round(n*ult.prop,0)} of them as cases,
#' and the rest as controls with first selecting \eqn{k_0=}\code{round(n*init.prop,0)} as cases initially, then selecting
#' a contagious case and then assigning the label case to the remaining points with infection probabilities 
#' inversely proportional to their position among the \code{k}NNs.
#' 
#' The initial and ultimate number of cases will be \eqn{k_0} and \eqn{n_1} on the average if the argument \code{poisson=TRUE}
#' (i.e., \eqn{k_0=}\code{rpois(1,round(n*init.prop,0)}) and \eqn{n_1=}\code{rpois(1,round(n*ult.prop,0))} ), otherwise
#' they will be exactly equal to \eqn{n_1=}\code{round(n*ult.prop,0)} and \eqn{k_0=}\code{round(n*init.prop,0)}.
#' More specifically, let \eqn{z_1,\ldots,z_{k_0}} be the initial cases. Then one of the cases is selected as a
#' contagious case, say \eqn{z_j} and then its \code{k}NNs (among the non-cases) are found.
#' Then label these \code{k}NN non-case points as cases with infection probabilities \code{prob} equal to the value
#' of the \code{rho*(1/(1:k))^pow} values at these points, where \code{rho} is a scaling parameter for
#' the infection probabilities and \code{pow} is a parameter in the power adjusting the \code{k}NN dependence.
#' We stop when we first exceed \eqn{n_1} cases. \code{rho} has to be in \eqn{(0,1)} for \code{prob} to be a \code{vector} of probabilities,
#' and for a given \code{rho}, \code{pow} must be \eqn{>  \ln(rho)/\ln(k)}.
#' If \code{rand.init=TRUE}, first \eqn{k_0} entries are chosen as the initial cases in the data set,
#' \code{dat}, otherwise, \eqn{k_0} initial cases are selected randomly among the data points.
#' 
#' Algorithmically, first all dat points are treated as non-cases (i.e. controls or healthy subjects).
#' Then the function follows the following steps for labeling of the points: 
#' 
#' step 0: \eqn{n_1} is generated randomly from a Poisson distribution with \code{mean = round(n*ult.prop,0)}, so that the 
#' average number of ultimate cases will be \code{round(n*ult.prop,0)} if the argument \code{poisson=TRUE}, else \eqn{n_1=}\code{round(n*ult.prop,0)}.
#' And \eqn{k_0} is generated randomly from a Poisson distribution with \code{mean = round(n*init.prop,0)}, so that the 
#' average number of initial cases will be \code{round(n*init.prop,0)} if the argument \code{poisson=TRUE}, else \eqn{k_0=}\code{round(n*init.prop,0)}.
#' 
#' step 1: Initially, \eqn{k_0} many points from dat are selected as cases.
#' The selection of initial cases are determined based on the argument \code{rand.init} (with default=\code{TRUE})
#' where if \code{rand.init=TRUE} then the initial cases are selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first \eqn{k_0} entries in the data set, \code{dat}, are selected as the cases.
#' 
#' step 2: Then it selects a contagious case among the cases, and randomly labels its \code{k} control NNs as cases with
#' decreasing infection probabilities \code{prob=rho*(1/(1:k))^pow}. See the description for the details
#' of the parameters in the \code{prob}.
#' 
#' step 3: The procedure ends when number of cases \eqn{n_c} exceeds \eqn{n_1}, and \eqn{n_c-n_1} of the cases (other than the
#' initial cases) are randomly selected and relabeled as controls, i.e. 0s,
#' so that the number of cases is exactly \eqn{n_1}.
#' 
#' Note that the infection probabilities of the \code{k}NNs of each initial case increase
#' with increasing rho; and probability of infection decreases as further NNs are considered from 
#' a contagious case (i.e. as \code{k} increases in the \code{k}NNs).
#' 
#' See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail where type II non-RL pattern is the 
#' case 2 of non-RL pattern considered in Section 6 with \eqn{n_1} is
#' fixed as a parameter rather than being generated from a Poisson distribution and \code{pow=1}.
#' 
#' Although the non-RL pattern is described for the case-control setting, it can be adapted for any two-class
#' setting when it is appropriate to treat one of the classes as cases or one of the classes behave like cases
#' and other class as controls.
#'
#' @param dat A set of points the non-RL procedure is applied to obtain cases and controls randomly in the 
#' type II fashion (see the description).
#' @param k An integer representing the number of NNs considered for each contagious case, i.e., 
#' \code{k}NNs of each contagious case are candidates to be infected to become cases.
#' @param rho A scaling parameter for the probabilities of labeling the points as cases
#' (see the description).
#' @param pow A parameter in the power adjusting the \code{k}NN dependence in the probabilities of labeling the
#' points as cases (see the description).
#' @param init.prop A real number between 0 and 1 representing the initial proportion of cases in the data set,
#' \code{dat}. The selection of the initial cases depends on the parameter \code{rand.init} (see the description).
#' @param ult.prop A real number between 0 and 1 representing the ultimate proportion of cases in the data set,
#' \code{dat} after the non-RL assignment.
#' @param rand.init A logical argument (default is \code{TRUE}) to determine the choice of the initial cases in the data set, \code{dat}.
#' If \code{rand.init=TRUE} then the initial cases are selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first \code{init.prop*n} entries in the data set, \code{dat}, are labeled as the cases.
#' @param poisson A logical argument (default is \code{FALSE}) to determine whether the number of initial and ultimate
#' cases, \eqn{k_0} and \eqn{n_1}, will be random or fixed. If \code{poisson=TRUE} then the \eqn{k_0} and \eqn{n_1} are from a Poisson distribution,
#' \eqn{k_0=}\code{rpois(1,round(n*init.prop,0)}) and \eqn{n_1=}\code{rpois(1,round(n*ult.prop,0))}
#' otherwise they are fixed, \eqn{k_0=}\code{round(n*init.prop,0)} and \eqn{n_1=}\code{round(n*ult.prop,0)}.
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{="cc"} for the case-control patterns for RL or non-RL of the given data points, \code{dat}}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{Number of NNs, \code{k}, a scaling parameter for the infection probabilities of \code{k}NNs, rho,
#' a parameter in the power adjusting the \code{k}NN dependence of the infection probabilities, initial proportion
#' of cases, \code{init.prop}, and the ultimate proportion of cases, \code{ult.prop}.}
#' \item{dat.points}{The set of points non-RL procedure is applied to obtain cases and controls randomly in the 
#' type II fashion}
#' \item{lab}{The labels of the points as 1 for cases and 0 for controls after the type II nonRL procedure is
#' applied to the data set, \code{dat}. Cases are denoted as red dots and controls as black circles in the plot.}
#' \item{init.cases}{The initial cases in the data set, \code{dat}. Denoted as red crosses in the plot of the points.}
#' \item{cont.cases}{The contagious cases in the data set, \code{dat}. Denoted as blue points in the plot of the points.}
#' \item{gen.points,ref.points}{Both are \code{NULL} for this function, as initial set of points, \code{dat}, are provided
#' for the non-RL procedure.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of cases and controls.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rnonRLI}}, \code{\link{rnonRLIII}}, \code{\link{rnonRLIV}}, and \code{\link{rnonRL}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-40;  #try also n<-20; n<-100;
#' #data generation
#' dat<-cbind(runif(n,0,1),runif(n,0,1))
#'
#' rho<-.8
#' pow<-2
#' knn<-5 #try 2 or 3
#' ip<-.3 #initial proportion
#' up<-.5 #ultimate proportion
#'
#' Xdat<-rnonRLII(dat,knn,rho,pow,ip,up,poisson=FALSE) #labeled data, try poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #normal original data
#' n<-40;  #try also n<-20; n<-100;
#' #data generation
#' dat<-cbind(rnorm(n,0,1),rnorm(n,0,1))
#'
#' rho<-0.8
#' pow<-2
#' knn<-5 #try 2 or 3
#' ip<-.3 #initial proportion
#' up<-.5 #ultimate proportion
#'
#' Xdat<-rnonRLII(dat,knn,rho,pow,ip,up,poisson=FALSE) #labeled data, try poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rnonRLII <- function(dat,k,rho,pow,init.prop,ult.prop,rand.init=TRUE,poisson=FALSE)
{
  n<-nrow(dat)
  if (k>=n)
  {stop('k must be less than the data size, n')}
  
  if (init.prop > ult.prop )
  {stop('initial proportion must be smaller than or equal to ultimate proportion')}
  
  if (init.prop<0 | ult.prop>1 )
  {stop("both init.prop and ult.prop must be in [0,1]")}
  
  if (rho < 0 | rho > 1)
  {stop(paste('rho must be between 0 and 1 for infection probabilities to be valid'))}
  
  lam.k0<-round(n*init.prop,0) #number of initial cases (on the average)
  lam.N1<-round(n*(ult.prop-init.prop),0) #number of new cases (on the average)
  
  ifelse(poisson==TRUE, {k0<-rpois(1,lam.k0);  N1<-rpois(1,lam.N1)},{k0<-lam.k0; N1<-lam.N1})
  
  lab<-rep(0,n)
  init.cases<-cont.cases<-NULL #initial and contagious cases
  
  cond<-which(c(k0==0 & N1==0, #1
                k0==0 & N1!=0 , #2
                k0!=0 & N1==0, #3
                k0!=0 & N1!=0 & N1+k0>=n, #4
                k0!=0 & N1>0 & N1+k0<n #5
  )==TRUE)
  
  ind.list<-1:n #index list
  
  switch(cond,
         "1" = { lab<-lab
         },
         "2" = { 
           N1<-min(N1,n)
           ind.cases<-sample(ind.list,N1)
           lab[ind.cases]<-1
         },
         "3" = {
           k0<-min(k0,n)
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           lab[ind0]<-1
           
           init.cases<-dat[ind0,]
           indc<-sample(ind0,1) #initial contagious case
           ind.cc<-indc #indices of contagious cases
           
           cont.cases<-dat[indc,]
         },
         "4" = { 
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           
           init.cases<-dat[ind0,]
           indc<-sample(ind0,1) #initial contagious case
           ind.cc<-indc #indices of contagious cases
           
           org.ind.list<-ind.list[-ind0] #original index \code{list} 
           cont.cases<-dat[indc,]
           lab<-rep(1,n)
         },
         "5" = { 
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           lab[ind0]<-1
           
           init.cases<-dat[ind0,]
           indc<-sample(ind0,1) #initial contagious case
           ind.cc<-indc #indices of contagious cases
           
           org.ind.list<-ind.list[-ind0] #original index \code{list} 
           cont.cases<-dat[indc,]
           
           ipd<-ipd.mat(dat)
           pr<-rho*(1/(1:k))^pow
           cnt<-k0
           while ( cnt <= N1+k0 )
           { 
             can.ord<-as.numeric(na.omit(order(ipd[indc,-ind.cases])[1:k]))
             ind.can0<-can.ord 
             ind.can<-org.ind.list[ind.can0]
             #labels of the candidates
             kn<-length(can.ord)
             lab.can<- rbinom(kn,1,prob=pr)
             lab[ind.can]<-lab.can
             cnt<-sum(lab==1)
             ind.cases<-ind.list[lab==1]
             org.ind.list<-ind.list[-ind.cases]
             indc<-sample(ind.cases,1) #contagious case
             ind.cc<-c(ind.cc,indc)
             cont.cases<-rbind(cont.cases,dat[indc,])
           }
           if (cnt>N1+k0) 
           {relab<-sample(ind.cases[-ind0],cnt-k0-N1)
           lab[relab]<-0
           }
         }
  )
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  n1<-sum(lab==1); n0<-n-n1
  npts<-c(n1,n0)
  names(npts)<-c("n1","n0")
  
  pname <-"parameters"
  init.prop<-k0/n; ult.prop<-n1/n
  param<-c(k,rho,pow,init.prop,ult.prop)
  names(param)<-c("Number of NNs","infection prob of kNNs","initial proportion of cases",
                  "ultimate proportion of cases")
  typ<-paste("Type II non-RL pattern with ",n1, " cases and ", n0," controls with ",pname, 
             " number of NNs = ",param[1],", scaling parameter = ",param[2],", power for NN dependence = ",param[3],
             ", initial proportion of cases = ",param[5]," and ultimate proportion of cases = ",param[5],sep="")
  
  rparam<-sub('^(-)?0[.]', '\\1.', round(param,2)) # to write decimal 0.1 as .1
  txt<-"type II non-RL pattern (for disease clustering)"
  main.txt<-paste("Type II Non-RL Pattern with Parameters\n ","  # of NNs=",param[1],", scaling param=",rparam[2],
                  ", power for NN dep=",rparam[3],", init prop=",rparam[4],", ult prop=",rparam[5],sep="")
  res<-list(
    pat.type="cc", #cc for case-control patterns for RL or non-RL of the given data points
    type=typ,
    parameters=param,
    dat.points=dat,
    lab=lab, #labels of the data points
    init.cases=init.cases, #initial cases
    cont.cases=cont.cases, #contagious cases
    gen.points=NULL, #generated points associated with Y points
    ref.points=NULL, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Type III Non-Random Labeling of a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Given the set of \eqn{n} points, \code{dat}, in a region, this function assigns \eqn{n_1=}\code{round(n*prop,0)} of them as cases,
#' and the rest as controls with first selecting a point, \eqn{Z_i}, as a case and assigning the
#' label case to the remaining points with infection probabilities \eqn{prob=rho (1-d_{ij}/d_{\max})^{pow}} where \eqn{d_{ij}} is the
#' distance from \eqn{Z_j} to \eqn{Z_i} for \eqn{j \ne i}, \eqn{d_{\max}} is the maximum of  \eqn{d_{ij}}  values, \code{rho} is a scaling parameter for
#' the infection probabilities and \code{pow} is a parameter in the power adjusting the distance dependence.
#' The number of cases will be \eqn{n_1} on the average if the argument \code{poisson=TRUE}
#' (i.e., \eqn{n_1=}\code{rpois(1,round(n*prop,0))} ), otherwise \eqn{n_1=}\code{round(n*prop,0)}.
#' We stop when we first exceed \eqn{n_1} cases. \code{rho} has to be positive for \code{prob} to be a \code{vector} of probabilities,
#' and for a given \code{rho}, \code{pow} must be \eqn{> - \ln(rho)/\ln(1-d_{ij}/d_{\max})},
#' also, when \code{pow} is given, \code{rho} must be \eqn{< (1-d_{ij}/d_{\max})^{-pow}}.
#' If \code{rand.init=TRUE}, initial case is selected randomly among the data points,
#' otherwise first entry is chosen as the case (or case is recorded as the first entry) in the data set, \code{dat}. 
#' 
#' Algorithmically, first all dat points are treated as non-cases (i.e. controls or healthy subjects).
#' Then the function follows the following steps for labeling of the points: 
#' 
#' step 0: \eqn{n_1} is generated randomly from a Poisson distribution with \code{mean = round(n*prop,0)}, so that the 
#' average number of cases will be round(n*prop,0) if the argument \code{poisson=TRUE}, else \eqn{n_1=}\code{round(n*prop,0)}.
#' 
#' step 1: Initially, one point from dat is selected as a case.
#' The selection of initial case is determined based on the argument \code{rand.init} (with default=\code{TRUE})
#' where if \code{rand.init=TRUE} then the initial case is selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first entry in the data set, \code{dat}, is selected as the case.
#' 
#' step 2: Then it assigns the label case to the remaining points
#' with infection probabilities \eqn{prob=rho (1-d_{ij}/d_{\max})^{pow}}, see the description for the details
#' of the parameters in the \code{prob}.
#' 
#' step 3: The procedure ends when number of cases \eqn{n_c} exceeds \eqn{n_1}, and \eqn{n_c-n_1} of the cases (other than the
#' initial contagious case) are randomly selected and relabeled as controls, i.e. 0s,
#' so that the number of cases is exactly \eqn{n_1}.
#' 
#' In the output cases are labeled as 1 and controls as 0, and initial contagious case is marked with a red cross
#' in the plot of the pattern.
#' Note that the infection probabilities of the points is inversely proportional to their distances to the
#' initial case and increase with increasing \code{rho}. 
#' This function might take a long time for certain choices of the arguments. For example, if \code{pow} is taken to be
#' too large, the infection probabilities would be too small, and case assignment will take a rather long time. 
#' 
#' See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail where type III non-RL pattern is the 
#' case 3 of non-RL pattern considered in Section 6 with \eqn{n_1} is
#' fixed as a parameter rather than being generated from a Poisson distribution and \eqn{k_{den}=1} and pow
#' is represented as \eqn{k_{pow}}.
#' 
#' Although the non-RL pattern is described for the case-control setting, it can be adapted for any two-class
#' setting when it is appropriate to treat one of the classes as cases or one of the classes behave like cases
#' and other class as controls.
#' 
#' @param dat A set of points the non-RL procedure is applied to obtain cases and controls randomly in the 
#' type III fashion (see the description).
#' @param prop A real number between 0 and 1 (inclusive) representing the proportion of new cases (on the average)
#' infected by the initial case, i.e., number of newly infected cases (in addition to the first case) is Poisson
#' with \code{mean=round(n*prop)} where \eqn{n} is the number of points in \code{dat}, if the argument \code{poisson=TRUE}, else it is \code{round(n*prop)}.
#' @param rho A scaling parameter for the probabilities of labeling the points as cases
#' (see the description).
#' @param pow A parameter in the power adjusting the distance dependence in the probabilities of labeling the
#' points as cases (see the description).
#' @param rand.init A logical argument (default is \code{TRUE}) to determine the choice of the initial case in the data set, \code{dat}.
#' If \code{rand.init=TRUE} then the initial case is selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first entry in the data set, \code{dat}, is labeled as the initial case.
#' @param poisson A logical argument (default is \code{FALSE}) to determine whether the number of cases \eqn{n_1},
#' will be random or fixed. If \code{poisson=TRUE} then the \eqn{n_1} is from a Poisson distribution, 
#' \eqn{n_1=}\code{rpois(1,round(n*prop,0))} otherwise it is fixed, \eqn{n_1=}\code{round(n*prop,0)}.
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{="cc"} for the case-control patterns for RL or non-RL of the given data points, \code{dat}}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{rho and pow, where \code{rho} is the scalign parameter and \code{pow} is the parameter in the power
#' adjusting the distance dependence in probabilities of labeling the points as cases.}
#' \item{dat.points}{The set of points non-RL procedure is applied to obtain cases and controls randomly in the 
#' type III fashion}
#' \item{lab}{The labels of the points as 1 for cases and 0 for controls after the type III nonRL procedure is
#' applied to the data set, \code{dat}. Cases are denoted as red dots and controls as black circles in the plot.}
#' \item{init.cases}{The initial case in the data set, \code{dat}. Marked with a red cross in the plot of the points.}
#' \item{cont.cases}{The contagious cases in the data set, \code{dat}. Denoted as blue points in the plot of the points.}
#' \item{gen.points,ref.points}{Both are \code{NULL} for this function, as initial set of points, \code{dat}, are provided
#' for the non-RL procedure.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of cases and controls.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rnonRLI}}, \code{\link{rnonRLII}}, \code{\link{rnonRLIV}}, and \code{\link{rnonRL}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-40;  #try also n<-20; n<-100;
#' prop<- .5; #try also .25, .75
#' #data generation
#' dat<-cbind(runif(n,0,1),runif(n,0,1))
#'
#' rho<-.8
#' pow<-2
#'
#' Xdat<-rnonRLIII(dat,prop,rho,pow,poisson=FALSE) #labeled data, try also poisson=TRUE
#'
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #normal original data
#' n<-40;  #try also n<-20; n<-100;
#' dat<-cbind(rnorm(n,0,1),rnorm(n,0,1))
#'
#' prop<- .5; #try also .25, .75
#' rho<-.8
#' pow<-2
#'
#' Xdat<-rnonRLIII(dat,prop,rho,pow,poisson=FALSE) #labeled data, try also poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rnonRLIII <- function(dat,prop,rho,pow,rand.init=TRUE,poisson=FALSE)
{
  n<-nrow(dat)
  
  if (prop<0 | prop>1)
  {stop('prop must be in [0,1]')}
  
  lam <- max(round(n*prop-1,0),0) #to avoid -1 for lam
  
  ifelse(poisson==TRUE, N1<-rpois(1,lam), N1<-lam)
  
  lab<-rep(0,n)
  init.cases<-NULL  #initial cases
  
  cond<-which(c(N1==0, #1
                N1+1>=n, #2
                N1>0 & N1+1<n #3
  )==TRUE)
  
  switch(cond,
         "1" = { lab<-lab },
         "2" = { 
           ind.list<-1:n
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,1), ind0<-1) #initial case
           ind.cases<-ind0
           init.cases<- dat[ind0,]
           
           lab<-rep(1,n) },
         "3" = { 
           ind.list<-1:n
           
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,1), ind0<-1) #initial case
           
           ind.cases<-ind0
           lab[ind0]<-1
           init.cases<- dat[ind0,]
           
           ipd<-ipd.mat(dat)
           max.dis<-max(ipd)
           
           cnt<-1
           while ( cnt <= N1+1 )
           {
             ind.can<-ind.list[-ind.cases]  #index of the candidates
             num.left<-n-cnt
             #labels of the candidates
             pr<- rho*((1-ipd[ind0,-ind.cases]/max.dis)^pow)
             lab.can<- rbinom(num.left,1,prob=pr) #to have cases as 1 and controls as 0
             lab[ind.can]<-lab.can
             cnt<-sum(lab==1)
             ind.cases<-ind.list[lab==1]
           }
           if (cnt>N1+1) 
           {relab<-sample(ind.cases[-ind0],cnt-1-N1)
           lab[relab]<-0
           }
         }
  )
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  n1<-sum(lab==1); n0<-n-n1 #n1 = # of cases, n0=# of controls
  npts<-c(n1,n0)
  names(npts)<-c("n1","n0")
  
  pname <-"parameters"
  prop<-n1/n
  param<-c(prop,rho,pow)
  names(param)<-c("scaling parameter","power for distance dependence")
  typ<-paste("Type III non-RL pattern with ",n1, " cases and ", n0," controls with proportion of cases =", param[1],
             ", scaling parameter =", param[2], " and power for distance dependence = ",param[3],sep="")
  
  rparam<-sub('^(-)?0[.]', '\\1.', round(param,2)) # to write decimal 0.1 as .1 
  
  txt<-"type III non-RL pattern (for disease clustering)"
  main.txt<-paste("Type III Non-RL Pattern with\n prop of cases=",rparam[1],", scaling param=",rparam[2],", power for dist dep=",rparam[3],sep="")
  
  res<-list(
    pat.type="cc", #cc for case-control patterns for RL or non-RL of the given data points
    type=typ,
    parameters=param,
    dat.points=dat,
    lab=lab, #labels of the data points
    init.cases=init.cases, #initial case
    gen.points=NULL,
    ref.points=NULL,
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title pdf of the Bivariate Normal Distribution
#'
#' @description
#' Computes the value of the probability density function (i.e. density) of the bivariate normal distribution
#' at the specified point \code{X}, with mean \code{mu} and standard deviations of the first and second components being \eqn{s_1}
#' and \eqn{s_2} (denoted as \code{s1} and \code{s2} in the arguments of the function, respectively) 
#' and correlation between them being \code{rho} (i.e., the covariance matrix is \eqn{\Sigma=S} where \eqn{S_{11}=s_1^2},
#' \eqn{S_{22}=s_2^2}, \eqn{S_{12}=S_{21}=s_1 s_2 rho}). 
#' 
#' @param X A set of 2D points of size \eqn{n} (i.e an \eqn{n \times 2} matrix or array) at which the density of the bivariate normal distribution
#' is to be computed.
#' @param mu A \eqn{1 \times 2} \code{vector} of real numbers representing the mean of the bivariate normal distribution,
#' default=\eqn{(0,0)}. 
#' @param s1,s2 The standard deviations of the first and second components of the bivariate normal distribution,
#' with default is \code{1} for both
#' @param rho The correlation between the first and second components of the bivariate normal distribution
#' with default=0.
#' 
#' @return 
#' The value of the probability density function (i.e. density) of the bivariate normal distribution
#' at the specified point \code{X}, with mean \code{mu} and standard deviations of the first and second components being \eqn{s_1}
#' and \eqn{s_2} and correlation between them being \code{rho}.
#'
#' @seealso \link[MASS]{mvrnorm}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' mu<-c(0,0)
#' s1<-1
#' s2<-1
#' rho<-.5
#'
#' n<-5
#' Xp<-cbind(runif(n),runif(n))
#' bvnorm.pdf(Xp,mu,s1,s2,rho)
#' @export 
bvnorm.pdf <- function(X,mu=c(0,0),s1=1,s2=1,rho=0)
{
  X<-matrix(X,ncol=2)
  x<-X[,1]
  y<-X[,2]
  mu1<-mu[1]
  mu2<-mu[2]
  z<-((x-mu1)/s1)^2-2*rho*(x-mu1)*(y-mu2)/(s1*s2)+((y-mu2)/s2)^2
  val<-(2*pi*s1*s2*sqrt(1-rho^2))^(-1)*exp(-z/(1-rho^2))
  val
} #end for the function
#'

#################################################################

#' @title Type IV Non-Random Labeling of a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Given the set of \eqn{n} points, \code{dat}, in a region, this function assigns \eqn{n_1=}\code{round(n*ult.prop,0)} of them as cases,
#' and the rest as controls with first selecting \eqn{k_0=}\code{round(n*init.prop,0)} as cases initially and assigning the
#' label case to the remaining points with infection probabilities equal to the scaled bivariate normal density values
#' at those points.
#' The initial and ultimate number of cases will be \eqn{k_0} and \eqn{n_1} on the average if the argument \code{poisson=TRUE}
#' (i.e., \eqn{k_0=}\code{rpois(1,round(n*init.prop,0)}) and \eqn{n_1=}\code{rpois(1,round(n*ult.prop,0))} ), otherwise
#' they will be exactly equal to \eqn{n_1=}\code{round(n*ult.prop,0)} and \eqn{k_0=}\code{round(n*init.prop,0)}.
#' More specifically, let \eqn{z_1,\ldots,z_{k_0}} be the initial cases and for \eqn{j=1,2,\ldots,k_0}
#' let \eqn{\phi_{G,j}(z_i)} be the value of the pdf of the \eqn{BVN(z_j,s_1,s_2,rho)}, which is the bivariate normal 
#' distribution mean=z_j and standard deviations of the first and second components being \eqn{s_1}
#' and \eqn{s_2} (denoted as \code{s1} and \code{s2} as arguments of the function) and 
#' correlation between them being \eqn{\rho} (denoted as \code{rho} as an argument of the function)
#' (i.e., the covariance matrix is \eqn{\Sigma=S} where \eqn{S_{11}=s_1^2},
#' \eqn{S_{22}=s_2^2}, \eqn{S_{12}=S_{21}=s_1 s_2 \rho}). Add these pdf values as
#' \eqn{p_j=\sum_{j=1}^{k_0} \phi_{G,j}(z_i)} for each \eqn{i=1,2,\ldots,n} and find \eqn{p_{\max}=\max p_j}. 
#' Then label the points (other than the initial cases) as cases with infection probabilities \code{prob} equal to the value
#' of the \eqn{p_j/p_{\max}} values at these points. 
#' We stop when we first exceed \eqn{n_1} cases. \eqn{\rho} has to be in (-1,1) for \code{prob} to be
#' a valid probability and \eqn{s_1} and \eqn{s_2} must be positive (actually these are required for the BVN density
#' to be nondegenerately defined).
#' If \code{rand.init=TRUE}, first \eqn{k_0} entries are chosen as the initial cases in the data set,
#' \code{dat}, otherwise, \eqn{k_0} initial cases are selected randomly among the data points.
#' 
#' Algorithmically, first all dat points are treated as non-cases (i.e. controls or healthy subjects).
#' Then the function follows the following steps for labeling of the points: 
#' 
#' step 0: \eqn{n_1} is generated randomly from a Poisson distribution with \code{mean = round(n*ult.prop,0)}, so that the 
#' average number of ultimate cases will be \code{round(n*ult.prop,0)} if the argument \code{poisson=TRUE}, else \eqn{n_1=}\code{round(n*ult.prop,0)}.
#' And \eqn{k_0} is generated randomly from a Poisson distribution with \code{mean = round(n*init.prop,0)}, so that the 
#' average number of initial cases will be round(n*init.prop,0) if the argument \code{poisson=TRUE}, else \eqn{k_0=}\code{round(n*init.prop,0)}. 
#' 
#' step 1: Initially, \eqn{k_0} many points from dat are selected as cases.
#' The selection of initial cases are determined based on the argument \code{rand.init} (with default=\code{TRUE})
#' where if \code{rand.init=TRUE} then the initial cases are selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first \eqn{k_0} entries in the data set, \code{dat}, are selected as the cases.
#' 
#' step 2: Then it assigns the label case to the remaining points
#' with infection probabilities \eqn{prob=\sum_{j=1}^{k_0} \phi_{G,j}(z_i)/p_{\max}},
#' which is the sum of the BVN densities scaled by the maximum of such sums.
#' See the description for the details of the parameters in the \code{prob}.
#' 
#' step 3: The procedure ends when number of cases \eqn{n_c} exceed \eqn{n_1}, and \eqn{n_c-n_1} of the cases (other than the initial
#' cases) are randomly selected and relabeled as controls, i.e. 0s, so that the number of cases is
#' exactly \eqn{n_1}.
#' 
#' In the output cases are labeled as 1 and controls as 0, and initial contagious case is marked with a red cross
#' in the plot of the pattern.
#' 
#' See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} for more detail where type IV non-RL pattern is the 
#' case 4 of non-RL pattern considered in Section 6 with \eqn{n_1} and \eqn{k_0} are
#' fixed as parameters and \code{rho} is represented as \eqn{k_{pow}} and \eqn{rho/k_{den}=1} in the article.
#' 
#' Although the non-RL pattern is described for the case-control setting, it can be adapted for any two-class
#' setting when it is appropriate to treat one of the classes as cases or one of the classes behave like cases
#' and other class as controls.
#' 
#' @param dat A set of points the non-RL procedure is applied to obtain cases and controls randomly in the 
#' type IV fashion (see the description).
#' @param init.prop A real number between 0 and 1 representing the initial proportion of cases in the data set,
#' \code{dat}. The selection of the initial cases depends on the parameter \code{rand.init} and the number of initial cases
#' depends on the parameter poisson (see the description).
#' @param ult.prop A real number between 0 and 1 representing the ultimate proportion of cases in the data set,
#' \code{dat} after the non-RL assignment. The number of ultimate cases depends on the parameter poisson
#' (see the description).
#' @param s1,s2 Positive real numbers representing the standard deviations of the first and second components
#' of the bivariate normal distribution.
#' @param rho A real number between -1 and 1 representing the correlation between the first and second components
#' of the bivariate normal distribution.
#' @param rand.init A logical argument (default is \code{TRUE}) to determine the choice of the initial case in the data set, \code{dat}.
#' If \code{rand.init=TRUE} then the initial case is selected randomly from the data points, and if \code{rand.init=}
#' \code{FALSE}, the first \eqn{k_0} entries in the data set, \code{dat}, is labeled as the initial case.
#' @param poisson A logical argument (default is \code{FALSE}) to determine whether the number of initial and ultimate
#' cases, \eqn{k_0} and \eqn{n_1}, will be random or fixed. If \code{poisson=TRUE} then the \eqn{k_0} and \eqn{n_1} are from a Poisson distribution,
#' \eqn{k_0=}\code{rpois(1,round(n*init.prop,0))} and \eqn{n_1=}\code{rpois(1,round(n*ult.prop,0))}
#' otherwise they are fixed, \eqn{k_0=}\code{round(n*init.prop,0)} and \eqn{n_1=}\code{round(n*ult.prop,0)}.
#'  
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{="cc"} for the case-control patterns for RL or non-RL of the given data points, \code{dat}}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{initial and ultimate proportion of cases after the non-RL procedure is applied to the data,
#' \code{s1}, \code{s2} and \code{rho} which are standard deviations and the correlation for the components of
#' the bivariate normal distribution.}
#' \item{dat.points}{The set of points non-RL procedure is applied to obtain cases and controls randomly in the 
#' type IV fashion}
#' \item{lab}{The labels of the points as 1 for cases and 0 for controls after the type IV nonRL procedure is
#' applied to the data set, \code{dat}. Cases are denoted as red dots and controls as black circles in the plot.}
#' \item{init.cases}{The initial cases in the data set, \code{dat}. Marked with red crosses in the plot of the points.}
#' \item{gen.points,ref.points}{Both are \code{NULL} for this function, as initial set of points, \code{dat}, are provided
#' for the non-RL procedure.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of cases and controls.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rnonRLI}}, \code{\link{rnonRLII}}, \code{\link{rnonRLIII}}, and \code{\link{rnonRL}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-40;  #try also n<-20; n<-100;
#' ult<-.5; #try also .25, .75
#' #data generation
#' dat<-cbind(runif(n,0,1),runif(n,0,1))
#'
#' int<-.1
#' s1<-s2<-.4
#' rho<- .1
#'
#' Xdat<-rnonRLIV(dat,int,ult,s1,s2,rho,poisson=FALSE) #labeled data, try also with poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #normal original data
#' n<-40;  #try also n<-20; n<-100;
#' dat<-cbind(rnorm(n,0,1),rnorm(n,0,1))
#' ult<-.5; #try also .25, .75
#'
#' int<-.1
#' s1<-s2<-.4
#' rho<-0.1
#'
#' Xdat<-rnonRLIV(dat,int,ult,s1,s2,rho,poisson=FALSE) #labeled data, try also with poisson=TRUE
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'  
#' @export 
rnonRLIV <- function(dat,init.prop,ult.prop,s1,s2,rho,rand.init=TRUE,poisson=FALSE)
{
  n<-nrow(dat)
  
  if (init.prop > ult.prop )
  {stop('initial proportion must be smaller than or equal to ultimate proportion')}
  
  if (init.prop<0 | ult.prop>1 )
  {stop("both init.prop and ult.prop must be in [0,1]")}
  
  if (rho <= -1 | rho >= 1)
  {stop(paste('rho must be between -1 and 1 for BVN pdf to be nondegenerate (i.e for infection probabilities to be valid)'))}
  
  if (min(s1,s2) <= 0)
  {stop(paste('both s1 and s2 must be positive for BVN pdf to be nondegenerate (i.e for infection probabilities to be valid)'))}
  
  lam.k0<-round(n*init.prop,0) #number of initial cases (on the average)
  lam.N1<-round(n*(ult.prop-init.prop),0) #number of new cases (on the average)
  
  ifelse(poisson==TRUE, {k0<-rpois(1,lam.k0);  N1<-rpois(1,lam.N1)},{k0<-lam.k0; N1<-lam.N1})
  
  lab<-rep(0,n)
  init.cases<-NULL #initial cases
  
  cond<-which(c(k0==0 & N1==0, #1
                k0==0 & N1!=0 , #2
                k0!=0 & N1==0, #3
                k0!=0 & N1!=0 & N1+k0>=n, #4
                k0!=0 & N1>0 & N1+k0<n #5
  )==TRUE)
  
  ind.list<-1:n #index list
  
  switch(cond,
         "1" = { lab<-lab
         },
         "2" = { 
           N1<-min(N1,n)
           ind.cases<-sample(ind.list,N1)
           lab[ind.cases]<-1
         },
         "3" = {
           k0<-min(k0,n)
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           lab[ind0]<-1
           
           init.cases<-dat[ind0,]
         },
         "4" = { 
           ifelse(rand.init==TRUE,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           
           init.cases<-dat[ind0,]
           lab<-rep(1,n)
         },
         "5" = { 
           ifelse(rand.init==TRUE & k0<=n,ind0<-sample(ind.list,k0), ind0<-ind.list[1:k0]) #initial cases
           ind.cases<-ind0
           lab[ind0]<-1
           
           init.cases<-dat[ind0,]
           
           pdf<-vector()
           for (i in 1:k0)
           {
             pdf<- cbind(pdf, bvnorm.pdf(dat,dat[ind0[i],],s1,s2,rho) )
           }
           
           prob<-as.vector(apply(pdf,1,sum))
           max.prob<-max(prob)
           
           cnt<-k0
           while ( cnt <= N1+k0 )
           {
             ind.can<-ind.list[-ind.cases] #index of the candidates
             num.left<-n-cnt
             #labels of the candidates
             max.prob<-max(prob[-ind.cases])
             pr<- prob[-ind.cases] /max.prob
             lab.can<- rbinom(num.left,1,prob=pr) #to have cases as 1 and controls as 0
             lab[ind.can]<-lab.can
             cnt<-sum(lab==1)
             ind.cases<-ind.list[lab==1]
           }
           if (cnt>N1+k0) 
           {relab<-sample(ind.cases[-ind0],cnt-k0-N1)
           lab[relab]<-0
           }
         }
  )
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  n1<-sum(lab==1);  n0<-n-n1 #n1 = # of cases, n0=# of controls
  npts<-c(n1,n0)
  names(npts)<-c("n1","n0")
  
  pname <-"parameters"
  init.prop<-k0/n; ult.prop<-n1/n
  param<-c(init.prop,ult.prop,s1,s2,rho)
  names(param)<-c("initial proportion of cases","ultimate proportion of cases",
                  "sigma1", "sigma2","corr. coeff.")
  
  typ<-paste("Type IV non-RL pattern with ",n1, " cases and ", n0," controls with ",pname, 
             " initial proportion of cases = ",param[1],", ultimate proportion of cases = ",param[2],
             ", sigma1 = ",param[3],", sigma2 = ",param[4],
             " and correlation coefficient = ",param[5],sep="")
  
  rparam<-rep(0,5)
  for (i in 1:5)
  {ifelse(param[i]%%1==0,rparam[i]<-param[i],rparam[i]<-round(param[i],2))}
  
  txt<-"type IV non-RL pattern (for disease clustering)"
  
  main.txt<-bquote(atop("Type IV Non-RL Pattern with Parameters",
                        "init prop" ==.(rparam[1]) *"," ~ "ult prop" ==.(rparam[2]) *"," ~sigma[1]== .(rparam[3]) *","~
                          sigma[2]== .(rparam[4]) *"," ~rho== .(rparam[5]) ))
  
  res<-list(
    pat.type="cc", #cc for case-control patterns for RL or non-RL of the given data points
    type=typ,
    parameters=param,
    dat.points=dat,
    lab=lab, #labels of the data points
    init.cases=init.cases, #initial cases
    gen.points=NULL, #generated points associated with Y points
    ref.points=NULL, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Non-Random Labeling of a Given Set of Points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Given the set of \eqn{n} points, \code{dat}, in a region, this function assigns some of them as cases,
#' and the rest as controls in a non-RL \code{type=type} fashion. 
#' 
#' Type I nonRL pattern assigns \eqn{n_1=}\code{round(n*prop,0)} of the data points as cases,
#' and the rest as controls with first selecting a point, \eqn{Z_i}, as a case and assigning the
#' label case to the remaining points with infection probabilities \code{prob=c(prop+((1-prop)*rho)/(1:k))} where \code{rho} is a
#' parameter adjusting the NN dependence of infection probabilities.
#' 
#' Type II nonRL pattern assigns \eqn{n_1=}\code{round(n*ult.prop,0)} of them as cases,
#' and the rest as controls with first selecting \eqn{k_0=}\code{round(n*init.prop,0)} as cases initially, then selecting
#' a contagious case and then assigning the label case to the remaining points with infection probabilities 
#' inversely proportional to their position in the \code{k}NNs.
#' 
#' Type III nonRL pattern assigns \eqn{n_1=}\code{round(n*prop,0)} of them as cases,
#' and the rest as controls with first selecting a point, \eqn{Z_i}, as a case and assigning the
#' label case to the remaining points with infection probabilities \eqn{prob=rho (1-d_{ij}/d_{\max})^{pow}} where \eqn{d_{ij}} is the
#' distance from \eqn{Z_j} to \eqn{Z_i} for \eqn{j \ne i}, \eqn{d_{\max}} is the maximum of  \eqn{d_{ij}}  values, \code{rho} is a scaling parameter for
#' the infection probabilities and \code{pow} is a parameter in the power adjusting the distance dependence.
#' 
#' Type IV nonRL pattern assigns \eqn{n_1=}\code{round(n*ult.prop,0)} of them as cases,
#' and the rest as controls with first selecting \eqn{k_0=}\code{round(n*init.prop,0)} as cases initially and assigning the
#' label case to the remaining points with infection probabilities equal to the scaled bivariate normal density values
#' at those points.
#' 
#' The number of cases in Types I and III will be \eqn{n_1} on the average if the argument \code{poisson=TRUE}
#' (i.e., \eqn{n_1=}\code{rpois(1,round(n*prop,0))} ), otherwise \eqn{n_1=}\code{round(n*prop,0)}.
#' The initial and ultimate number of cases in Types II and IV will be \eqn{k_0} and \eqn{n_1} on the average if the argument
#' \code{poisson=TRUE} (i.e., \eqn{k_0=}\code{rpois(1,round(n*init.prop,0)}) and \eqn{n_1=}\code{rpois(1,round(n*ult.prop,0))}), otherwise
#' they will be exactly equal to \eqn{n_1=}\code{round(n*ult.prop,0)} and \eqn{k_0=}\code{round(n*init.prop,0)}.
#' 
#' At each type, we stop when we first exceed \eqn{n_1} cases. That is, the procedure ends when number of cases \eqn{n_c}
#' exceed \eqn{n_1}, and \eqn{n_c-n_1} of the cases (other than the initial case(s)) are randomly selected and relabeled as
#' controls, i.e. 0s, so that the number of cases is exactly \eqn{n_1}.
#' 
#' In the output cases are labeled as 1 and controls as 0, and initial contagious case is marked with a red cross
#' in the plot of the pattern.
#' 
#' See \insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat} and the functions \code{\link{rnonRLI}},
#' \code{\link{rnonRLII}}, \code{\link{rnonRLIII}}, and \code{\link{rnonRLIV}} for more detail on each type of
#' non-RL pattern.
#' 
#' Although the non-RL pattern is described for the case-control setting, it can be adapted for any two-class
#' setting when it is appropriate to treat one of the classes as cases or one of the classes behave like cases
#' and other class as controls.
#' 
#' The parameters of the non-RL patterns are specified in the argument \code{par.vec}, and the logical arguments \code{rand.init}
#' and poisson pass on to the types where required. \code{rand.init} is not used in type I but used in all other types,
#' poisson is used in all types, and init.from.cases is used in type I non-RL only.
#' 
#' @param dat A set of points the non-RL procedure is applied to obtain cases and controls randomly in the 
#' \code{type=type} fashion (see the description).
#' @param par.vec The parameter vector. It is \code{c(prop,k,rho)} for type I, \code{c(k,rho,pow,init.prop,ult.prop)}
#' for type II, \code{c(prop,rho,pow)} for type III, and \code{c(init.prop,ult.prop,s1,s2,rho)} for type IV non-RL patterns.
#' The parameters must be entered in this order in \code{par.vec} as a vector.
#' See the respective functions for more detail on the parameters.
#' @param type The type of the non-RL pattern. Takes on values \code{"I"}-\code{"IV"} for types I-IV non-RL
#' patterns (see the description above).
#' @param rand.init A logical argument (default is \code{TRUE}) to determine the choice of the initial case(s)
#' in the data set, \code{dat} for types II-IV non-RL pattern. If \code{rand.init=TRUE} then the initial case(s) is (are)
#' selected randomly from the data points, and if \code{rand.init=FALSE}, the first one is labeled as a case for type III
#' and the first \code{init.prop*n} entries in the data set, \code{dat}, are labeled as the cases types II and IV.
#' @param poisson A logical argument (default is \code{FALSE}) to determine whether the number of cases is random or fixed.
#' In types II and IV initial and ultimate number of cases, \eqn{k_0} and \eqn{n_1}, will be random if \code{poisson=TRUE} and fixed
#' otherwise. In types I and III the number of cases, \eqn{n_1}, will be random if poisson=TRUEURE and fixed otherwise.
#' See the description.
#' @param init.from.cases A logical argument (default is \code{TRUE}) to determine whether the initial cases at each
#' round will be take from cases or controls in type I non-RL pattern. 
#' The initial cases are taken from cases if \code{init.from.cases=TRUE}, and from controls otherwise.
#' See the function \code{\link{rnonRLI}}.
#'   
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{="cc"} for the case-control patterns for RL or non-RL of the given data points, \code{dat}}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{\code{par.vec}, the parameters required for each type of non-RL pattern. See the description
#' in the parameter list.}
#' \item{lab}{The labels of the points as 1 for cases and 0 for controls after the nonRL procedure is
#' applied to the data set, \code{dat}. Cases are denoted as red dots and controls as black circles in the plot.}
#' \item{init.cases}{The initial cases in the data set, \code{dat}. Marked with red crosses in the plot of the points.}
#' \item{cont.cases}{The contagious cases in the data set, \code{dat} in type II non-RL pattern.
#' Denoted as blue points in the plot of the points.}
#' \item{gen.points,ref.points}{Both are \code{NULL} for this function, as initial set of points, \code{dat}, are provided
#' for all of the non-RL procedures.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of cases and controls.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and the reference points}
#'
#' @seealso \code{\link{rnonRLI}}, \code{\link{rnonRLII}}, \code{\link{rnonRLIII}}, and \code{\link{rnonRLIV}}
#'
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' #data generation
#' n<-40;  #try also n<-20; n<-100;
#' dat<-cbind(runif(n,0,1),runif(n,0,1))
#'
#' #Type I non-RL pattern
#' #c(prop,k,rho) for type I
#' prop<-.5; knn<-3; rho<- .3
#' prv<-c(prop,knn,rho)
#'
#' Xdat<-rnonRL(dat,type="I",prv) #labeled data 
#' # or try Xdat<-rnonRL(dat,type="I",prv) for type I non-RL
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #Type II non-RL pattern
#' #c(k,rho,pow,init.prop,ult.prop) for type II
#' rho<-.8; pow<-2; knn<-5; ip<-.3; up<-.5
#' prv<-c(knn,rho,pow,ip,up)
#'
#' Xdat<-rnonRL(dat,type="II",prv) #labeled data 
#' # or try Xdat<-rnonRL(dat,type="I",prv) for type I non-RL
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #Type III non-RL pattern
#' #c(prop,rho,pow) for type III
#' prop<- .5; rho<-.8; pow<-2
#' prv<-c(prop,rho,pow)
#'
#' Xdat<-rnonRL(dat,type="III",prv) #labeled data 
#' # or try Xdat<-rnonRL(dat,type="I",prv) for type I non-RL
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #Type IV non-RL pattern
#' #c(init.prop,ult.prop,s1,s2,rho) for type IV
#' ult<-.5; int<- .1; s1<-s2<-.4; rho<- .1
#' prv<-c(int,ult,s1,s2,rho)
#'
#' Xdat<-rnonRL(dat,type="IV",prv) #labeled data 
#' # or try Xdat<-rnonRL(dat,type="I",prv) for type I non-RL
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat) 
#'   
#' @export 
rnonRL <- function(dat,par.vec,type,rand.init=TRUE,poisson=FALSE,init.from.cases=TRUE)
{
  res <- switch(type,
         I = { res <- rnonRLI(dat,par.vec[1],par.vec[2],par.vec[3],poisson=poisson,init.from.cases = init.from.cases) },
         II = { res <- rnonRLII(dat,par.vec[1],par.vec[2],par.vec[3],par.vec[4],par.vec[5],rand.init=rand.init,poisson=poisson) },
         III = { res <- rnonRLIII(dat,par.vec[1],par.vec[2],par.vec[3],rand.init=rand.init,poisson=poisson) },
         IV = { res <- rnonRLIV(dat,par.vec[1],par.vec[2],par.vec[3],par.vec[4],par.vec[5],rand.init=rand.init,poisson=poisson)  }
  )
  
  if (is.null(res)) stop("Enter numbers 1-4 or I-IV in quotes for type")
  
  res
} #end for the function
#'

#################################################################

#' @title Generation of Points under Segregation of Two Classes
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n_i} 2D points from class \eqn{j} with parameters  \eqn{r_j} for \eqn{j=1,2}.
#' The generated points are from two different classes which are segregated from each other.
#' The pattern generation starts with the initial points \code{X1.init} and \code{X2.init} (with default=\code{NULL} for both).
#' If both \code{X1.init=NULL} and \code{X2.init=NULL}, both \code{X1.init} and \code{X2.init} are generated uniformly in the unit square.
#' If only \code{X1.init=NULL}, \code{X1.init} is the sum of a point uniformly generated in the unit square and \code{X2.init} and
#' if only \code{X2.init=NULL}, \code{X2.init} is the sum of a point uniformly generated in the unit square and \code{X1.init}.
#' After the initial points from each class are available, \eqn{n_j} points from class \eqn{j} are generated
#' as \code{Xj[i,]<-Xj[(i-1),]+ru*c(cos(tu),sin(tu))}  
#' where \code{ru<-runif(1,0,rj)} and \code{tu<-runif(1,0,2*pi)} for \eqn{i=2,\ldots,n_j}
#' with \code{Xj[1,]=Xj.init} for \eqn{j=1,2}.
#' That is, at each step the new point in class \eqn{j} is generated within a circle with radius equal to \eqn{r_j}
#' (uniform in the polar coordinates).
#' Note that, the level of segregation is stronger if the initial points are further apart, and the level
#' of segregation increases as the radius values gets smaller.
#' 
#' @param n1,n2 Positive integers representing the number of class 1 and class 2 (i.e. \eqn{X_1} and \eqn{X_2}) points
#' to be generated under the segregation pattern.
#' @param r1,r2 Positive real numbers representing the radius of attraction within class, i.e. radius of the
#' circle center and generated points are from the same class.  
#' @param X1.init,X2.init 2D points representing the initial points for the segregated classes, default=\code{NULL}
#' for both. If both \code{X1.init=NULL} and \code{X2.init=NULL}, both \code{X1.init} and \code{X2.init} are generated uniformly in the
#' unit square.
#' If only \code{X1.init=NULL}, \code{X1.init} is the sum of a point uniformly generated in the unit square and \code{X2.init} and
#' if only \code{X2.init=NULL}, \code{X2.init} is the sum of a point uniformly generated in the unit square and \code{X1.init}.
#' The initial points are
#' marked with crosses in the plot of the points.
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{"2c"} for the 2-class pattern of segregation of the two classes}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{Radial (i.e. circular) within class radii of segregation, \code{r1} and \code{r2},
#' controlling the level of segregation}
#' \item{lab}{The class labels of the generated points, it is 1 class 1 or \eqn{X_1} points and 
#' 2 for class 2 or \eqn{X_2} points}
#' \item{init.cases}{The initial points for class 1 and class 2, one initial point for each class.}
#' \item{gen.points}{The output set of generated points (i.e. class 1 and class 2 points) segregated 
#' from each other.}
#' \item{ref.points}{The input set of reference points, it is \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers, which are the number of generated class 1 and class 2 points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated and
#' the initial points}
#'
#' @seealso \code{\link{rassoc}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20;  #try also n1<-10; n1<-100;
#' n2<-20; #try also n1<-40; n2<-50
#'
#' r1<-.3; r2<-.2
#'
#' #data generation
#' Xdat<-rseg(n1,n2,r1,r2) #labeled data
#' Xdat
#'
#' table(Xdat$lab)
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #with one initial point
#' X1init<-c(3,2)
#'
#' Xdat<-rseg(n1,n2,r1,r2,X1.init=X1init)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #with two initial points
#' X1init<-c(3,2)
#' X2init<-c(4,2)
#'
#' Xdat<-rseg(n1,n2,r1,r2,X1init,X2init)
#' Xdat
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat) 
#' 
#' @export 
rseg <- function(n1,n2,r1,r2,X1.init=NULL,X2.init=NULL)
{
  X1<-matrix(0,n1,2)
  X2<-matrix(0,n2,2)
  
  #initial class 1 and class 2 points
  cond<-which(c(is.null(X1.init)==TRUE & is.null(X2.init)==TRUE, #1
                is.null(X1.init)==TRUE & is.null(X2.init)==FALSE, #2
                is.null(X1.init)==FALSE & is.null(X2.init)==TRUE, #3
                is.null(X1.init)==FALSE & is.null(X2.init)==FALSE #4
  )==TRUE)
  
  switch(cond,
         "1" = { X1[1,]<-runif(2,0,1); X2[1,]<-runif(2,0,1) },
         "2" = { X1[1,]<-X2.init+runif(1,0,1); X2[1,]<-X2.init },
         "3" = { X1[1,]<-X1.init; X2[1,]<-X1.init+runif(1,0,1) },
         "4" = { X1[1,]<-X1.init; X2[1,]<-X2.init }
  )
  
  if (n1>1)
  {
    for (i in 2:n1)
    {
      ru<-runif(1,0,r1) 
      tu<-runif(1,0,2*pi)
      X1[i,]<-X1[(i-1),]+ru*c(cos(tu),sin(tu))
    }
  }
  if (n2>1)
  {
    for (j in 2:n2)
    {
      ru<-runif(1,0,r2) 
      tu<-runif(1,0,2*pi)
      X2[j,]<-X2[(j-1),]+ru*c(cos(tu),sin(tu))
    }
  }
  
  dat<-rbind(X1,X2)
  lab<-rep(c(1,2),c(n1,n2))
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  pname <-"radii of segregation"
  param<-c(r1,r2)
  names(param)<-c("r1","r2")
  typ<-paste("Segregation pattern with ",n1, " class 1 points and ", n2," class 2 points with ",pname," r1=",
             param[1]," and r2=",param[2],sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  rparam<-rep(0,2)
  ifelse(param[1]%%1==0,rparam[1]<-param[1],rparam[1]<-round(param[1],2))
  ifelse(param[2]%%1==0,rparam[2]<-param[2],rparam[2]<-round(param[2],2))
  
  txt<-"Segregation of Two Classes"
  main.txt<-bquote(atop("Segregation Pattern with Radii of Segregation",
                        ~r[1]== .(rparam[1]) ~"and" ~r[2]==.(rparam[2])))
  
  res<-list(
    pat.type="2c", #2c for 2-class pattern of segregation
    type=typ,
    parameters=param,
    lab=lab, #labels of the data points
    init.cases=rbind(X1[1,],X2[1,]), #initial points
    gen.points=dat, #generated points under segregation
    ref.points=NULL, 
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Uniform Points in a Circle
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \code{n} 2D points uniformly in the circle with center=\code{cent} and radius=\code{rad} using the rejection 
#' sampling approach (i.e., the function generates points in the smallest square containing the circle, keeping
#' only the points inside the circle until \eqn{n} points are generated).
#' The defaults for \code{cent=c(0,0)} and \code{rad=1}.
#' 
#' @param n A positive integer representing the number of points to be generated uniformly in the circle
#' @param cent A 2D point representing the center of the circle, with default=\code{c(0,0)}
#' @param rad A positive real number representing the radius of the circle.  
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{"1c"} for the 1-class pattern of the uniform data in the circle}
#' \item{type}{The type of the point pattern}
#' \item{parameters}{center of the circle, \code{cent}, and the radius of the circle, \code{rad}}
#' \item{lab}{The class labels of the generated points, \code{NULL} for this function, since points belong to the same 
#' class}
#' \item{init.cases}{The initial points, \code{NULL} for this function}
#' \item{gen.points}{The output set of generated points uniform in the circle.}
#' \item{ref.points}{The input set of reference points, it is \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The number of generated points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link[stats]{runif}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1);  #try also 10, 100, or 1000;
#' r<-.1; #try also r<-.3 or .5
#' cent<-c(1,2)
#'
#' #data generation
#' Xdat<-runif.circ(n,cent,r) #generated data
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
runif.circ <- function(n,cent=c(0,0),rad=1)
{  x1<-y1<-vector()
i<-1
while (i <= n)
{
  x<-runif(1,-rad,rad)
  y<-runif(1,-rad,rad)
  rad2<-rad^2
  if
 (x^2+y^2<=rad2)
  {
    x1<-cbind(x1,x)
    y1<-cbind(y1,y)
    i<-i+1
  }
}

pts0<-cbind(as.numeric(x1),as.numeric(y1))

dat<-t(t(pts0)+cent) #to translate the points to have center = cent

Xlim<-range(dat[,1])
Ylim<-range(dat[,2])

pname <-"center and radius"
param<-c(cent,rad)
names(param)<-c("c1","c2","radius")
typ<-paste(n, " Uniform Points in the Circle with Center (c1,c2) = (",cent[1],",",cent[2], ") and Radius = ", rad,sep="")

npts<-n

rparam<-rep(0,3)
for (i in 1:3)
  ifelse(param[i]%%1==0,rparam[i]<-param[i],rparam[i]<-round(param[i],2))

txt<-"Uniform Distribution in a Circle"
main.txt<-bquote(atop("Uniform Points in the Circle with",
                      "center" == ~(.(rparam[1]) *"," ~.(rparam[2])) ~"and "~r[0]== .(rparam[3])))
res<-list(
  pat.type="1c", #1c for 1-class pattern from a distribution
  type=typ,
  parameters=param,
  lab=NULL, #labels of the data points
  init.cases=NULL, #initial points
  gen.points=dat, #generated points uniformly in the circle
  ref.points=NULL, 
  desc.pat=txt, #description of the pattern
  mtitle=main.txt,
  num.points=npts,
  xlimit=Xlim,
  ylimit=Ylim
)

class(res)<-"Patterns"
res$call <-match.call()
res
} #end for the function
#'

#################################################################

#' @title Generation of Points with Clusters along the First Diagonal
#'
#' @description
#' An object of class \code{"Clusters"}.
#' 
#' Generates \code{n} 2D points with \code{k} (\eqn{k \ge 2}) clusters along the first diagonal 
#' where about \eqn{n/k} points belongs to each cluster.
#' 
#' If \code{distribution="uniform"}, the points are uniformly generated in their square
#' supports where one square is the unit square (i.e., with vertices \eqn{(0,0), (1,0), (1,1),(0,1)}), and 
#' the others are unit squares translated \eqn{j \sqrt{2} d}  units along the first diagonal for \eqn{j=1,2,\ldots,k-1}
#' (i.e. with vertices \eqn{(j d,j d), (1+j d,j d), (1+j d,1+j d),(j d,1+j d)}). 
#' 
#' If \code{distribution="bvnormal"}, the points are generated from the bivariate normal distribution with means equal to the
#' centers of the above squares (i.e. for each cluster with \code{mean=}\eqn{((1+j d)/2,(1+j d)/2)} for \eqn{j=0,1,\ldots,k-1}
#' and the covariance matrix \eqn{sd I_2}, where \eqn{I_2} is the \eqn{2 \times 2} identity matrix.
#' 
#' Notice that the clusters are more separated, i.e., generated data indicates more clear clusters as \eqn{d} increases
#' in either positive or negative direction with \eqn{d=0} indicating one cluster in the data. For a fixed \eqn{d}, when \code{distribution="bvnormal"},
#' the clustering gets stronger if the variance of each component, \eqn{sd^2}, gets smaller, and clustering gets weaker
#' as the variance of each component gets larger where default is \eqn{sd=1/6}.
#' 
#' @param n A positive integer representing the number of points to be generated from the two clusters
#' @param k A positive integer representing the number of clusters to be generated
#' @param d Shift in the first diagonal indicating the level of clustering in the data. Larger absolute values in
#' either direction (i.e. positive or negative) would yield stronger clustering.
#' @param sd The standard deviation of the components of the bivariate normal distribution with default \eqn{sd=1/6}, 
#' used only when \code{distribution="bvnormal"}.
#' @param distribution The argument determining the distribution of each cluster. Takes on values \code{"uniform"} and
#' \code{"bvnormal"} whose centers are \eqn{d} units apart along the first diagonal direction.
#' 
#' @return A \code{list} with the elements
#' \item{type}{The type of the clustering pattern}
#' \item{parameters}{The number of clusters, \code{k}, the diagonal shift d representing the level of clustering
#' (for both distribution types) and standard deviation, \code{sd}, for the bivariate normal distribution only}
#' \item{gen.points}{The output set of generated points from the clusters.}
#' \item{desc.pat}{Description of the clustering pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The number of generated points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link{rhor.clust}} and \code{\link{rrot.clust}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1);  #try also n<-50; n<-1000;
#' d<-.5 #try also -75,.75, 1
#' k<-3 #try also 5
#'
#' #data generation
#' Xdat<-rdiag.clust(n,k,d)
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #data generation (bvnormal)
#' n<-20  #or try sample(1:20,1);  #try also n<-50; n<-1000;
#' d<-.5 #try also -.75,.75, 1
#' k<-3 #try also 5
#' Xdat<-rdiag.clust(n,k,d,distr="bvnormal") #try also Xdat<-rdiag.clust(n,k,d,sd=.09,distr="bvnormal")
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rdiag.clust <- function(n,k,d,sd=1/6,distribution=c("uniform","bvnormal"))
{
  distribution <-match.arg(distribution)
  if (length(distribution) > 1 || is.na(distribution))
    stop("distribution must be one \"uniform\", \"bvnormal\"")
  
  ifelse(distribution=="uniform",
         X1<-matrix(runif(2*n),ncol=2),
         X1<-matrix(rnorm(2*n,1/2,sd),ncol=2))
  
  kl<-sample(1:k, n, replace=T) #sample to choose the clusters for the points
  
  kl.add<-cbind((kl-1)*d,(kl-1)*d)
  X1<-X1+kl.add #s12=k points are translated (k-1)*d units in x and y directions
  
  Xlim<-range(X1[,1])
  Ylim<-range(X1[,2])
  
  pname <-"parameters"
  
  ifelse(distribution=="uniform",
         {param<-c(k,d)
         names(param)<-c("number of clusters","diagonal shift")},
         {param<-c(k,d,sd)
         names(param)<-c("number of clusters","diagonal shift","std dev")})
  
  rparam<-rep(0,3)
  for (i in 2:3)
  {  ifelse(param[i]%%1==0,rparam[i]<-param[i],rparam[i]<-round(param[i],2))  }
  
  if (distribution=="uniform")
  {typ<-paste("Clustering pattern with uniform data on diagonally shifting squares with ",pname," k = ",param[1]," and d = ",param[2],sep="")
  txt<-"uniform clustering (with diagonally shifting squares)"
  main.txt<-bquote("Uniform Clusters with Parameters k" ==.(k) *"," ~delta ==.(rparam[2]) )
  } else
  {
    typ<-paste("Clustering pattern with bivariate normal data with diagonally shifting centers with ",pname," k = ",param[1]," d = ",param[2]," and sd = ",param[3],sep="")
    txt<-"bivariate normal clustering (with diagonally shifting centers)"
    main.txt<-bquote("Bivariate Normal Clusters with" ~ .(pname)~ "k" ==.(k) *"," ~delta ==.(rparam[2]) *", sd = " ~.(rparam[3]) )  
  }
  
  res<-list(
    type=typ,
    parameters=param,
    gen.points=X1, #generated points according to the clustering pattern
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=n,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<- "Clusters"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Points with Clusters along the Horizontal Axis
#'
#' @description
#' An object of class \code{"Clusters"}.
#' 
#' Generates \code{n} 2D points with \code{k} (\eqn{k \ge 2}) clusters along the horizontal axis
#' where about \eqn{n/k} points belongs to each cluster.
#' 
#' If \code{distribution="uniform"}, the points are uniformly generated in their square
#' supports where one square is the unit square (i.e., with vertices \eqn{(0,0), (1,0), (1,1),(0,1)}), and 
#' the others are \eqn{d} units shifted horizontally from each other so that their lower end vertices are
#' \eqn{(j-1)+(j-1) d} for \eqn{j=1,2,\ldots,k}.
#' 
#' If \code{distribution="bvnormal"}, the points are generated from the bivariate normal distribution with means equal to the
#' centers of the above squares (i.e. for each cluster with mean=(j+(j-1)d-1/2,1/2) for \eqn{j=1,2,\ldots,k}
#' and the covariance matrix \eqn{sd I_2}, where \eqn{I_2} is the \eqn{2 \times 2} identity matrix.
#' 
#' Notice that the clusters are more separated, i.e., generated data indicates more clear clusters as \eqn{d} increases
#' in either direction with \eqn{d=0} indicating one cluster in the data. For a fixed \eqn{d}, when \code{distribution="bvnormal"},
#' the clustering gets stronger if the variance of each component, \eqn{sd^2}, gets smaller, and clustering gets weaker
#' as the variance of each component gets larger where default is \eqn{sd=1/6}.
#' 
#' @param n A positive integer representing the number of points to be generated from all the clusters
#' @param k A positive integer representing the number of clusters to be generated
#' @param d Horizontal shift indicating the level of clustering in the data. Larger absolute values in either
#' direction (i.e. positive or negative) would yield stronger clustering.
#' @param sd The standard deviation of the components of the bivariate normal distribution with default \eqn{sd=1/6}, 
#' used only when \code{distribution="bvnormal"}.
#' @param distribution The argument determining the distribution of each cluster. Takes on values \code{"uniform"} and
#' \code{"bvnormal"} whose centers are \eqn{d} units apart along the horizontal direction.
#' 
#' @return A \code{list} with the elements
#' \item{type}{The type of the clustering pattern}
#' \item{parameters}{The number of clusters, \code{k}, and the horizontal shift, \code{d}, representing the level of clustering
#' (for both distribution types) and standard deviation, \code{sd}, for the bivariate normal distribution only.}
#' \item{gen.points}{The output set of generated points from the \code{k} clusters.}
#' \item{desc.pat}{Description of the clustering pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The number of generated points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link{rdiag.clust}} and \code{\link{rrot.clust}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-100;  #try also n<-50; or n<-1000;
#' d<-.5 #try also -.5,.75, 1
#' k<-3 #try also 5
#'
#' #data generation
#' Xdat<-rhor.clust(n,k,d)
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #data generation (bvnormal)
#' n<-100;  #try also n<-50; n<-1000;
#' d<-.1 #try also -.1, .75, 1
#' k<-3 #try also 5
#' Xdat<-rhor.clust(n,k,d,distr="bvnormal") #try also Xdat<-rhor.clust(n,k,d,sd=.15,distr="bvnormal")
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rhor.clust <- function(n,k,d,sd=1/6,distribution=c("uniform","bvnormal"))
{
  distribution <-match.arg(distribution)
  if (length(distribution) > 1 || is.na(distribution))
    stop("distribution must be one \"uniform\", \"bvnormal\"")
  
  if (d<0 & distribution=="uniform")
  {stop('d must be nonnegative for uniform square clusters')}
  
  ifelse(distribution=="uniform",X1<-matrix(runif(2*n),ncol=2),
         X1<-matrix(rnorm(2*n,1/2,sd),ncol=2))
  
  kl<-sample(1:k, n, replace=T)  #sample to choose the clusters for the points
  
  kl.add<-cbind((kl-1)+(kl-1)*d,0)
  X1<-X1+kl.add
  
  Xlim<-range(X1[,1])
  Ylim<-range(X1[,2])
  
  pname <-"parameters"
  
  ifelse(distribution=="uniform",
         {param<-c(k,d)
         names(param)<-c("number of clusters","horizontal shift")},
         {param<-c(k,d,sd)
         names(param)<-c("number of clusters","horizontal shift","std dev")})
  
  rparam<-rep(0,3)
  for (i in 2:3)
  {  ifelse(param[i]%%1==0,rparam[i]<-param[i],rparam[i]<-round(param[i],2))  }
  
  if (distribution=="uniform")
  {typ<-paste("Clustering pattern with uniform data on horizontally shifting squares with ",pname," k = ",param[1]," and d = ",param[2],sep="")
  txt<-"uniform clustering (with horizontally shifting squares)"
  main.txt<-bquote("Uniform Clusters with Parameters k" ==.(k) *"," ~delta ==.(rparam[2]) )
  } else
  {
    typ<-paste("Clustering pattern with bivariate normal data with horizontally shifting centers with ",pname," k = ",param[1]," d = ",param[2]," and sd = ",param[3],sep="")
    txt<-"bivariate normal clustering (with horizontally shifting centers)"
    main.txt<-bquote("Bivariate Normal Clusters with" ~ .(pname)~ "k" ==.(k) *"," ~delta ==.(rparam[2]) *", sd = " ~.(rparam[3]) )  
  }
  
  res<-list(
    type=typ,
    parameters=param,
    gen.points=X1, #generated points according to the clustering pattern
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=n,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<- "Clusters"
  res$call <-match.call()
  res
} #end of the function
#'
 
#################################################################

#' @title Generation of Points with Rotational Clusters
#'
#' @description
#' An object of class \code{"Clusters"}.
#' 
#' Generates \code{n} 2D points with \code{k} (\eqn{k \ge 2}) clusters with centers d unit away from origin and angles
#' between the rays joining successive centers and origin is \eqn{2 \pi/k} where about \eqn{n/k} points belongs to each cluster.
#' 
#' If \code{distribution="uniform"}, the points are uniformly generated in their square
#' supports with unit edge lengths and centers at \eqn{(d \cos(j 2 \pi/k),d \cos(j 2\pi/k))} for \eqn{j=1,2,\ldots,k}.
#' 
#' If \code{distribution="bvnormal"}, the points are generated from the bivariate normal distribution with means equal to the
#' centers of the above squares (i.e. for each cluster with \code{mean=}\eqn{(d \cos(j 2 \pi/k),d \cos(j 2\pi/k))}
#' for \eqn{j=1,2,\ldots,k} and the covariance matrix \eqn{sd I_2}, where \eqn{sd=d\sqrt{2 (1-cos(2 \pi/k))}/3}
#' and \eqn{I_2} is the \eqn{2 \times 2} identity matrix.
#' 
#' Notice that the clusters are more separated, i.e., generated data indicates more clear clusters as \eqn{d} increases
#' in either direction with \eqn{d=0} indicating one cluster in the data. For a fixed \eqn{d}, when \code{distribution="bvnormal"},
#' the clustering gets stronger if the variance of each component, \eqn{sd^2}, gets smaller, and clustering gets weaker
#' as the variance of each component gets larger where default is \eqn{sd=d\sqrt{2 (1-cos(2 \pi/k))}/3}.
#' 
#' @param n A positive integer representing the number of points to be generated from all the clusters
#' @param k A positive integer representing the number of clusters to be generated
#' @param d Radial shift indicating the level of clustering in the data. Larger absolute values in either
#' direction (i.e. positive or negative) would yield stronger clustering.
#' @param sd The standard deviation of the components of the bivariate normal distribution with default 
#' \eqn{sd=d\sqrt{2 (1-cos(2 \pi/k))}/3}, used only when \code{distribution="bvnormal"}.
#' @param distribution The argument determining the distribution of each cluster. Takes on values \code{"uniform"} and
#' \code{"bvnormal"} whose centers are \eqn{d} units apart along the horizontal direction.
#' 
#' @return A \code{list} with the elements
#' \item{type}{The type of the clustering pattern}
#' \item{parameters}{The number of clusters, \code{k}, and the radial shift, \code{d}, representing the level of clustering
#' (for both distribution types) and standard deviation, \code{sd}, for the bivariate normal distribution only.}
#' \item{gen.points}{The output set of generated points from the \code{k} clusters.}
#' \item{desc.pat}{Description of the clustering pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The number of generated points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link{rdiag.clust}} and \code{\link{rhor.clust}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-100;  #try also n<-50; n<-1000;
#' d<- 1.5 #try also -1, 1, 1.5, 2
#' k<-3 #try also 5
#' #data generation
#' Xdat<-rrot.clust(n,k,d)
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #data generation (bvnormal)
#' n<-100;  #try also n<-50; n<-1000;
#' d<- 1.5 #try also -1, 1, 1.5, 2
#' k<-3 #try also 5
#' Xdat<-rrot.clust(n,k,d,distr="bvnormal") #also try Xdat<-rrot.clust(n,k,d,sd=.5,distr="bvnormal")
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rrot.clust <- function(n,k,d,sd=d*sqrt(2*(1-cos(2*pi/k)))/3,distribution=c("uniform","bvnormal"))
{
  distribution <-match.arg(distribution)
  if (length(distribution) > 1 || is.na(distribution))
    stop("distribution must be one \"uniform\", \"bvnormal\"")
  
  ifelse(distribution=="uniform",X1<-matrix(runif(2*n,0,1)-1/2,ncol=2),
         X1<-matrix(rnorm(2*n,0,sd),ncol=2))
  
  kl<-sample(1:k, n, replace=T) #sample to choose the clusters for the points
  
  theta<-2*pi/k
  kl.add<-d*c(cos(kl*theta),sin(kl*theta))
  
  X1<-X1+kl.add #s12=k points are translated rotationally around (1/2,1/2)
  
  Xlim<-range(X1[,1])
  Ylim<-range(X1[,2])
  
  pname <-"parameters"
  
  ifelse(distribution=="uniform",
         {param<-c(k,d)
         names(param)<-c("number of clusters","radial shift")},
         {param<-c(k,d,sd)
         names(param)<-c("number of clusters","radial shift","std dev")})
  
  rparam<-rep(0,3)
  for (i in 2:3)
  {  ifelse(param[i]%%1==0,rparam[i]<-param[i],rparam[i]<-round(param[i],2))  }
  
  
  if (distribution=="uniform")
  {typ<-paste("Clustering pattern with uniform data on rotationally shifting squares with ",pname," k = ",param[1]," and d = ",param[2],sep="")
  txt<-"uniform clustering (with rotationally shifting squares)"
  main.txt<-bquote("Uniform Clusters with Parameters k" ==.(k) *"," ~delta ==.(rparam[2]) )
  } else
  {
    typ<-paste("Clustering pattern with bivariate normal data with rotationally shifting centers with ",pname," k = ",param[1]," d = ",param[2]," and sd = ",param[3],sep="")
    txt<-"bivariate normal clustering (with rotationally shifting centers)"
    main.txt<-bquote("Bivariate Normal Clusters with" ~ .(pname)~ "k" ==.(k) *"," ~delta ==.(rparam[2]) *", sd = " ~.(rparam[3]) )  
  }
  
  res<-list(
    type=typ,
    parameters=param,
    gen.points=X1, #generated points according to the clustering pattern
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=n,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<- "Clusters"
  res$call <-match.call()
  res
} #end of the function
#'

################################################

#' @title Generation of Points from Self Correspondence Pattern
#'
#' @description
#' An object of class \code{"Patterns"}.
#' 
#' Generates \eqn{n_1} 2D points from class 1 and  \eqn{n_2} (denoted as \code{n2} as an argument)
#' 2D points from class 2 in such a way that
#' self-reflexive pairs are more frequent than expected under CSR independence.
#'   
#' If \code{distribution="uniform"}, the points from class 1, say \eqn{X_i} are generated as follows: 
#' \eqn{X_i \stackrel{iid}{\sim} Uniform(S_1)} for \eqn{S_1=(c1r[1],c1r[2])^2} for \eqn{i=1,2,\ldots,n_{1h}}
#' where \eqn{n_{1h}=\lfloor n_1/2 \rfloor},
#' and for \eqn{k=n_{1h},+1,\ldots,n_1}, \eqn{X_k=X_{k-n_{1h}}+r (\cos(T_k), \sin(T_k))} where \eqn{r \sim Uniform(0,r_0)}
#' and \eqn{T_k} are iid \eqn{\sim Uniform(0,2 \pi)}.
#' Similarly, the points from class 2, say \eqn{Y_j} are generated as follows: 
#' \eqn{Y_j \stackrel{iid}{\sim} Uniform(S_2)} for \eqn{S_2=(c2r[1],c2r[2])^2} for \eqn{j=1,2,\ldots,n_{2h}} where \eqn{n_{2h}=\lfloor n_2/2\rfloor)},
#' and for \eqn{l=n_{2h},+1,\ldots,n_2}, \eqn{Y_l=Y_{l-n_{2h}}+r (\cos(T_l), \sin(T_l))} where \eqn{r \sim Uniform(0,r_0)} and
#' \eqn{T_l \stackrel{iid}{\sim} Uniform(0,2 \pi)}.
#' This version is the case IV in the article (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat}).
#' 
#' If \code{distribution="bvnormal"}, the points from class 1, say \eqn{X_i} are generated as follows: 
#' \eqn{X_i \stackrel{iid}{\sim} BVN(CM(S_1),I_{2x})} where \eqn{CM(S_1)} is the center of mass of \eqn{S_1} and I_{2x} is a \eqn{2 \times 2} matrix with diagonals
#' equal to \eqn{s_1^2} with \eqn{s_1=(c1r[2]-c1r[1])/3} and off-diagonals are 0 for \eqn{i=1,2,\ldots,n_{1h}} where \eqn{n_{1h}=\lfloor{n_1/2\rfloor}},
#' and for \eqn{k=n_{1h}+1,\ldots,n_1}, \eqn{X_k = Z_k+r (\cos(T_k), \sin(T_k))} where \eqn{Z_k \sim BVN(X_{k-n_{1h}}, I_2(r_0))}
#' with \eqn{I_2(r_0)} being the \eqn{2 \times 2} matrix with diagonals \eqn{r_0/3} and 0 off-diagonals, \eqn{r \sim Uniform(0,r_0)} and
#' \eqn{T_k} are iid \eqn{\sim Uniform(0,2 \pi)}.
#' Similarly, the points from class 2, say \eqn{Y_j} are generated as follows: 
#' \eqn{Y_j \stackrel{iid}{\sim} BVN(CM(S_2),I_{2y})} where \eqn{CM(S_1)} is the center of mass of \eqn{S_1} and I_{2y} is a \eqn{2 \times 2} matrix with diagonals
#' equal to \eqn{s_2^2} with \eqn{s_2=(c2r[2]-c2r[1])/3} and off-diagonals are 0 for \eqn{j=1,2,\ldots,n_{2h}} where \eqn{n_{2h}=\lfloor n_2/2\rfloor)},
#' and for \eqn{l=n_{2h},+1,\ldots,n_2}, \eqn{Y_l = W_k+r (\cos(T_l), \sin(T_l))} where \eqn{W_l \sim BVN(Y_{l-n_{2h}}, I_2(r_0))}
#' with \eqn{I_2(r_0)} being the \eqn{2 \times 2} matrix with diagonals \eqn{r_0/3} and 0 off-diagonals, \eqn{r \sim Uniform(0,r_0)} and
#' \eqn{T_l \stackrel{iid}{\sim} Uniform(0,2 \pi)}.
#' 
#' Notice that the classes will be segregated if the supports \eqn{S_1} and \eqn{S_2} are separated, with more separation
#' implying stronger segregation. Furthermore, \eqn{r_0} (denoted as \code{r0} as an argument) determines the level of self-reflexivity or self correspondence,
#' i.e. smaller \eqn{r_0} implies a higher level of self correspondence and vice versa for higher \eqn{r_0} .
#' 
#' See also (\insertCite{ceyhan:NNCorrespond2018;textual}{nnspat})
#' and the references therein.
#' 
#' @param n1,n2 Positive integers representing the numbers of points to be generated from the two classes
#' @param c1r,c2r Ranges of the squares which constitute the supports of the two classes
#' @param r0 The radius of attraction which determines the level of self-reflexivity (or self correspondence) in 
#' both the uniform and bvnormal distributions for the two classes
#' @param distribution The argument determining the distribution of each class. Takes on values \code{"uniform"} and
#' \code{"bvnormal"} (see the description for the details).
#' 
#' @return A \code{list} with the elements
#' \item{pat.type}{\code{"2c"} for the 2-class pattern of self-correspondence of the two classes}
#' \item{type}{The type of the spatial pattern}
#' \item{parameters}{The radius of attraction \eqn{r_0} which determines the level of self-correspondence.}
#' \item{lab}{The class labels of the generated points, it is 1 class 1 or X1 points and 
#' 2 for class 2 or \eqn{X_2} points}
#' \item{init.cases}{The initial points for class 1 and class 2, one initial point for each class, marked
#' with a cross in the plot.}
#' \item{gen.points}{The output set of generated points from the self-correspondence pattern.}
#' \item{ref.points}{The input set of reference points, it is \code{NULL} for this function.}
#' \item{desc.pat}{Description of the species correspondence pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{num.points}{The number of generated points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}- and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link{Zself.ref}} and \code{\link{Xsq.spec.cor}}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-50;  #try also n1<-50; n1<-1000;
#' n2<-50;  #try also n2<-50; n2<-1000;
#'
#' c1r<-c(0,1) #try also c(0,5/6), C(0,3/4), c(0,2/3)
#' c2r<-c(0,1) #try also c(1/6,1), c(1/4,1), c(1/3,1)
#' r0<-1/9 #try also 1/7, 1/8
#'
#' #data generation
#' Xdat<-rself.ref(n1,n2,c1r,c2r,r0)
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#'
#' #data generation (bvnormal)
#' Xdat<-rself.ref(n1,n2,c1r,c2r,r0,distr="bvnormal")
#' Xdat
#'
#' summary(Xdat)
#' plot(Xdat,asp=1)
#' plot(Xdat)
#' 
#' @export 
rself.ref <- function(n1,n2,c1r,c2r,r0,distribution=c("uniform","bvnormal")) 
{
  distribution <-match.arg(distribution)
  if (length(distribution) > 1 || is.na(distribution))
    stop("distribution must be one \"uniform\", \"bvnormal\"")
  
  if (n1<=1 || n2<=1)
  {stop('class sizes n1 and n2 must both be larger than 1')}
  
  X1<-matrix(0,n1,2)
  X2<-matrix(0,n2,2)
  
  x1<-c1r[1]; x2<-c1r[2]
  y1<-c2r[1]; y2<-c2r[2]
  
  n1f<-floor(n1/2); #floor of n1/2
  n2f<-floor(n2/2); #floor of n2/2
  
  if (distribution=="uniform")
  {
    X1[1:n1f,]<-matrix(runif(2*n1f,x1,x2),ncol=2)
    X2[1:n2f,]<-matrix(runif(2*n2f,y1,y2),ncol=2)
  } else
  { m1<-(x1+x2)/2; 
  sd1<-(x2-x1)/3
  
  m2<-(y1+y2)/2; 
  sd2<-(y2-y1)/3
  
  X1[1:n1f,]<-matrix(rnorm(2*n1f,m1,sd1),ncol=2)
  X2[1:n2f,]<-matrix(rnorm(2*n2f,m2,sd2),ncol=2)  
  }
  
  if (distribution=="uniform")
  {
    for (i in (n1f+1):n1)
    {
      ru<-runif(1,0,r0) 
      tu<-runif(1,0,2*pi)
      X1[i,]<-X1[i-n1f,]+ru*c(cos(tu),sin(tu))
    }
    for (j in (n2f+1):n2)
    {
      ru<-runif(1,0,r0) 
      tu<-runif(1,0,2*pi)
      X2[j,]<-X2[j-n2f,]+ru*c(cos(tu),sin(tu))
    }
  } else
  {
    for (i in (n1f+1):n1)
    {
      mi<-X1[i-n1f,]; sdi<-r0/3
      X1[i,]<-c(rnorm(1,mi[1],sdi),rnorm(1,mi[2],sdi))
    }
    for (j in (n2f+1):n2)
    {
      mj<-X2[j-n2f,]; sdj<-r0/3
      X2[j,]<-c(rnorm(1,mj[1],sdj),rnorm(1,mj[2],sdj))
    } 
  }
  
  dat<-rbind(X1,X2)
  lab<-rep(c(1,2),c(n1,n2))
  
  Xlim<-range(dat[,1])
  Ylim<-range(dat[,2])
  
  pname <-"radius of attraction"
  param<-c(r0)
  names(param)<-c("r0")
  typ<-paste("Self-Reflexivity Pattern with ",n1, " Class 1 Points and ", n2," Class 2 Points with Radius of Attraction r0 = ",
             param,sep="")
  
  npts<-c(n1,n2)
  names(npts)<-c("n1","n2")
  
  ifelse(param%%1==0,rparam<-param,rparam<-round(param,2))
  
  txt<-"self-reflexivity for two classes"
  main.txt<-bquote("Self-Reflexivity with Radius of Attraction," ~r[0]== .(rparam))
  
  res<-list(
    pat.type="2c", #2c for bivariate pattern from a distribution
    type=typ,
    parameters=param,
    lab=lab, #labels of the data points
    init.cases=rbind(X1[1,],X2[1,]), #initial points
    gen.points=dat, #generated points associated with Y points
    ref.points=NULL, #attraction points, i.e., points to which generated points are associated
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )
  
  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Tocher's randomized correction to the exact \eqn{p}-value
#'
#' @description 
#' Tocher's modification is used for the Fisher's exact test on the contingency tables making it less conservative,
#' by including the probability for the current table based on a randomized test
#' (\insertCite{tocher:1950;textual}{nnspat}). It is applied When table-inclusive version of the \eqn{p}-value,
#' \eqn{p^>_{inc}}, is larger, but table-exclusive version, \eqn{p^>_{exc}}, is less than the level of the test \eqn{\alpha},
#' a random number, \eqn{U}, is generated from uniform distribution in \eqn{(0,1)}, and if \eqn{U \leq (\alpha-p^>_{exc})/p_t},
#' \eqn{p^>_{exc}} is used, otherwise \eqn{p_{inc}} is used as the \eqn{p}-value.
#' 
#' Table-inclusive and exclusive \eqn{p}-values are defined as follows. 
#' Let the probability of the contingency table itself
#' be \eqn{p_t=f(n_{11}|n_1,n_2,c_1;\theta)} where \eqn{\theta} is the odds ratio
#' under the null hypothesis (e.g. \eqn{\theta=1} under independence) and 
#' \eqn{f} is the probability mass function of the hypergeometric distribution.
#' In testing the one-sided alternative \eqn{H_o:\,\theta=1} versus \eqn{H_a:\,\theta>1},
#' let \eqn{p=\sum_S f(t|n_1,n_2,c_1;\theta=1)}, then
#' with \eqn{S=\{t:\,t \geq n_{11}\}}, we get the \emph{table-inclusive version} which is denoted as \eqn{p^>_{inc}}
#' and with \eqn{S=\{t:\,t> n_{11}\}}, we get the \emph{table-exclusive version}, denoted as \eqn{p^>_{exc}}.
#' 
#' See (\insertCite{ceyhan:exact-NNCT;textual}{nnspat}) for more details.
#' 
#' @param ptable Probability of the contingency table under the null hypothesis using the hypergeometric 
#' distribution for Fisher's exact test.
#' @param pval Table inclusive \eqn{p}-value for Fisher's exact test on the contingency table.
#'
#' @return A modified \eqn{p}-value based on the Tocher's randomized correction.
#'
#' @seealso \code{\link{prob.nnct}}, \code{\link{exact.pval1s}}, and \code{\link{exact.pval2s}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' ptab<-.03
#' pval<-.06
#' tocher.cor(ptab,pval)
#'
#' @export 
tocher.cor <- function(ptable,pval)
{
  pexcl<-pval-ptable
  PV<-pval
  if (pexcl<.05 || pval >= .05)
  {
    ifelse(runif(1)<= (.05-pexcl)/ptable,PV<-pexcl,PV<-pval)
  }
  PV
} #end for the function
#'

#################################################################

#' @title Probability of the current nearest neighbor contingency table
#'
#' @description 
#' Computes the probability of the observed \eqn{2 \times 2} nearest neighbor contingency table (NNCT) 
#' \eqn{p_t=f(n_{11}|n_1,n_2,c_1;\theta)} where \eqn{\theta=(n_1-1)(n_2-1)/(n_1 n_2)} which is the odds ratio
#' under RL or CSR independence and
#' \eqn{f} is the probability mass function of the hypergeometric distribution.
#' That is, given the margins of the current NNCT, the probability of obtaining the current table with the odds
#' ratio \eqn{\theta} being the value under the null hypothesis.
#' This value is used to compute the table-inclusive and exclusive \eqn{p}-values for the exact inference on NNCTs.
#' 
#' See (\insertCite{ceyhan:exact-NNCT;textual}{nnspat}) for more details.
#' 
#' @param ct A NNCT
#'
#' @return The probability of getting the observed NNCT, \code{ct} , under the null hypothesis.
#'
#' @seealso \code{\link{exact.pval1s}} and \code{\link{exact.pval2s}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' ct<-matrix(sample(20:40,4),ncol=2)
#' prob.nnct(ct)
#'
#' ct<-matrix(sample(20:40,4),ncol=2)
#' prob.nnct(ct)
#'
#' @export 
prob.nnct <- function(ct)
{
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr>2 || kc>2)
  {stop('number of classes must be 2 for this function')}
  
  n11<-ct[1,1]
  rs<-row.sum(ct)
  cs<-col.sum(ct)
  n<-sum(rs)
  
  n1<-rs[1]; n2<-rs[2]
  c1<-cs[1]
  
  or<-(n1-1)*(n2-1)/(n1*n2) #odds ratio under CSR or RL
  
  lo <- max(0,n1+c1-n)
  hi <- min(n1,c1)
  
  support <- lo:hi
  
  logdc <- dhyper(support, n1, n2, c1, log = TRUE)
  
  dnhyper <- function(ncp) { #density of non-central hyp geo distribution with ncp(noncent para)=ncp
    d <- logdc + log(ncp) * support
    d <- exp(d - max(d))
    d/sum(d)
  }
  
  ptab<-dnhyper(or)[support==n11]
  ptab
} #end for the function
#'

#################################################################

#' @title \eqn{p}-value correction to the one-sided version of exact NNCT test
#'
#' @description 
#' In using Fisher's exact test on the \eqn{2 \times 2} nearest neighbor contingency tables (NNCTs) a correction
#' may be needed for the \eqn{p}-value. For the one-sided alternatives, the probabilities of 
#' more extreme tables are summed up, including or excluding the 
#' probability of the table itself (or some middle way). 
#' Let the probability of the contingency table itself be \eqn{p_t=f(n_{11}|n_1,n_2,c_1;\theta_0)}
#' where \eqn{\theta_0=(n_1-1)(n_2-1)/(n_1 n_2)} which is the odds ratio
#' under RL or CSR independence and
#' \eqn{f} is the probability mass function of the hypergeometric distribution.
#' For testing the one-sided alternative \eqn{H_o:\,\theta=\theta_0} versus \eqn{H_a:\,\theta>\theta_0},
#' we consider the following four methods in calculating the \eqn{p}-value:
#' \itemize{
#' \item [(i)] with \eqn{S=\{t:\,t \geq n_{11}\}}, we get the
#' \emph{table-inclusive version} which is denoted as \eqn{p^>_{inc}},
#' \item [(ii)] with \eqn{S=\{t:\,t> n_{11}\}}, we get the
#' \emph{table-exclusive version}, denoted as \eqn{p^>_{exc}}.
#' \item [(iii)] Using \eqn{p=p^>_{exc}+p_t/2}, we get the \emph{mid-\eqn{p} version}, denoted as \eqn{p^>_{mid}}.
#' \item [(iv)] We can also use \emph{Tocher corrected version} which is denoted as \eqn{p^>_{Toc}}
#' (see \code{\link{tocher.cor}} for details).
#' }
#' 
#' See (\insertCite{ceyhan:exact-NNCT;textual}{nnspat}) for more details.
#' 
#' @param ptable Probability of the observed \eqn{2 \times 2} NNCT under the null hypothesis using the hypergeometric distribution
#' for Fisher's exact test.
#' @param pval Table inclusive \eqn{p}-value for Fisher's exact test on the NNCT.
#' @param type The type of the \eqn{p}-value correction for the one-sided exact test on the NNCT, default=\code{"inc"}.
#' Takes on values \code{"inc"}, \code{"exc"}, \code{"mid"}, \code{"tocher"} (or equivalently \code{1-4}, respectively) for table inclusive,
#' table-exclusive, mid-\eqn{p}-value, and Tocher corrected \eqn{p}-value, respectively.
#'
#' @return A modified \eqn{p}-value based on the correction specified in \code{type}.
#'
#' @seealso \code{\link{exact.pval2s}} and \code{\link{tocher.cor}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' ct<-matrix(sample(20:40,4),ncol=2)
#' ptab<-prob.nnct(ct)
#' pv<-.3
#' exact.pval1s(ptab,pv)
#' exact.pval1s(ptab,pv,type="exc")
#' exact.pval1s(ptab,pv,type="mid")
#'
#' @export 
exact.pval1s <- function(ptable,pval,type="inc") #exact pvalue for one-sided tests
{
  Pv <- switch(type,
         inc = { pv <- pval ;
         names(pv)<-"table-inclusive p-value (with more extreme tables in one direction than the current table)"},
         exc = { pv <- pval-ptable;
         names(pv)<-"table-exclusive p-value (with more extreme tables in one direction than the current table)"},
         mid = { pv <- pval-ptable/2;
         names(pv)<-"mid-p-value (with more extreme tables in one direction than the current table)"},
         tocher = { pv <- tocher.cor(ptable,pval);
         names(pv)<-"Tocher randomized p-value (with more extreme tables in one direction than the current table)"}
  )
  
  if (is.null(Pv)) stop("Enter numbers 1-4 or inc, exc, mid, tocher in quotes for type")
  
  pv
} #end for the function
#'

################################################

#' @title \eqn{p}-value correction to the two-sided version of exact NNCT test
#'
#' @description 
#' In using Fisher's exact test on the \eqn{2 \times 2} nearest neighbor contingency tables (NNCTs) a correction may be needed
#' for the \eqn{p}-value. For the one-sided alternatives, the probabilities of 
#' more extreme tables are summed up, including or excluding the 
#' probability of the table itself (or some middle way). 
#' 
#' There is additional complexity in \eqn{p}-values for the two-sided alternatives.
#' A recommended method is adding up probabilities of the same
#' size and smaller than the probability associated with the current table.
#' Alternatively, one can double the one-sided \eqn{p}-value (see (\insertCite{agresti:1992;textual}{nnspat}).
#' 
#' Let the probability of the contingency table itself be \eqn{p_t=f(n_{11}|n_1,n_2,c_1;\theta_0)}
#' where \eqn{\theta_0=(n_1-1)(n_2-1)/(n_1 n_2)} which is the odds ratio
#' under RL or CSR independence and
#' \eqn{f} is the probability mass function of the hypergeometric distribution.
#' 
#' **Type (I):** For double the one-sided \eqn{p}-value, we propose the following four variants:
#'   \itemize{
#' \item [(i)] twice the minimum of \eqn{p_{inc}} for the one-sided tests, which is
#' table-inclusive version for this type of two-sided test, and denoted as \eqn{p^I_{inc}},
#' \item [(ii)] twice the minimum of \eqn{p_{inc}} minus twice the table
#' probability \eqn{p_t}, which is table-exclusive version of this type of
#' two-sided test, and denoted as \eqn{p^I_{exc}},
#' \item [(iii)] table-exclusive version of this type of
#' two-sided test plus \eqn{p_t}, which is mid-\eqn{p}-value for
#' this test, and denoted as \eqn{p^I_{midd}},
#' \item [(iv)]Tocher corrected version (see \code{\link{tocher.cor}} for details).
#' }
#' 
#' **Type (II):** For summing the \eqn{p}-values of more extreme ---than that of the table--- cases
#' in both directions, the following variants are obtained.
#' The \eqn{p}-value is \eqn{p=\sum_S f(t|n_1,n_2,c_1;\theta=1)} with
#' \itemize{
#' \item [(i)] \eqn{S=\{t:\,f(t|n_1,n_2,c_1;\theta=1) \leq p_t\}}, which is
#' called \emph{table-inclusive version}, \eqn{p^{II}_{inc}},
#' \item [(ii)] the probability of the observed table is included twice, once for each side;
#' that is \eqn{p=p^{II}_{inc}+p_t}, which is called \emph{twice-table-inclusive version}, \eqn{p^{II}_{tinc}},
#' \item [(iii)] table-inclusive minus \eqn{p_t}, which is referred as \emph{table-exclusive version}, \eqn{p^{II}_{exc}},
#' \item [(iv)] table-exclusive plus one-half
#' the \eqn{p_t}, which is called \emph{mid-\eqn{p} version}, \eqn{p^{II}_{mid}} and,
#' \item [(v)]\emph{Tocher corrected version}, \eqn{p^{II}_{Toc}}, is obtained as before.
#' }
#' 
#' See (\insertCite{ceyhan:exact-NNCT;textual}{nnspat}) for more details.
#' 
#' @param ptable Probability of the observed \eqn{2 \times 2} NNCT under the null hypothesis using the hypergeometric
#' distribution for Fisher's exact test.
#' @param pval Table inclusive \eqn{p}-value for Fisher's exact test on the NNCT.
#' @param type The type of the \eqn{p}-value correction for the two-sided exact test on the NNCT, default=\code{"inc"}.
#' Takes on values \code{"inc"}, \code{"exc"}, \code{"mid"}, \code{"tocher"} (or equivalently \code{1-4}, respectively) for table inclusive,
#' table-exclusive, mid-\eqn{p}-value, and Tocher corrected \eqn{p}-value, respectively.
#' @param double A logical argument (default is \code{FALSE}) to determine whether type I or II correction should be 
#' applied to the two-sided \eqn{p}-value. If \code{TRUE} type I correction (for doubling the minimum of the one-sided \eqn{p}-value) 
#' is applied, otherwise, type II correction (using the probabilities for the more extreme tables) is applied.
#'
#' @return A modified \eqn{p}-value based on the correction specified in \code{type}.
#'
#' @seealso \code{\link{exact.pval1s}} and \code{\link{tocher.cor}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' ct<-matrix(sample(20:40,4),ncol=2)
#' ptab<-prob.nnct(ct)
#' pv<-.23
#' exact.pval2s(ptab,pv)
#' exact.pval2s(ptab,pv,type="exc")
#' exact.pval2s(ptab,pv,type="mid")
#'
#' @export 
exact.pval2s <- function(ptable,pval,type="inc",double=FALSE) #exact pvalue for two-sided tests
{
  if (double==FALSE)
  {
    Pv <- switch(type,
           inc = { pv <- pval;
           names(pv)<-"table-inclusive p-value (with more extreme tables in both directions than the current table)"},
           twice.inc = { pv <- pval+ptable;
           names(pv)<-"twcie-table-inclusive p-value (with more extreme tables in both directions than the current table)"},
           exc = { pv <- pval-ptable;
           names(pv)<-"table-exclusive p-value (with more extreme tables in both directions than the current table)"},
           mid = { pv <- pval-ptable/2;
           names(pv)<-"mid-p-value (with more extreme tables in both directions than the current table)"},
           tocher = { pv <- tocher.cor(ptable,pval);
           names(pv)<-"Tocher randomized p-value (with more extreme tables in both directions than the current table)"}
    )
  } else 
  {
    Pv <- switch(type,
           inc = { pv <- 2*min(pval,1-pval);
           names(pv)<-"table-inclusive p-value (with doubling the one-sided p-value)"},
           exc = { pv <- 2*min(pval,1-pval)-2*ptable;
           names(pv)<-"table-exclusive p-value (with doubling the one-sided p-value)"},
           mid = { pv <- 2*min(pval,1-pval)-ptable;
           names(pv)<-"mid-p-value (with doubling the one-sided p-value)"},
           tocher = { pv <- tocher.cor(ptable,pval);
           names(pv)<-"Tocher randomized p-value (with doubling the one-sided p-value)"}
    )
  }
  if (is.null(Pv)) stop("Enter numbers 1-4 or inc, twice.inc, exc, mid, tocher in quotes for type")
  
  pv
} #end for the function
#'

#################################################################

#' @title Exact version of Pearson's chi-square test on NNCTs 
#' 
#' @description
#' An object of class \code{"htest"} performing exact version of Pearson's chi-square test on nearest neighbor contingency
#' tables (NNCTs) for the RL or CSR independence for 2 classes.
#' Pearson's \eqn{\chi^2} test is based on the test statistic 
#' \eqn{\mathcal X^2=\sum_{j=1}^2\sum_{i=1}^2 (N_{ij}-\mu_{ij})^2/\mu_{ij}},
#' which has \eqn{\chi^2_1} distribution in the limit provided
#' that the contingency table is constructed under the independence null hypothesis.
#' The exact version of Pearson's test uses the exact distribution of \eqn{\mathcal X^2} rather than large sample 
#' \eqn{\chi^2} approximation.
#' That is, for the one-sided alternative, we calculate
#' the \eqn{p}-values as in the function \code{\link{exact.pval1s}};
#' and for the two-sided alternative, we calculate
#' the \eqn{p}-values as in the function \code{\link{exact.pval2s}} with double argument determining
#' the type of the correction. 
#' 
#' This test would be equivalent to Fisher's exact test \code{\link[stats]{fisher.test}} if the odds ratio=1
#' (which can not be specified in the current version), and the odds ratio for the RL or CSR independence null
#' hypothesis is \eqn{\theta_0=(n_1-1)(n_2-1)/(n_1 n_2)} which is used in the function and
#' the \eqn{p}-value and confidence interval computations are are adapted from \code{\link[stats]{fisher.test}}.
#' 
#' See \insertCite{ceyhan:SWJ-spat-sym2014;textual}{nnspat} for more details.
#' 
#' @param ct A \eqn{2 \times 2} NNCT
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for the odds ratio
#' @param pval.type The type of the \eqn{p}-value correction for the exact test on the NNCT, default=\code{"inc"}.
#' Takes on values \code{"inc"}, \code{"exc"}, \code{"mid"}, \code{"tocher"} (or equivalently \code{1-4}, respectively) for table inclusive,
#' table-exclusive, mid-\eqn{p}-value, and Tocher corrected \eqn{p}-value, respectively.
#' @param double A logical argument (default is \code{FALSE}) to determine whether type I or II correction should be 
#' applied to the two-sided \eqn{p}-value. Used only when \code{alternative="two.sided"}. 
#' If \code{TRUE} type I correction (for doubling the minimum of the one-sided \eqn{p}-value) 
#' is applied, otherwise, type II correction (using the probabilities for the more extreme tables) is applied.
#'   
#' @return A \code{list} with the elements
#' \item{statistic}{The test statistic, it is \code{NULL} for this function}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the odds ratio in the \eqn{2 \times 2} NNCT
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate, i.e., the observed odds ratio the \eqn{2 \times 2} NNCT.}
#' \item{null.value}{Hypothesized null value for the odds ratio in the \eqn{2 \times 2} NNCT, which is
#' \eqn{\theta_0=(n_1-1)(n_2-1)/(n_1 n_2)} for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the contingency table, \code{ct}}
#'  
#' @seealso \code{\link[stats]{fisher.test}}, \code{\link{exact.pval1s}}, and \code{\link{exact.pval2s}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' cls<-sample(1:2,n,replace = TRUE)  #or try cls<-rep(1:2,c(10,10))
#' ct<-nnct(ipd,cls)
#' ct
#'
#' exact.nnct(ct)
#' fisher.test(ct)
#'
#' exact.nnct(ct,alt="g")
#' fisher.test(ct,alt="g")
#'
#' exact.nnct(ct,alt="l",pval.type = "mid")
#'
#' #############
#' ct<-matrix(sample(10:20,9),ncol=3)
#' fisher.test(ct) #here exact.nnct(ct) gives error message, since number of classes > 2
#' 
#' @export
exact.nnct <- function(ct,alternative="two.sided",conf.level = 0.95,pval.type="inc",double=FALSE)
{
  kr<-nrow(ct); kc<-ncol(ct);
  if (kr!=2 || kc!=2)
  {stop('number of classes must be 2 for this function')}
  
  alternative <- char.expand(alternative, c("two.sided","less", "greater"))
  
  if (alternative!="two.sided" & double==TRUE)
  {stop('argument mismatch: double=TRUE can only be used with two.sided alternative only')}
  
  n11<-ct[1,1]
  rs<-row.sum(ct)
  cs<-col.sum(ct)
  n<-sum(rs)
  
  n1<-rs[1]; n2<-rs[2]
  c1<-cs[1]
  
  or<-(n1-1)*(n2-1)/(n1*n2) #odds ratio under CSR or RL
  
  lo <- max(0,n1+c1-n)
  hi <- min(n1,c1)
  
  support <- lo:hi
  
  logdc <- dhyper(support, n1, n2, c1, log = TRUE)
  
  dnhyper <- function(ncp) { #density of non-central hyp geo distribution with ncp(noncent para)=ncp
    d <- logdc + log(ncp) * support
    d <- exp(d - max(d))
    d/sum(d)
  }
  ptab<-dnhyper(or)[support==n11]
  
  if (is.na(ptab))
  {stop('The probability of the current table is NaN, so the exact test is not defined')}
  
  pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
    if (ncp == 1) {
      return(if (upper.tail) phyper(n11 - 1, n1, n2, c1,lower.tail = FALSE) else phyper(n11,n1, n2, c1))
    }
    if (ncp == 0) {
      return(as.numeric(if (upper.tail) q <= lo else q >= lo))
    }
    if (ncp == Inf) {
      return(as.numeric(if (upper.tail) q <= hi else q >= hi))
    }
    sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= q])
  }
  
  switch(alternative, 
         less = {pv<-pnhyper(n11, or);
         pval<-exact.pval1s(ptab,pv,pval.type)}, 
         greater ={ pv<-pnhyper(n11, or, upper.tail = TRUE);
         pval<-exact.pval1s(ptab,pv,pval.type)}, 
         two.sided = {
           if (or == 0) as.numeric(n11 == lo) 
           else if (or == Inf) as.numeric(n11 == hi) 
           else { relErr <- 1 + 10^(-7)
           d <- dnhyper(or)
           pv<-sum(d[d <= d[n11 - lo + 1] * relErr])
           pval<-exact.pval2s(ptab,pv,pval.type,double)
           }
         })
  
  ncp.U <- function(x, alpha) {
    if (x == hi) 
      return(Inf)
    p <- pnhyper(x, or)
    if (p < alpha) 
      uniroot(function(t) pnhyper(x, t) - alpha, c(0, 1))$root
    else if (p > alpha) 
      1/uniroot(function(t) pnhyper(x, 1/t) - alpha, c(.Machine$double.eps, 1))$root
    else 1
  }
  ncp.L <- function(x, alpha) {
    if (x == lo) 
      return(0)
    p <- pnhyper(x, or, upper.tail = TRUE)
    if (p > alpha) 
      uniroot(function(t) pnhyper(x, t, upper.tail = TRUE) - 
                alpha, c(0, 1))$root
    else if (p < alpha) 
      1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = TRUE) - 
                  alpha, c(.Machine$double.eps, 1))$root
    else 1
  }
  
  cint <- switch(alternative, 
                 less = c(0, ncp.U(n11,1 - conf.level)), 
                 greater = c(ncp.L(n11, 1 - conf.level),Inf),
                 two.sided = {alpha <- (1 - conf.level)/2
                 c(ncp.L(n11, alpha), ncp.U(n11, alpha))}
  )
  attr(cint, "conf.level") <-conf.level 
  
  method <-c("Fisher's Exact Test for NNCT with",names(pval))
  estimate<-(ct[1,1]*ct[2,2])/(ct[1,2]*ct[2,1])
  names(estimate) <-c("odds ratio")
  
  null.val<- or #Expected odds ratio under CSR or RL
  names(null.val) <-"odds ratio in NNCT"
  
  dname <-deparse(substitute(ct))
  
  rval <-list(
    statistic=NULL,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'

#################################################################

#' @title Interpoint Distance Matrix for Standardized Data
#'
#' @description 
#' This function computes and returns the distance matrix computed by using the specified distance measure 
#' to compute the distances between the rows of a data matrix which is standardized row or column-wise.
#' That is, the output is the interpoint distance (IPD) matrix of the rows of the given set of points \code{x}
#' \code{\link[stats]{dist}} function in the \code{stats} package of the standard R distribution.
#' The argument column is the logical argument (default=\code{TRUE}) to determine row-wise or column-wise standardization.
#' If \code{TRUE} each column is divided by its standard deviation, else each row is divided by its standard deviation.
#' This function is different from the \code{\link[stats]{dist}} function in the \code{stats} package.
#' \code{dist} returns the distance matrix in a lower triangular form, and dist.std.data returns in a full matrix
#' of distances of standardized data set.
#' \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @param x A set of points in matrix or data frame form where points correspond to the rows.
#' @param column A logical argument (default is \code{TRUE}) to determine whether standardization is row-wise or
#' column-wise. If \code{TRUE} it is column-wise else row-wise standardization.
#' @param \dots Additional parameters to be passed on the \code{dist} function.
#'
#' @return A distance matrix whose i,j-th entry is the distance between rows \eqn{i} and \eqn{j} of \code{x}, which is
#' standardized row-wise or column-wise.
#'
#' @seealso \code{\link[stats]{dist}}, \code{\link{ipd.mat}}, and \code{\link{ipd.mat.euc}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' range(ipd)
#'
#' ipd2<-dist.std.data(Y) #distance of standardized data
#' range(ipd2)
#'
#' ipd2<-dist.std.data(Y,method="max") #distance of standardized data
#' range(ipd2)
#'
#' #############
#' Y<-matrix(runif(60,0,100),ncol=3)
#' ipd<-ipd.mat(Y)
#' range(ipd)
#'
#' ipd2<-dist.std.data(Y) #distance of standardized data
#' range(ipd2)
#'
#' @export
dist.std.data <- function(x,column=TRUE,...)
{
  if (column==TRUE)
  {std<-apply(x,2,FUN=sd)
  x.std<-sweep(x,2,std,FUN="/")
  } else
  {
    std<-apply(x,1,FUN=sd)
    x.std<-sweep(x,1,std,FUN="/")
  }
  if (all(is.na(x.std)))
  {stop('The standardized data is not defined')}
  
  dis<-dist(x.std,...)
  dist2full(dis) #uses dist2full function to complete it to full matrix distance
} #end for the function
#' 

#################################################################

#' @title Smallest and Largest Distances in a Distance Matrix
#'
#' @description 
#' This function finds and returns the \code{k} smallest and \code{k} largest distances in a distance matrix or distance object,
#' and also provides pairs of objects these distances correspond to.
#' The code is adapted from 
#' \url{http://people.stat.sc.edu/Hitchcock/chapter1_R_examples.txt}.
#' 
#' @param ds A distance matrix or a distance object
#' @param k A positive integer representing the number of (min and max) distances to be presented, default is \eqn{k=1}
#'
#' @return A \code{list} with the elements
#' \item{min.dis}{The \code{k} smallest distances in \code{ds}}
#' \item{ind.min.dis}{The indices (i.e. row numbers) of the \code{k} pairs of object which has the
#' \code{k} smallest distances in \code{ds}}
#' \item{max.dis}{The \code{k} largest distances in \code{ds}}
#' \item{ind.max.dis}{The indices (i.e. row numbers) of the \code{k} pairs of object which has the
#' \code{k} largest distances in \code{ds}}
#'
#' @seealso \code{\link[stats]{dist}}, \code{\link{ipd.mat}}, and \code{\link{ipd.mat.euc}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' ipd<-ipd.mat(Y)
#' pick.min.max(ipd)
#' #or
#' pick.min.max(dist(Y))
#'
#' pick.min.max(ipd,2)
#'
#' @export
pick.min.max <- function(ds, k=1)
{
  ds <- as.matrix(ds)
  n<-nrow(ds)
  if (k>=n)
  {stop('k must be less than the sample size, n')}
  
  if (n<=1)
  {stop('sample size, n, must be >= 2')}
  
  if (is.matrix(ds)) {
    ds2 <- ds
    disLT <- ds[lower.tri(ds)]
  } else {
    ds2 <- dist2full(ds)
    disLT <- ds2[lower.tri(ds2)]
  }
  
  kmin <- matrix(sort(disLT)[1:k],ncol=k)
  kmax <- matrix(sort(disLT,decreasing=T)[1:k],ncol=k)
  
  kmin.pairs <- kmax.pairs <- matrix(0,nrow=k,ncol=2)
  for (i in 1:k){
    kmin.pairs[i,] <- which(ds2==kmin[i],arr.ind=T)[2,]
  }
  
  for (i in 1:k){
    kmax.pairs[i,] <- which(ds2==kmax[i],arr.ind=T)[2,]
  }
  
  if (k==1)
  {rownames(kmin)<-"smallest distance="
  rownames(kmax)<-"largest distance="
  colnames(kmin.pairs) <- colnames(kmax.pairs) <-c("corresponding to pair of objects","")
  } else
  {rownames(kmin)<-paste("smallest", k, "distances=")
  rownames(kmax)<-paste("largest", k, "distances=")
  colnames(kmin.pairs) <- colnames(kmax.pairs) <-c("corresponding to pairs of objects","")
  }
  
  list(
    min.dis=kmin,
    ind.min.dis=kmin.pairs,
    max.dis=kmax,
    ind.max.dis=kmax.pairs
  )
} #end for the function
#' 

################################################
###FUNCTIONS for CUZICK and EDWARDS k-NN TESTS###
################################################

#' @title Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' This function computes Cuzick and Edwards \eqn{T_k} test statistic based on the number of cases within \code{k}NNs of the cases
#' in the data.
#' 
#' For disease clustering, \insertCite{cuzick:1990;textual}{nnspat} suggested a \code{k}-NN test based on number of cases
#' among \code{k} NNs of the case points. 
#' Let \eqn{z_i} be the \eqn{i^{th}} point and \eqn{d_i^k} be the number cases among \code{k} NNs of \eqn{z_i}.
#' Then Cuzick-Edwards' \code{k}-NN test is \eqn{T_k=\sum_{i=1}^n \delta_i d_i^k}, where \eqn{\delta_i=1} 
#' if \eqn{z_i} is a case, and 0 if \eqn{z_i} is a control.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly.
#' Also, \eqn{T_1} is identical to the count for cell \eqn{(1,1)} in the nearest neighbor contingency table (NNCT)
#' (See the function \code{\link{nnct}} for more detail on NNCTs).
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014,cuzick:1990;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param k Integer specifying the number of NNs (of subject \eqn{i}), default is \code{1}.
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return Cuzick and Edwards \eqn{T_k} test statistic for disease clustering
#'
#' @seealso \code{\link{Tcomb}}, \code{\link{seg.ind}}, \code{\link{Pseg.coeff}} and \code{\link{ceTkinv}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#'
#' ceTk(Y,cls)
#' ceTk(Y,cls,method="max")
#' ceTk(Y,cls,k=3)
#' ceTk(Y,cls+1,case.lab = 2)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ceTk(Y,fcls,case.lab="a") #try also ceTk(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  # here ceTk(Y,cls) gives an error message
#' 
#' @export 
ceTk <- function(dat,cc.lab,k=1,case.lab=NULL,...)
{
  n<-nrow(dat)
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  if (length(table(cc.lab))!=2)
  {stop('cc.lab must have two types/labels for the case-control setting')}
  
  if (!is.null(case.lab))
  {lab2<-rep(0,length(cc.lab))
  lab2[cc.lab==case.lab]<-1
  cc.lab<-lab2}
  
  if (sum(cc.lab==1)==0)
  {stop('case-control labeling is incorrect, either label the cases as 1 and controls as 0 or specify the case.lab')}
  
  ipd<-ipd.mat(dat,...)
  Tk<-0
  nlist<-which(cc.lab==1) #indices of cases in cc.lab
  for (i in nlist)
  {
    Tk<-Tk+sum(cc.lab[order(ipd[i,])[2:(k+1)]])
  }
  Tk
} #end for the function
#'

#################################################################

# funsExpTk
#'
#' @title Expected Value for Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' Two functions: \code{EV.Tk} and \code{EV.Tkaij}.
#' 
#' Both functions compute the expected value of Cuzick and Edwards \eqn{T_k} test statistic based on the number of cases 
#' within \code{k}NNs of the cases in the data under RL or CSR independence.
#' 
#' The number of cases are denoted as \eqn{n_1} (denoted as \code{n1} as an argument)
#' for both functions and number of controls as \eqn{n_0} (denoted as \code{n0} as an argument) in \code{EV.Tk},
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The function \code{EV.Tkaij} uses Toshiro Tango's moments formulas based on the \eqn{A=(a_{ij})} matrix
#' (and is equivalent to the function \code{EV.Tk}, see \insertCite{tango:2007;textual}{nnspat},
#' where \eqn{a_{ij}(k) = 1} if \eqn{z_j} is among the \code{k}NNs of \eqn{z_i} and 0 otherwise.
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014;textual}{nnspat}).
#' 
#' @param k Integer specifying the number of NNs (of subject \eqn{i}).
#' @param n1,n0 The number of cases and controls, \eqn{n_1} used for both functions, and \eqn{n_0} used in \code{EV.Tk} only.
#' @param a The \eqn{A=(a_{ij})} matrix
#'  
#' @return The expected value of Cuzick and Edwards \eqn{T_k} test statistic for disease clustering
#'  
#' @seealso \code{\link{ceTk}} and \code{\link{EV.Tcomb}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsExpTk
NULL
#'
#' @rdname funsExpTk
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20
#' n0<-25
#' k<-1 #try also 3, 5, sample(1:5,1)
#'
#' EV.Tk(k,n1,n0)
#'
#' ###
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#' n1<-sum(cls==1)
#' n0<-sum(cls==0)
#' a<-aij.mat(Y,k)
#'
#' EV.Tk(k,n1,n0)
#' EV.Tkaij(k,n1,a)
#'
#' @export
EV.Tk <- function(k,n1,n0)
{
  n<-n1+n0
  k*n1*(n1-1)/(n-1)
} #end for the function
#'
#' @rdname funsExpTk
#'
#' @export
EV.Tkaij <- function(k,n1,a) 
{
  n<-nrow(a)
  p1<-n1*(n1-1)/(n*(n-1))
  S<-sum(a)
  p1*S
} #end for the function
#'

#################################################################

# funsAijmat
#'
#' @title Aij matrices for computation of Moments of Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' Two functions: \code{aij.mat} and \code{aij.nonzero}.
#' 
#' The function \code{aij.mat} yields the \eqn{A=(a_{ij}(k))} matrix where \eqn{a_{ij}(k) = 1} if \eqn{z_j} is among the \code{k}NNs of \eqn{z_i}
#' and 0 otherwise due to \insertCite{tango:2007;textual}{nnspat}.
#' This matrix is useful in calculation of the moments of Cuzick-Edwards \eqn{T_k} tests.
#'    
#' The function \code{aij.nonzero} keeps only nonzero entries, i.e., row and column entries where 
#' in each row, for the entry \eqn{(r_1,c_1)} \eqn{r_1} is the row entry and \eqn{c_1} is the column entry. Rows are from
#' 1 to n, which stands for the data point or observation, and column entries are from 1 to \code{k}, where \code{k} is specifying
#' the number of \code{k}NNs (of each observation) considered. This function saves in storage memory, but needs to be
#' carefully unfolded in the functions to represent the actual the \eqn{A} matrix.
#' 
#' See also (\insertCite{tango:2007;textual}{nnspat}).
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param k Integer specifying the number of NNs (of subject \eqn{i}), default is \code{1}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'  
#' @return The function \code{aij.mat} returns the \eqn{A_{ij}} matrix for computation of moments of Cuzick and Edwards \eqn{T_k}
#' Test statistic while the function \code{aij.nonzero} returns the (locations of the) non-zero entries in the \eqn{A_{ij}}
#' matrix
#'  
#' @seealso \code{\link{aij.theta}} and \code{\link{EV.Tkaij}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsAijmat
NULL
#'
#' @rdname funsAijmat
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' k<-3 #try also 2,3
#'
#' Aij<-aij.mat(Y,k)
#' Aij
#' Aij2<-aij.mat(Y,k,method="max")
#' range(Aij,Aij2)
#'
#' apply(Aij,2,sum) #row sums of Aij
#'
#' aij.nonzero(Y,k)
#' aij.nonzero(Y,k,method="max")
#'
#' @export
aij.mat <- function(dat,k,...)
{
  ipd<-ipd.mat(dat,...)
  n<-nrow(dat)
  a<-matrix(0,n,n)
  
  if (n<=1)
  { return(a)}
  
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  ord<-t(apply(ipd,1,order))
  
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      a[i,j]<-sum(j==ord[i,2:(k+1)]) #order(ipd[i,])[2:(k+1)])
    }
  }
  a
}#end for the function
#'
#' @rdname funsAijmat
#'
#' @export
aij.nonzero <- function(dat,k,...)
{
  ipd<-ipd.mat(dat,...)
  n<-nrow(dat)
  
  a<-matrix(0,n,n)
  
  if (n<=1)
  {return(NA)}
  
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  ord<-t(apply(ipd,1,order))
  
  cmat<-matrix(0,nrow=n,ncol=k) #for nonzero indices
  for (i in 1:n)
  {
    clist<-vector()
    for (j in 1:n)
    {
      alpha<-sum(j==ord[i,2:(k+1)]) #order(ipd[i,])[2:(k+1)])
      if (alpha==1)
      {
        clist<-c(clist,j)
      }
    }
    cmat[i,]<-clist
  }
  res<-cbind(1:n,cmat)
  col.nm<-"rlist"
  for (i in 1:k)
  {col.nm<-c(col.nm,paste("clist",i,sep=""))}
  colnames(res)<-col.nm
  res
} #end for the function
#'

#################################################################

#' @title Asymptotic Variance of Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' This function computes the asymptotic variance of Cuzick and Edwards \eqn{T_k} test statistic based on the number
#' of cases within \code{k}NNs of the cases in the data. 
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) for computing \eqn{N_s} and \eqn{N_t}, which are required in the computation of the
#' asymptotic variance. \eqn{N_s} and \eqn{N_t} are defined on page 78 of (\insertCite{cuzick:1990;textual}{nnspat}) as follows.
#' \eqn{N_s=\sum_i\sum_j a_{ij} a_{ji}} (i.e., number of ordered pairs for which \code{k}NN relation is symmetric)
#' and \eqn{N_t= \sum \sum_{i \ne l}\sum a_{ij} a_{lj}} (i.e, number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \code{k}NNs of \eqn{i} and \eqn{j} is among \code{k}NNs of \eqn{l}).
#' For the \eqn{A} matrix, see the description of the functions \code{aij.mat} and \code{aij.nonzero}.
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param n1 Number of cases
#' @param k Integer specifying the number of NNs (of subject \eqn{i})
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t}.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'   
#' @return A \code{list} with the elements
#' \item{asy.var}{The asymptotic variance of Cuzick and Edwards \eqn{T_k} test statistic for disease clustering}
#' \item{Ns}{The \eqn{N_s} value standing for the number of ordered pairs for which \code{k}NN relation is symmetric,
#' see the description.}
#' \item{Nt}{The \eqn{N_t} value standing for the number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \code{k}NNs of \eqn{i} and \eqn{j} is among \code{k}NNs of \eqn{l} see the description.}
#'
#' @seealso \code{\link{ceTk}}, \code{\link{varTk}}, and \code{\link{varTkaij}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#' k<-3 #try also 2,3
#'
#' asyvarTk(Y,n1,k)
#' asyvarTk(Y,n1,k,nonzero.mat=FALSE)
#' asyvarTk(Y,n1,k,method="max")
#' 
#' @export 
asyvarTk <- function(dat,n1,k,nonzero.mat=TRUE,...)
{
  n<-nrow(dat)
  
  if (n<=1)
  {return(NA)}
  
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  a<-aij.mat(dat,k,...) #needed for Ns and Nt
  
  if (nonzero.mat)
  {
    ak<-aij.nonzero(dat,k,...)
    row.ak<-ak[,1]; col.ak<-ak[,-1]
    
    Ns<-0
    for (i in 1:n)
    { for (el in 1:k)
    { ifelse(k==1, col.aki<-col.ak[i],col.aki<-col.ak[i,])
      Ns<-Ns+a[row.ak[i],col.aki[el]]*a[col.aki[el],row.ak[i]]
    }
    }
    
    Nt<-0
    for (j in 1:n)
    {
      tauj<-sum(a[,j])
      Nt<-Nt+tauj*(tauj-1)
    } 
  } else
  { 
    Ns<-Nt<-0 #for Ns and Nt see C&E 1990 page 78
    for (j in 1:n)
    {
      tauj<-sum(a[,j])
      for (i in 1:n)
      {
        Ns<-Ns+a[i,j]*a[j,i]
      }
      Nt<-Nt+tauj*(tauj-1)
    }
  }
  
  p1<-n1/n
  varTk<-n1*p1*(1-p1)*(k+Ns/n-p1*(k^2+Ns/n-Nt/n))
  
  list(
    asy.var=varTk,
    Ns=Ns,
    Nt=Nt)
} #end for the function
#'

#################################################################

#' @title \eqn{N_t} Value (found with the definition formula)
#'
#' @description
#' This function computes the \eqn{N_t} value which is required in the computation of the asymptotic variance
#' of Cuzick and Edwards \eqn{T_k} test. Nt is defined on page 78 of (\insertCite{cuzick:1990;textual}{nnspat}) as follows.
#' \eqn{N_t= \sum \sum_{i \ne l}\sum a_{ij} a_{lj}} (i.e, number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \eqn{k}NNs of \eqn{i} and \eqn{j} is among \eqn{k}NNs of \eqn{l}).
#' 
#' This function yields the same result as the \code{asyvarTk} and \code{varTk} functions with \code{$Nt} inserted at the
#' end.
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param a The \eqn{A=(a_{ij})} matrix. The argument \code{a} is the \eqn{A} matrix, obtained as output fromm \code{aij.mat}.
#'   
#' @return Returns the \eqn{N_t} value standing for the number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \eqn{k}NNs of \eqn{i} and \eqn{j} is among \eqn{k}NNs of \eqn{l}. See the description.
#'
#' @seealso \code{\link{asyvarTk}}, \code{\link{varTk}}, and \code{\link{varTkaij}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' k<-2 #try also 2,3
#' a<-aij.mat(Y,k)
#' Nt.def(a)
#' 
#' @export
Nt.def <- function(a)
{
  n<-nrow(a)
  Nt<-0
  lseq1<-1:n
  
  for (i in 1:n)
  {
    lseq<-lseq1[-i]
    for (l in lseq)
      for (j in 1:n)
        Nt<-Nt+a[i,j]*a[l,j]
  }
  Nt
} #end for the function
#'

#################################################################

#' @title \eqn{p_k} Value (used in the computation of the variance of Cuzick and Edwards \eqn{T_k} statistic)
#'
#' @description
#' This function computes the \eqn{p_k} value which is required in the computation of the variance
#' of Cuzick and Edwards \eqn{T_k} test. 
#' \eqn{p_k} is defined as the ratio \eqn{n_1(n_1-1)\cdots (n_1-(k-1))/(n(n-1)\cdots (n-(k-1))}.
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#'  
#' @param n A positive integer representing the number of points in the data set
#' @param n1 Number of cases
#' @param k Integer specifying the number of NNs (of subject \eqn{i})
#'   
#' @return Returns the \eqn{p_k} value. See the description.
#'
#' @seealso \code{\link{asyvarTk}}, \code{\link{varTk}}, and \code{\link{varTkaij}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#' 
pk <- function(n,n1,k)
{
  r<-prod(n1:(n1-(k-1)))/prod(n:(n-(k-1)))
  r
} #end for the function
#'

#################################################################

# funsVarTk
#'
#' @title Variance of Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' Two functions: \code{VarTk} and \code{VarTkaij}.
#' 
#' Both functions compute the (finite sample) variance of Cuzick and Edwards \eqn{T_k} test statistic based on the 
#' number of cases within \code{k}NNs of the cases in the data under RL or CSR independence.
#' 
#' The common arguments for both functions are \code{n1}, representing the number of cases and \code{k}.
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) for computing \eqn{N_s} and \eqn{N_t}, which are required in the computation of the
#' variance. \eqn{N_s} and \eqn{N_t} are defined on page 78 of (\insertCite{cuzick:1990;textual}{nnspat}) as follows.
#' \eqn{N_s=\sum_i\sum_j a_{ij} a_{ji}} (i.e., number of ordered pairs for which \code{k}NN relation is symmetric)
#' and \eqn{N_t= \sum \sum_{i \ne l}\sum a_{ij} a_{lj}} (i.e, number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \code{k}NNs of \eqn{i} and \eqn{j} is among \code{k}NNs of \eqn{l}).
#' 
#' The function \code{VarTkaij} uses Toshiro Tango's moments formulas based on the \eqn{A=(a_{ij})} matrix
#' (and is equivalent to the function \code{VarTk}, see \insertCite{tango:2007;textual}{nnspat},
#' where \eqn{a_{ij}(k) = 1} if \eqn{z_j} is among the \code{k}NNs of \eqn{z_i} and 0 otherwise.
#' 
#' The function \code{varTkaij} is equivalent to \code{varTk} (with \code{$var} extension).
#' 
#' See (\insertCite{cuzick:1990,tango:2007;textual}{nnspat}).
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point, used in \code{VarTk} only.
#' @param n1 Number of cases
#' @param k Integer specifying the number of NNs (of subject \eqn{i})
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t}.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used. Used in \code{VarTk} only.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. Used in \code{VarTk} only.
#' @param a The \eqn{A=(a_{ij})} matrix, used in \code{VarTkaij} only.
#' 
#' @return The function \code{VarTk} returns a \code{list} with the elements
#' \item{var.Tk}{The (finite sample) variance of Cuzick and Edwards \eqn{T_k} test statistic for disease clustering}
#' \item{Ns}{The \eqn{N_s} value standing for the number of ordered pairs for which \code{k}NN relation is symmetric,
#' see the description.}
#' \item{Nt}{The \eqn{N_t} value standing for the number of triplets \eqn{(i,j,l)} \eqn{i,j}, and \eqn{l} distinct so that
#' \eqn{j} is among \code{k}NNs of \eqn{i} and \eqn{j} is among \code{k}NNs of \eqn{l} see the description.}
#' 
#' The function \code{VarTkaij} returns only \code{var.Tk} as above.
#'  
#' @seealso \code{\link{asyvarTk}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsVarTk
NULL
#'
#' @rdname funsVarTk
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#' k<-2 #try also 2,3
#'
#' a<-aij.mat(Y,k)
#'
#' varTk(Y,n1,k)
#' varTk(Y,n1,k,nonzero.mat=FALSE)
#' varTk(Y,n1,k,method="max")
#' 
#' @export
varTk <- function(dat,n1,k,nonzero.mat=TRUE,...)
{ 
  n<-nrow(dat)
  
  if (n<=1)
  { return(list(
    var.Tk=NA,
    Ns=NA,
    Nt=NA))}
  
  if (k>=n)
  {k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  a<-aij.mat(dat,k,...)
  if (nonzero.mat)
    {
    ak<-aij.nonzero(dat,k,...)
    row.ak<-ak[,1]; col.ak<-ak[,-1]
    
    Ns<-0
    for (i in 1:n)
    { for (el in 1:k)
    { ifelse(k==1, col.aki<-col.ak[i],col.aki<-col.ak[i,])
      Ns<-Ns+a[row.ak[i],col.aki[el]]*a[col.aki[el],row.ak[i]]
    }
    }
    
    Nt<-0
    for (j in 1:n)
    {
      tauj<-sum(a[,j])
      Nt<-Nt+tauj*(tauj-1)
    }
  } else
  {
    Ns<-0
    Nt<-0
    for (j in 1:n)
    {
      tauj<-sum(a[,j])
      for (i in 1:n)
      {
        Ns<-Ns+a[i,j]*a[j,i]
      }
      Nt<-Nt+tauj*(tauj-1)
    }
  }
  
  p2<-pk(n,n1,2); p3<-pk(n,n1,3); p4<-pk(n,n1,4)
  
  var.Tk<-(k*n+Ns)*p2*(1-p2)+((3*k^2-k)*n+Nt-2*Ns)*(p3-p2^2)-(k^2*(n^2-3*n)+Ns-Nt)*(p2^2-p4)
  
  list(
    var.Tk=var.Tk,
    Ns=Ns,
    Nt=Nt)
} #end for the function
#'
#' @rdname funsVarTk
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#' k<-1 #try also 2,3, sample(1:5,1)
#'
#' a<-aij.mat(Y,k)
#'
#' varTkaij(n1,k,a)
#' varTk(Y,n1,k)$var
#' 
#' @export
varTkaij <- function(n1,k,a)
{
  n<-nrow(a); 
  if (n<=1)
  {return(NA)}
  
  if (k>=n)
  {k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  n0<-n-n1
  b<-1/2*(a+t(a))
  
  p2<-pk(n,n1,2); p3<-pk(n,n1,3); p4<-pk(n,n1,4)
  
  S1<-sum(b)
  S2<-2*sum(b^2)
  bi.<-apply(b,1,sum)
  S31<-sum(bi.*b) #this part avoids the two for loops 
  
  S3<-4*S31-2*S2
  S4<-S1^2-(S2+S3)
  
  ETk<-EV.Tk(k,n1,n0)
  var.Tk<-p2*S2+p3*S3+p4*S4-ETk^2
  var.Tk
} #end for the function
#'

#################################################################

#' @title \eqn{Z}-test for Cuzick and Edwards \eqn{T_k} statistic 
#' 
#' @description
#' An object of class \code{"htest"} performing a \eqn{z}-test for Cuzick and Edwards \eqn{T_k} test statistic based on the 
#' number of cases within \code{k}NNs of the cases in the data.
#' 
#' For disease clustering, \insertCite{cuzick:1990;textual}{nnspat} suggested a \code{k}-NN test \eqn{T_k} based on number of cases
#' among \code{k} NNs of the case points.
#' Under RL of \eqn{n_1} cases and \eqn{n_0} controls to the given locations in the study region,
#' \eqn{T_k} approximately has \eqn{N(E[T_k],Var[T_k]/n_1)} distribution for large \eqn{n_1}.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. 
#' Also, \eqn{T_1} is identical to the count for cell \eqn{(1,1)} in the nearest neighbor contingency table (NNCT)
#' (See the function \code{\link{nnct}} for more detail on NNCTs).
#' Thus, the \eqn{z}-test for \eqn{T_k} is same as the cell-specific \eqn{z}-test for cell \eqn{(1,1)} in the NNCT (see
#' \code{\link{cell.spec}}).
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#' 
#' The logical argument \code{asy.var} (default=\code{FALSE}) is for using the asymptotic variance or the exact (i.e. finite
#' sample) variance for the variance of \eqn{T_k} in its standardization.
#' If \code{asy.var=TRUE}, the asymptotic variance is used for \eqn{Var[T_k]} (see \code{asyvarTk}), otherwise the exact
#' variance (see \code{varTk}) is used.
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014,cuzick:1990;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param k Integer specifying the number of NNs (of subject \eqn{i}).
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for Cuzick and Edwards \eqn{T_k} statistic
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t} (argument is passed on to
#' \code{asyvarTk}). If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param asy.var A logical argument (default is \code{FALSE}) to determine whether the asymptotic variance or 
#' the exact (i.e. finite sample) variance for the variance of \eqn{T_k} in its standardization. 
#' If \code{TRUE}, the asymptotic variance is used for \eqn{Var[T_k]}, otherwise the exact variance is used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the Cuzick and Edwards \eqn{T_k} test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the Cuzick and Edwards \eqn{T_k} value
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the Cuzick and Edwards \eqn{T_k} value.}
#' \item{null.value}{Hypothesized null value for the Cuzick and Edwards \eqn{T_k} value
#' which is \eqn{k n_1 (n_1-1)/(n-1)} for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}}
#'  
#' @seealso \code{\link{ceTk}}, \code{\link{cell.spec}}, and \code{\link{Xsq.ceTk}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' k<-1 #try also 2,3, sample(1:5,1)
#'
#' ZceTk(Y,cls,k)
#' ZceTk(Y,cls,k,nonzero.mat=FALSE)
#' ZceTk(Y,cls,k,method="max")
#'
#' ZceTk(Y,cls+1,k,case.lab = 2,alt="l")
#' ZceTk(Y,cls,k,asy.var=TRUE,alt="g")
#' 
#' @export
ZceTk <- function(dat,cc.lab,k,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,
                  case.lab=NULL,nonzero.mat=TRUE,asy.var=FALSE,...)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  n<-nrow(dat)
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tk<-ceTk(dat,cc.lab,k,case.lab,...)
  ETk<-EV.Tk(k,n1,n0)
  ifelse(asy.var,var.Tk<-asyvarTk(dat,n1,k,nonzero.mat,...)$asy.var,
         var.Tk<-varTk(dat,n1,k,nonzero.mat,...)$var)
  stderr <-sqrt(var.Tk)
  ts<-(Tk-ETk)/stderr #TS.CE
  
  if (is.na(ts))
  {stop("The test statistic is NaN, so z-test for Cuzick-Edwards Tk is not defined")}
  
  names(ts) <-ifelse(asy.var==TRUE,paste("Test statistic for Cuzick-Edwards Tk for k=",k," (with asymptotic variance), Z ",sep=""),
                     paste("Test statistic for Cuzick-Edwards Tk for k=",k," (with finite sample variance), Z ",sep=""))
  
  estimate<-Tk
  names(estimate) <-c("Cuzick-Edwards Tk")
  method <-ifelse(asy.var==TRUE,paste("Z-Test for Cuzick-Edwards kNN test, Tk, with asymptotic variance ",sep=""),
                  paste("Z-Test for Cuzick-Edwards kNN test, Tk, with finite sample variance ",sep=""))
  
  null.val<- ETk
  names(null.val) <-"(expected) value of Cuzick-Edwards Tk under the null hypothesis"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )

  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(dat))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'

################################################
###FUNCTIONS for CUZICK and EDWARDS Tcomb TESTS###
################################################
#see page 87 of C&E 1990 for Tcomb, which combines multiple Tk's in one test

#' @title \eqn{N_{tkl}} Value
#' 
#' @description
#' This function computes the \eqn{N_{tkl}} value which is required in the computation of the exact and asymptotic variance
#' of Cuzick and Edwards \eqn{T_{comb}} test, which is a linear combination of some \eqn{T_k} tests. 
#' \eqn{N_{tkl}} is defined on page 80 of (\insertCite{cuzick:1990;textual}{nnspat}) as follows.
#' Let \eqn{a_{ij}(k)} be 1 if \eqn{j} is a \code{k} NN of \eqn{i} and zero otherwise and 
#' \eqn{N_t(k,l) = \sum \sum_{i \ne m}\sum a_{ij}(k) a_{mj}(l)}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#'   
#' See (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param k,l Integers specifying the number of NNs (of subjects \eqn{i} and \eqn{m} in \eqn{a_{ij}(k) a_{mj}(l)}).
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t} (argument is passed on to
#' \code{asycovTkTl} and \code{covTkTl}).
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'   
#' @return Returns the \eqn{N_{tkl}} value. See the description.
#' 
#' @seealso \code{\link{asycovTkTl}}, and \code{\link{covTkTl}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' k<-1 #try also 2,3 or sample(1:5,1)
#' l<-1 #try also 2,3 or sample(1:5,1)
#' c(k,l)
#'
#' Ntkl(Y,k,l)
#' Ntkl(Y,k,l,nonzero.mat = FALSE)
#' Ntkl(Y,k,l,method="max")
#' 
#' @export
Ntkl <- function(dat,k,l,nonzero.mat=TRUE,...)
{
  n<-nrow(dat)
  Nt<-0
  mseq1<-1:n
  if (nonzero.mat) #if nonzero.mat=T
  {
    ak<-aij.nonzero(dat,k,...)
    al<-aij.nonzero(dat,l,...)
    for (i in 1:n)
    {
      mseq<-mseq1[-i]
      for (m in mseq)
        Nt<-Nt+length(intersect(ak[i,-1],al[m,-1]))
    }
  } else #i.e., nonzero.mat=FALSE
  {
    ak<-aij.mat(dat,k,...)
    al<-aij.mat(dat,l,...)
    for (i in 1:n)
    {
      mseq<-mseq1[-i]
      for (m in mseq)
        Nt<-Nt+sum(ak[i,]*al[m,]) #this is short for "for (j in 1:n)"
    }  
  }
  Nt
} #end for the function
#'

#################################################################

#' @title Asymptotic Covariance between \eqn{T_k} and \eqn{T_l} Values
#' 
#' @description
#' This function computes the asymptotic covariance between \eqn{T_k} and \eqn{T_l} values
#' which is used in the computation of the asymptotic variance
#' of Cuzick and Edwards \eqn{T_{comb}} test, which is a linear combination of some \eqn{T_k} tests. 
#' The limit is as \eqn{n_1} goes to infinity.
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#'   
#' See page 80 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param n1 Number of cases
#' @param k,l Integers specifying the number of NNs (of subjects \eqn{i} and \eqn{m} in \eqn{a_{ij}(k) a_{mj}(l)}).
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t}.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'   
#' @return Returns the asymptotic covariance between \eqn{T_k} and \eqn{T_l} values.
#' 
#' @seealso \code{\link{covTkTl}}, \code{\link{covTcomb}}, and \code{\link{Ntkl}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#'
#' k<-1 #try also 2,3 or sample(1:5,1)
#' l<-1 #try also 2,3 or sample(1:5,1)
#' c(k,l)
#'
#' asycovTkTl(Y,n1,k,l)
#' asycovTkTl(Y,n1,k,l,nonzero.mat = FALSE)
#' asycovTkTl(Y,n1,k,l,method="max")
#' 
#' @export
asycovTkTl <- function(dat,n1,k,l,nonzero.mat=TRUE,...)
{
  n<-nrow(dat)
  
  if (n<=1)
  {return(NA)}
  
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  if (l>=n)
  { l<-min(l,n-1)
  warning("l is greater than or equal to data size, n, so n-1 is used for l")}
  
  Nt<-Ntkl(dat,k,l,nonzero.mat,...)
  
  if (nonzero.mat)
  { 
    ak<-aij.nonzero(dat,k,...)
    al<-aij.nonzero(dat,l,...)
    
    Ns<-0
    {for (i in 1:n)
      Ns<-Ns+sum(al[ak[i,-1],-1]==i) 
    }
  } else
  {
    ak<-aij.mat(dat,k,...)
    al<-aij.mat(dat,l,...)
    Ns<-0
    {for (i in 1:n)
      Ns<-Ns+sum(ak[i,]*al[,i]) #this is short for "for (j in 1:n) Ns<-Ns+ak[i,j]*al[j,i]"
    }
  }
  
  p1<-n1/n
  asycova<-n1*p1*(1-p1)*(k+Ns/n-p1*(k*l+Ns/n-Nt/n))
  
  asycova
} #end for the function
#'

#################################################################

#' @title Finite Sample Covariance between \eqn{T_k} and \eqn{T_l} Values
#' 
#' @description
#' This function computes the exact (i.e., finite sample) covariance between \eqn{T_k} and \eqn{T_l} values
#' which is used in the computation of the exact variance
#' of Cuzick and Edwards \eqn{T_{comb}} test, which is a linear combination of some \eqn{T_k} tests.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#'   
#' See page 80 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @inheritParams asycovTkTl
#'   
#' @return Returns the exact covariance between \eqn{T_k} and \eqn{T_l} values.
#' 
#' @seealso \code{\link{asycovTkTl}}, \code{\link{covTcomb}}, and \code{\link{Ntkl}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#'
#' k<-1 #try also 2,3 or sample(1:5,1)
#' l<-1 #try also 2,3 or sample(1:5,1)
#' c(k,l)
#'
#' covTkTl(Y,n1,k,l)
#' covTkTl(Y,n1,k,l,method="max")
#' asycovTkTl(Y,n1,k,l)
#'
#' covTkTl(Y,n1,k,l,nonzero.mat = FALSE)
#' asycovTkTl(Y,n1,k,l,nonzero.mat = FALSE)
#' 
#' @export
covTkTl <- function(dat,n1,k,l,nonzero.mat=TRUE,...)
{
  n<-nrow(dat)
  
  if (n<=1)
  {return(NA)}
  
  if (k>=n)
  { k<-min(k,n-1)
  warning("k is greater than or equal to data size, n, so n-1 is used for k")}
  
  if (l>=n)
  { l<-min(l,n-1)
  warning("l is greater than or equal to data size, n, so n-1 is used for l")}
  
  Nt<-Ntkl(dat,k,l,nonzero.mat,...)
  
  if (nonzero.mat)
  { 
    ak<-aij.nonzero(dat,k,...)
    al<-aij.nonzero(dat,l,...)
    
    Ns<-0
    {for (i in 1:n)
      Ns<-Ns+sum(al[ak[i,-1],-1]==i) 
    }
  } else
  {
    ak<-aij.mat(dat,k,...)
    al<-aij.mat(dat,l,...)
    Ns<-0
    {for (i in 1:n)
      Ns<-Ns+sum(ak[i,]*al[,i]) #this is short for "for (j in 1:n) Ns<-Ns+ak[i,j]*al[j,i]"
    }
  }
  
  p2<-n1*(n1-1)/(n*(n-1))
  p3<-n1*(n1-1)*(n1-2)/(n*(n-1)*(n-2))
  p4<-n1*(n1-1)*(n1-2)*(n1-3)/(n*(n-1)*(n-2)*(n-3))
  
  cova<-(k*n+Ns)*p2*(1-p2)+((3*k*l-k)*n+Nt-2*Ns)*(p3-p2^2)-(k*l*(n^2-3*n)+Ns-Nt)*(p2^2-p4)
  
  cova
} #end for the function
#'

#################################################################

#' @title Covariance matrix for \eqn{T_k} values in \code{Tcomb}
#' 
#' @description
#' This function computes the covariance matrix for the \eqn{T_k} values used in the \eqn{T_{comb}} test statistics,
#' which is a linear combination of some \eqn{T_k} tests. 
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases is denoted as \eqn{n_1} to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The argument \code{klist} is the \code{vector} of integers specifying the indices of the \eqn{T_k} values used
#' in obtaining the \eqn{T_{comb}}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#' 
#' The logical argument \code{asy.cov} (default=\code{FALSE}) is for using the asymptotic covariance or the exact (i.e. finite
#' sample) covariance for the vector of \eqn{T_k} values used in \code{Tcomb}.
#' If \code{asy.cov=TRUE}, the asymptotic covariance is used, otherwise the exact covariance is used.
#'   
#' See page 87 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param n1 Number of cases
#' @param klist \code{list} of integers specifying the indices of the \eqn{T_k} values used in obtaining the \eqn{T_{comb}}.
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t}.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param asy.cov A logical argument (default is \code{FALSE}) to determine whether asymptotic or exact (i.e., finite
#' sample) covariances between \eqn{T_k} and \eqn{T_l} values are to be used to obtain the entries of the covariance matrix.
#' If \code{TRUE} the asymptotic covariance values are used, otherwise exact covariance values are used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'   
#' @return Returns the covariance matrix for the \eqn{T_k} values used in \code{Tcomb}.
#' 
#' @seealso \code{\link{asycovTkTl}}, \code{\link{covTcomb}}, and \code{\link{Ntkl}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#'
#' kl<-sample(1:5,3) #try also sample(1:5,2)
#' kl
#' covTcomb(Y,n1,kl)
#' covTcomb(Y,n1,kl,method="max")
#' covTcomb(Y,n1,kl,nonzero.mat = FALSE)
#'
#' covTcomb(Y,n1,kl,asy=TRUE)
#' 
#' @export 
covTcomb <- function(dat,n1,klist,nonzero.mat=TRUE,asy.cov=FALSE,...) 
{
  n<-nrow(dat)
  if (n<=1)
  {return(NA)}
  
  if (max(klist)>=n)
  { klist<-unique(c(klist[klist<=n-1],n-1))
  warning("klist has elements greater than or equal to data size, n, so klist= {", toString(klist),"} is used instead")}
  
  r<-length(klist) #must be renewed when max(klist)>=n
  
  varTc<-vector()
  for (k in klist)
  {
    ifelse(asy.cov,VarTc<-asyvarTk(dat,n1,k,nonzero.mat)$asy.var,
           VarTc<-varTk(dat,n1,k,nonzero.mat)$var)
    varTc<-c(varTc,VarTc)
  }
  
  if (r==1)
  {return(varTc)}
  
  sigma<-matrix(0,r,r)
  for (i in 1:(r-1))
    for (j in (i+1):r)
    {
      k<-klist[i]; l<-klist[j]
      ifelse(asy.cov==TRUE,sigma[i,j]<-asycovTkTl(dat,n1,k,l,nonzero.mat),#,...),
             sigma[i,j]<-covTkTl(dat,n1,k,l,nonzero.mat))#,...))
      sigma[j,i]<-sigma[i,j]
    }
  
  diag(sigma)<-varTc
  sigma
} #end for the function
#'

#################################################################

#' @title Square root of a matrix
#' 
#' @description
#' Computes the square root of the matrix \eqn{A}, where \eqn{A} does not have to be a square matrix, 
#' when the square root exists.
#' See https://people.orie.cornell.edu/davidr/SDAFE2/Rscripts/SDAFE2.R
#' 
#' @param A A matrix, not necessarily square
#'   
#' @return Returns the square root of \eqn{A}, if exists, otherwise gives an error message.
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' A<-matrix(sample(20:40,4),ncol=2)
#' matrix.sqrt(A)
#'
#' A<-matrix(sample(20:40,16),ncol=4)
#' matrix.sqrt(A)
#' #sqrt of inverse of A, or sqrt inverse of A
#' matrix.sqrt(solve(A))
#'
#' #non-square matrix
#' A<-matrix(sample(20:40,20),ncol=4)
#' matrix.sqrt(A)
#' 
#' @export
matrix.sqrt <- function(A) 
{ 
  sva <- svd(A) 
  if (min(sva$d)>=0) 
  { Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d))) 
  } else 
  {stop("Matrix square root is not defined") }
  return(Asqrt) 
} #end for the function
#'

#################################################################

#' @title Cuzick & Edwards Tcomb Test Statistic
#' 
#' @description
#' This function computes the value of Cuzick & Edwards \eqn{T_{comb}} test statistic in disease clustering, where \eqn{T_{comb}}
#' is a linear combination of some \eqn{T_k} tests.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. 
#' 
#' The argument \code{klist} is the \code{vector} of integers specifying the indices of the \eqn{T_k} values used
#' in obtaining the \eqn{T_{comb}}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#' 
#' The logical argument \code{asy.cov} (default=\code{FALSE}) is for using the asymptotic covariance or the exact (i.e. finite
#' sample) covariance for the vector of \eqn{T_k} values used in \code{Tcomb} in the standardization of \eqn{T_{comb}}.
#' If \code{asy.cov=TRUE}, the asymptotic covariance is used, otherwise the exact covariance is used. 
#'   
#' See page 87 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param klist \code{list} of integers specifying the indices of the \eqn{T_k} values used in obtaining the \eqn{T_{comb}}.
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of \eqn{N_s} and \eqn{N_t}.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param asy.cov A logical argument (default is \code{FALSE}) to determine whether asymptotic or exact (i.e., finite
#' sample) covariances between \eqn{T_k} and \eqn{T_l} values are to be used to obtain the entries of the covariance matrix.
#' If \code{TRUE} the asymptotic covariance values are used, otherwise exact covariance values are used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'   
#' @return Returns the value of the \eqn{T_{comb}} test statistic
#' 
#' @seealso \code{\link{ceTk}}, \code{\link{EV.Tcomb}}, and \code{\link{ZTcomb}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1) #try also n<-50, 100
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#'
#' kl<-sample(1:5,3) #try also sample(1:5,2)
#' kl
#' Tcomb(Y,cls,kl)
#' Tcomb(Y,cls,kl,method="max")
#' Tcomb(Y,cls+1,kl,case.lab=2)
#' Tcomb(Y,cls,kl,nonzero.mat = FALSE)
#' Tcomb(Y,cls,kl,asy.cov = TRUE)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' Tcomb(Y,fcls,kl,case.lab="a")
#' 
#' @export 
Tcomb <- function(dat,cc.lab,klist,case.lab=NULL,nonzero.mat=TRUE,asy.cov=FALSE,...)
{
  n<-nrow(dat)
  
  if (max(klist)>=n)
  { klist<-unique(c(klist[klist<=n-1],n-1))
  warning("klist has elements greater than or equal to data size, n, so klist= {", toString(klist),"} is used instead")}
  
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tv<-vector()
  for (k in klist)
  {Tv<-c(Tv,ceTk(dat,cc.lab,k,case.lab,...))}
  
  ifelse(asy.cov==TRUE,Sigma<-covTcomb(dat,n1,klist,nonzero.mat,asy.cov = TRUE,...) ,
         Sigma<-covTcomb(dat,n1,klist,nonzero.mat,...) ) #cov matrix for Tcomb
  Sinvrt<-matrix.sqrt(solve(Sigma)) #sqrt of inverse of Sigma
  
  Tcomb<-sum(Sinvrt %*% Tv)
  Tcomb
} #end for the function
#'

#################################################################

#' @title Expected Value for Cuzick & Edwards \eqn{T_{comb}} Test Statistic
#' 
#' @description
#' This function computes the expected value of Cuzick & Edwards \eqn{T_{comb}} test statistic in disease clustering,
#' where \eqn{T_{comb}} is a linear combination of some \eqn{T_k} tests. 
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases is denoted as \eqn{n_1} to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The argument \code{klist} is the \code{vector} of integers specifying the indices of the \eqn{T_k} values used
#' in obtaining the \eqn{T_{comb}}.
#' 
#' The argument \code{sig} is the covariance matrix of the vector of \eqn{T_k} values used in \code{Tcomb}, and can be computed
#' via the the \code{\link{covTcomb}} function.
#'   
#' See page 87 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' @param n1 Number of cases
#' @param n A positive integer representing the number of points in the data set
#' @param klist \code{list} of integers specifying the indices of the \eqn{T_k} values used in obtaining the \eqn{T_{comb}}.
#' @param sig The covariance matrix of the vector of \eqn{T_k} values used in \code{Tcomb}
#'   
#' @return Returns the expected value of the \eqn{T_{comb}} test statistic
#' 
#' @seealso \code{\link{Tcomb}}, and \code{\link{ZTcomb}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1) #try also n<-50, 100
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' n1<-sum(cls==1)
#'
#' kl<-sample(1:5,3) #try also sample(1:5,2)
#' kl
#' sig<-covTcomb(Y,n1,kl)
#' EV.Tcomb(n1,n,kl,sig)
#' 
#' @export 
EV.Tcomb <- function(n1,n,klist,sig)
{
  if (max(klist)>=n)
  { klist<-unique(c(klist[klist<=n-1],n-1))
  warning("klist has elements greater than or equal to data size, n, so klist= {", toString(klist),"} is used instead")}
  
  p2<-n1*(n1-1)/(n*(n-1))
  sig.inv.rt<-matrix.sqrt(solve(sig)) #sqrt of inverse of Sigma
  EV<-n*p2*sum(klist*sig.inv.rt) #this is incorrect in C&E, and fixed here
  EV
} #end for the function
#'

#################################################################

#' @title \eqn{Z}-test for Cuzick and Edwards \eqn{T_{comb}} statistic 
#' 
#' @description
#' An object of class \code{"htest"} performing a \eqn{z}-test for Cuzick and Edwards \eqn{T_{comb}} test statisticin disease clustering,
#' where \eqn{T_{comb}} is a linear combination of some \eqn{T_k} tests. 
#' 
#' For disease clustering, \insertCite{cuzick:1990;textual}{nnspat} developed a \eqn{k}-NN test \eqn{T_k} based on 
#' number of cases among \eqn{k} NNs of the case points, and also proposed a test combining various \eqn{T_k} tests,
#' denoted as \eqn{T_{comb}}.
#' 
#' See page 87 of (\insertCite{cuzick:1990;textual}{nnspat}) for more details.
#' 
#' Under RL of \eqn{n_1} cases and \eqn{n_0} controls to the given locations in the study region,
#' \eqn{T_{comb}} approximately has \eqn{N(E[T_{comb}],Var[T_{comb}])} distribution for large \eqn{n_1}.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. 
#' 
#' The argument \code{klist} is the \code{vector} of integers specifying the indices of the \eqn{T_k} values used
#' in obtaining the \eqn{T_{comb}}.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{TRUE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}) in the computations.
#' 
#' The logical argument \code{asy.cov} (default=\code{FALSE}) is for using the asymptotic covariance or the exact (i.e. finite
#' sample) covariance for the vector of \eqn{T_k} values used in \code{Tcomb} in the standardization of \eqn{T_{comb}}.
#' If \code{asy.cov=TRUE}, the asymptotic covariance is used, otherwise the exact covariance is used.
#' 
#' See also (\insertCite{ceyhan:SiM-seg-ind2014,cuzick:1990;textual}{nnspat})
#' and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param klist \code{list} of integers specifying the indices of the \eqn{T_k} values used in obtaining the \eqn{T_{comb}}.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for Cuzick and Edwards \eqn{T_{comb}} statistic
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computation of covariance of \eqn{T_k} values forming the
#' \code{T_{comb}} statistic (argument is passed on to \code{covTcomb}). If \code{TRUE} the nonzero location matrix is used,
#' otherwise the \eqn{A} matrix itself is used.
#' @param asy.cov A logical argument (default is \code{FALSE}) to determine whether asymptotic or exact (i.e., finite
#' sample) covariances between \eqn{T_k} and \eqn{T_l} values are to be used to obtain the entries of the covariance matrix.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the Cuzick and Edwards \eqn{T_{comb}} test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the Cuzick and Edwards \eqn{T_{comb}} value
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the Cuzick and Edwards \eqn{T_{comb}} value.}
#' \item{null.value}{Hypothesized null value for the Cuzick and Edwards \eqn{T_{comb}} value
#' which is \eqn{E[T_{comb}]} for this function, which is the output of \code{EV.Tcomb} function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}}
#'  
#' @seealso \code{\link{Tcomb}}, \code{\link{EV.Tcomb}}, and \code{\link{covTcomb}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#'
#' kl<-sample(1:5,3) #try also sample(1:5,2)
#' ZTcomb(Y,cls,kl)
#' ZTcomb(Y,cls,kl,method="max")
#'
#' ZTcomb(Y,cls,kl,nonzero.mat=FALSE)
#' ZTcomb(Y,cls+1,kl,case.lab = 2,alt="l")
#' ZTcomb(Y,cls,kl,conf=.9,alt="g")
#' ZTcomb(Y,cls,kl,asy=TRUE,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ZTcomb(Y,fcls,kl,case.lab="a")
#' 
#' @export
ZTcomb <- function(dat,cc.lab,klist,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,
                   case.lab=NULL,nonzero.mat=TRUE,asy.cov=FALSE,...)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  n<-nrow(dat)
  klist0<-klist
  if (max(klist)>=n)
  { klist<-unique(c(klist[klist<=n-1],n-1))
  warning("klist has elements greater than or equal to data size, n, so klist= {", toString(klist),"} is used instead")}
  
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tv<-vector()
  for (k in klist)
  {Tv<-c(Tv,ceTk(dat,cc.lab,k,case.lab,...))}
  
  ifelse(asy.cov==TRUE,Sigma<-covTcomb(dat,n1,klist,nonzero.mat,asy.cov = TRUE,...) ,
         Sigma<-covTcomb(dat,n1,klist,nonzero.mat,...) ) #cov matrix for Tcomb
  Sinvrt<-matrix.sqrt(solve(Sigma)) #sqrt of inverse of Sigma
  p2<-n1*(n1-1)/(n*(n-1))
  ETcomb<-n*p2*sum(klist*Sinvrt) #or use EV.Tcomb(n1,n,klist,Sigma)
  
  Tcomb<-sum(Sinvrt %*% Tv)
  m<-length(klist) #also Var[Tcomb]
  stderr<-sqrt(m)
  ts<-(Tcomb-ETcomb)/stderr
  names(ts) <-ifelse(asy.cov==TRUE,
                     "Test statistic for Cuzick-Edwards Tcomb which combines Tk for k in klist with asymptotic covariance, Z",
                     "Test statistic for Cuzick-Edwards Tcomb which combines Tk for k in klist with finite sample covariance, Z")
  
  estimate<-Tcomb
  names(estimate) <-c("Cuzick-Edwards Tcomb")
  
  method <- ifelse(asy.cov==TRUE,
                   paste(c("Z-Test for Cuzick-Edwards combined Tk for k in klist={", klist0,"} with asymptotic covariance"), collapse=" "),
                   paste(c("Z-Test for Cuzick-Edwards combined Tk for k in klist={", klist0,"} with finite sample covariance"), collapse=" "))
  
  null.val<- ETcomb
  names(null.val) <-"(expected) value of Cuzick-Edwards combined Tk under the null hypothesis"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(dat))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'

################################################ 
###FUNCTIONS for CUZICK and EDWARDS Trun TESTS###
################################################

#' @title Cuzick and Edwards \eqn{T_{run}} Test statistic
#'
#' @description
#' This function computes Cuzick and Edwards \eqn{T_{run}} test statistic based on the sum of the number of successive 
#' cases from each cases until a control is encountered in the data for detecting rare large clusters.
#' 
#' \eqn{T_{run}} test statistic is defined as \eqn{T_{run}=\sum_{i=1}^n \delta_i d_i^r} where \eqn{\delta_i=1} 
#' if \eqn{z_i} is a case, and 0 if \eqn{z_i} is a control and \eqn{d_i^r} is the number successive cases encountered beginning
#' at \eqn{z_i} until a control is encountered. 
#'  
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL},
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's and 1's
#' accordingly.
#' 
#' See also (\insertCite{cuzick:1990;textual}{nnspat}) and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'
#' @return A \code{list} with two elements
#' \item{Trun}{Cuzick and Edwards \eqn{T_{run}} test statistic for disease clustering}
#' \item{run.vec}{The \code{vector} of number of consecutive cases till the first control for each point in the data set}
#'
#' @seealso \code{\link{ceTk}}, \code{\link{Tcomb}} and \code{\link{ceTkinv}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#'
#' ceTrun(Y,cls)
#' ceTrun(Y,cls,method="max")
#' ceTrun(Y,cls+1,case.lab = 2)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ceTrun(Y,fcls,case.lab="a") #try also ceTrun(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #here ceTrun(Y,cls) #gives an error message
#' 
#' @export 
ceTrun <- function(dat,cc.lab,case.lab=NULL,...)
{
  if (length(table(cc.lab))!=2)
  {stop('cc.lab must have two types/labels for the case-control setting')}
  
  if (!is.null(case.lab))
  {
    lab2<-rep(0,length(cc.lab))
    lab2[cc.lab==case.lab]<-1
    cc.lab<-lab2}
  
  if (sum(cc.lab==1)==0)
  {stop('case-control labeling is incorrect, either label the cases as 1 and controls as 0 or specify the case.lab')}
  
  ipd<-ipd.mat(dat,...)
  n<-nrow(dat)
  Trun<-0; rv<-rep(0,n) #runvector
  nlist<-which(cc.lab==1) #indices of cases in cc.lab
  for (i in nlist)
  {
    crit<-1; k<-1
    while (crit>0)
    {
      crit<-prod(cc.lab[order(ipd[i,])[2:(k+1)]])
      k<-k+1
    }
    Trun<-Trun+(k-2) #k-2 to exclude the original case point and the first ctrl reached
    rv[i]<-k-2
  }
  list(Trun=Trun,
       run.vec=rv
  )
} #end for the function
#'

#################################################################

# funsExpTrun
#'
#' @title Expected Value for Cuzick and Edwards \eqn{T_{run}} Test statistic
#'
#' @description
#' Two functions: \code{EV.Trun} and \code{EV.Trun.alt}.
#' 
#' Both functions compute the expected value of Cuzick and Edwards \eqn{T_{run}} test statistic based on the number of 
#' consecutive cases from the cases in the data under RL or CSR independence.
#' 
#' The number of cases are denoted as \eqn{n_1} (denoted as \code{n1} as an argument)
#' and number of controls as \eqn{n_0} for both functions (denoted as \code{n0} as an argument),
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The function \code{EV.Trun.alt} uses a loop and takes slightly longer than the function \code{EV.Trun},
#' hence \code{EV.Trun} is used in other functions. 
#' 
#' See also (\insertCite{cuzick:1990;textual}{nnspat}).
#' 
#' @param n1,n0 The number of cases and controls used as arguments for both functions.
#'  
#' @return The expected value of Cuzick and Edwards \eqn{T_{run}} test statistic for disease clustering
#'  
#' @seealso \code{\link{ceTrun}} and \code{\link{EV.Tk}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsExpTrun
NULL
#'
#' @rdname funsExpTrun
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20
#' n0<-25
#'
#' EV.Trun(n1,n0)
#'
#' @export
EV.Trun <- function(n1,n0)
{
  n1*(n1-1)/(n0+1)
} #end for the function
#'
#' @rdname funsExpTrun
#'
EV.Trun.alt <- function(n1,n0)
{
  ev<-0
  n<-n1+n0
  for (j in 1:(n1-1))
  {
    pj<-prod(n1:(n1-j))/prod(n:(n-j))
    ev<-ev+j*pj*n0/(n-j-1)
  }
  ev*n
} #end for the function
#'

#################################################################

# funsVarTrun
#'
#' @title Variance of Cuzick and Edwards \eqn{T_{run}} Test statistic
#'
#' @description
#' Two functions: \code{varTrun} and \code{varTrun.sim}.
#' 
#' The function \code{varTrun} computes the (finite sample) variance of Cuzick and Edwards \eqn{T_{run}} test statistic 
#' which is based on the number of consecutive cases from the cases in the data under RL or CSR independence.
#' And the function \code{varTrun.sim} estimates this variance based on simulations under the RL hypothesis.
#' 
#' The only common argument for both functions is \code{dat}, the data set used in the functions.
#' 
#' \eqn{n_1} is an argument for \code{varTrun} and is the number of cases (denoted as \code{n1} as an argument).
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. The argument \code{Nsim} represents the number of resamplings (without replacement) in the
#' RL scheme, with default being \code{1000}. \code{cc.lab}, \code{case.lab} and \code{Nsim} are arguments for \code{varTrun.sim} only.
#' 
#' The function \code{varTrun} might take a very long time when data size is large (even larger than 50),
#' hence the need for the \code{varTrun.sim} function. 
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}).
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point, 
#' used in both functions.
#' @param n1 Number of cases, used in \code{varTrun} only.
#' @param cc.lab Case-control labels, 1 for case, 0 for control, used in \code{varTrun.sim} only.
#' @param Nsim The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
#' variance of \eqn{T_{run}}, used in \code{varTrun.sim} only.
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}, used in \code{varTrun.sim} only.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' Used in \code{varTrun} only.
#' 
#' @return The function \code{varTrun} returns the variance of Cuzick and Edwards \eqn{T_{run}} test statistic
#' under RL or CSR independence.
#' And the function \code{varTrun.sim} estimates the same variance based on simulations under the RL hypothesis.
#' 
#' @seealso \code{\link{ceTrun}} and \code{\link{EV.Trun}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsVarTrun
NULL
#'
#' @rdname funsVarTrun
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1) #try also 40, 50, 60
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#' n1<-sum(cls==1)
#' n0<-sum(cls==0)
#' c(n1,n0)
#'
#' varTrun(Y,n1)
#' varTrun(Y,n1,method="max")
#' 
#' @export
varTrun <- function(dat,n1,...)
{
  n<-nrow(dat)
  if (n<=1)
  {return(NA)}
  
  n0<-n-n1
  ipd<-ipd.mat(dat,...)
  ord<-t(apply(ipd,1,order))
  S<-0; 
  for (i in 1:n)
  {for (j in (1:n))
  {
    for (k in 1:(n1-1))
    {
      set1<-ord[i,1:(k+1)]
      for (l in 1:(n1-1))
      {
        set2<-ord[j,1:(l+1)]
        nijkl<-length(intersect(set1,set2))
        jnt<-(1+k+l-nijkl)
        pjnt<-prod(n1:(n1-jnt))/prod(n:(n-jnt))
        S<-S+pjnt 
      }
    }
  }
  }
  ev<-EV.Trun(n1,n0)
  S<-S-ev^2
  S
} #end for the function
#'
#' @rdname funsVarTrun
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1) #try also 40, 50, 60
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#' n1<-sum(cls==1)
#' varTrun(Y,n1) #the actual value (might take a long time if \eqn{n} is large)
#'
#' Nmc<-1000
#' varTrun.sim(Y,cls,Nsim=Nmc)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' varTrun.sim(Y,fcls,Nsim=Nmc,case.lab="a")
#' 
#' @export
varTrun.sim <- function(dat,cc.lab,Nsim=1000,case.lab=NULL)
{
  if (length(table(cc.lab))!=2)
  {stop('cc.lab must have two types/labels for the case-control setting')}
  
  if (!is.null(case.lab))
  {lab2<-rep(0,length(cc.lab))
  lab2[cc.lab==case.lab]<-1
  cc.lab<-lab2}
  
  ce.vec<-vector()
  for (i in 1:Nsim)
  {
    slab<-sample(cc.lab) #random labeling of cases and ctrls
    ce<-ceTrun(dat,slab)
    ce.vec<-c(ce.vec,ce$Tr)
  }
  var(ce.vec)
} #end for the function
#'

#################################################################

#' @title \eqn{Z}-test for Cuzick and Edwards \eqn{T_{run}} statistic 
#' 
#' @description
#' An object of class \code{"htest"} performing a \eqn{z}-test for Cuzick and Edwards \eqn{T_{run}} test statistic 
#' which is based on the number of consecutive cases from the cases in the data under RL or CSR independence.
#' 
#' Under RL of \eqn{n_1} cases and \eqn{n_0} controls to the given locations in the study region,
#' \eqn{T_{run}} approximately has \eqn{N(E[T_{run}],Var[T_{run}])} distribution for large \eqn{n}.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. 
#' 
#' The logical argument var.sim (default=\code{FALSE}) is for using the simulation estimated variance or the exact 
#' variance for the variance of \eqn{T_{run}} in its standardization.
#' If \code{var.sim=TRUE}, the simulation estimated variance is used for \eqn{Var[T_{run}]} (see \code{varTrun.sim}), 
#' otherwise the exact variance (see \code{varTrun}) is used.
#' Moreover, when \code{var.sim=TRUE}, the argument \code{Nvar.sim} represents the number of resamplings 
#' (without replacement) in the RL scheme, with default being \code{1000}.
#' 
#' The function \code{varTrun} might take a very long time when data size is large (even larger than 50);
#' in this case, it is recommended to use \code{var.sim=TRUE} in this function.
#'  
#' See also (\insertCite{cuzick:1990;textual}{nnspat}) and the references therein.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for Cuzick and Edwards \eqn{T_{run}} statistic
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param var.sim A logical argument (default is \code{FALSE}) to determine whether the simulation estimated variance or
#' the exact variance be used for the variance of \eqn{T_{run}} in its standardization.
#' If \code{var.sim=TRUE}, the simulation estimated variance is used for \eqn{Var[T_{run}]} (see \code{varTrun.sim}), 
#' otherwise the exact variance (see \code{varTrun}) is used.
#' @param Nvar.sim The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
#' variance of \eqn{T_{run}}, used only when \code{var.sim=TRUE}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#'  
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the Cuzick and Edwards \eqn{T_{run}} test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative}
#' \item{conf.int}{Confidence interval for the Cuzick and Edwards \eqn{T_{run}} value
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.} 
#' \item{estimate}{Estimate of the parameter, i.e., the Cuzick and Edwards \eqn{T_{run}} value.}
#' \item{null.value}{Hypothesized null value for the Cuzick and Edwards \eqn{T_{run}} value
#' which is \eqn{n_1 (n_1-1)/(n_0+1)} for this function.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}}
#'  
#' @seealso \code{\link{ceTrun}}, \code{\link{ZceTk}}, and \code{\link{ZTcomb}}
#' 
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1) #try also 40, 50, 60
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#'
#' ZTrun(Y,cls)
#' ZTrun(Y,cls,method="max")
#' ZTrun(Y,cls,var.sim=TRUE)
#' ZTrun(Y,cls+1,case.lab = 2,alt="l")
#' ZTrun(Y,cls,conf=.9,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ZTrun(Y,fcls,case.lab="a")
#' 
#' @export
ZTrun <- function(dat,cc.lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,
                  case.lab=NULL,var.sim=FALSE,Nvar.sim=1000,...)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  n<-nrow(dat)
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tr<-ceTrun(dat,cc.lab,case.lab,...)$Tr
  
  ifelse(var.sim==TRUE,varTr<-varTrun.sim(dat,cc.lab,Nsim=Nvar.sim,case.lab),
         varTr<-varTrun(dat,n1,...) ) #var of Trun
  ETr<-EV.Trun(n1,n0)
  
  stderr<-sqrt(varTr)
  ts<-(Tr-ETr)/stderr
  names(ts) <-ifelse(var.sim==TRUE,
                     "Test statistic for Cuzick-Edwards Trun with variance estimated via simulation, Z ",
                     "Test statistic for Cuzick-Edwards Trun with the actual variance, Z ")
  
  estimate<-Tr
  names(estimate) <-c("Cuzick-Edwards Trun")
  
  method <- ifelse(var.sim==TRUE,
                   "Z-Test for Cuzick-Edwards Trun for clustering with variance estimated via simulation",
                   "Z-test for Cuzick-Edwards Trun for clustering with the actual variance")
  
  null.val<- ETr
  names(null.val) <-"(expected) value of Cuzick-Edwards Trun under the null hypothesis"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(dat))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'

################################################
###FUNCTIONS for CUZICK and EDWARDS Tkinv TESTS###
################################################

#' @title Cuzick and Edwards \eqn{T_k^{inv}} Test statistic
#'
#' @description
#' This function computes Cuzick and Edwards \eqn{T_k^{inv}} test statistic based on the sum of number of cases closer to 
#' each case than the \code{k}-th nearest control to the case.
#' 
#' \eqn{T_k^{inv}} test statistic is an extension of the run length test allowing a fixed number of controls in the run 
#' sequence. 
#' 
#' \eqn{T_k^{inv}} test statistic is defined as \eqn{T_k^{inv}=\sum_{i=1}^n \delta_i \nu_i^k} where \eqn{\delta_i=1} 
#' if \eqn{z_i} is a case, and 0 if \eqn{z_i} is a control and \eqn{\nu_i^k} is the number of cases closer
#' to the index case than the \code{k} nearest control, i.e., number of cases encountered beginning
#' at \eqn{z_i} until \code{k}-th control is encountered. 
#'  
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL},
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's and 1's
#' accordingly.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param k Integer specifying the number of the closest controls to subject \eqn{i}.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return A \code{list} with two elements
#' \item{Tkinv}{Cuzick and Edwards \eqn{T_k^{inv}} test statistic for disease clustering}
#' \item{run.vec}{The \code{vector} of number of cases till the \code{k}-th control for each point in the data set}
#'
#' @seealso \code{\link{ceTrun}}, \code{\link{ceTk}}, and \code{\link{Tcomb}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' cls
#' k<-2 #also try 3,4
#'
#' ceTkinv(Y,k,cls)
#' ceTkinv(Y,k,cls+1,case.lab = 2)
#' ceTkinv(Y,k,cls,method="max")
#'
#' ceTrun(Y,cls)
#' ceTkinv(Y,k=1,cls)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ceTkinv(Y,k,fcls,case.lab="a") #try also ceTrun(Y,fcls)
#'
#' #############
#' n<-40
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(1:4,n,replace = TRUE)  #here ceTkinv(Y,k,cls) #gives error
#' 
#' @export 
ceTkinv <- function(dat,k,cc.lab,case.lab=NULL,...)
{
  if (length(table(cc.lab))!=2)
  {stop('cc.lab must have two types/labels for the case-control setting')}
  
  if (!is.null(case.lab))
  {lab2<-rep(0,length(cc.lab))
  lab2[cc.lab==case.lab]<-1
  cc.lab<-lab2}
  
  n0<-sum(cc.lab==0)
  if (n0==0)
  {stop('case-control labeling is incorrect, either label the cases as 1 and controls as 0 or specify the case.lab')}
  
  if (k>n0)
  {stop('k can not be larger than the number of controls in the sample')}
  
  ipd<-ipd.mat(dat,...)
  n<-nrow(dat)
  Tkinv<-0; rv<-rep(0,n) #runvector
  nlist<-which(cc.lab==1) #indices of cases in cc.lab
  for (i in nlist)
  {
    crit<-0; j<-k
    while (crit<k)
    {
      crit<-sum(cc.lab[order(ipd[i,])[2:(j+1)]]==0)
      j<-j+1
    }
    Tkinv<-Tkinv+(j-k-1) #j-k-1 to exclude the original case point and the \eqn{k^{th}} ctrl reached
    rv[i]<-j-k-1
  }
  list(Tkinv=Tkinv,
       run.vec=rv
  )
} #end for the function
#'

#################################################################

#' @title Expected Value of Cuzick and Edwards \eqn{T_k^{inv}} Test statistic
#'
#' @description
#' This function computes the expected value of Cuzick and Edwards \eqn{T_k^{inv}} test statistic which is based on the 
#' sum of number of cases closer to each case than the \code{k}-th nearest control to the case.
#' 
#' The number of cases are denoted as \eqn{n_1} (denoted as \code{n1} as an argument)
#' and number of controls as \eqn{n_0} for both functions (denoted as \code{n0} as an argument),
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' See the function \code{\link{ceTkinv}} for the details of the \eqn{T_k^{inv}} test.
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}) and references therein.
#' 
#' @param n1,n0 The number of cases and controls
#' @param k Integer specifying the number of the closest controls to subject \eqn{i}.
#' 
#' @return The expected value of Cuzick and Edwards \eqn{T_k^{inv}} test statistic for disease clustering
#'
#' @seealso \code{\link{ceTkinv}}, \code{\link{ceTrun}}, and \code{\link{EV.Trun}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n1<-20
#' n0<-25
#' k<-2 #try also 2, 3
#'
#' EV.Tkinv(n1,n0,k)
#'
#' EV.Tkinv(n1,n0,k=1)
#' EV.Trun(n1,n0)
#' 
#' @export 
EV.Tkinv <- function(n1,n0,k)
{
  k*n1*(n1-1)/(n0+1)
} #end for the function
#'

#################################################################

#' @title Simulated Variance of Cuzick and Edwards \eqn{T_k^{inv}} Test statistic
#'
#' @description
#' This function estimates the variance of Cuzick and Edwards \eqn{T_k^{inv}} test statistic by Monte Carlo simulations
#' under the RL hypothesis.
#' 
#' The exact variance of \eqn{T_k^{inv}} is currently not available and (\insertCite{cuzick:1990;textual}{nnspat}) say
#' that "The permutational variance of \eqn{T_k^{inv}} becomes unwieldy for \eqn{k > 1} and is more easily simulated", hence
#' we estimate the variance of \eqn{T_k^{inv}} by RL of cases and controls to the given point data.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly. The argument \code{Nsim} represents the number of resamplings (without replacement) in the
#' RL scheme, with default being \code{1000}. 
#' 
#' See (\insertCite{cuzick:1990;textual}{nnspat}).
#' 
#' See the function \code{\link{ceTkinv}} for the details of the \eqn{T_k^{inv}} test.
#' 
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point, 
#' @param k Integer specifying the number of the closest controls to subject \eqn{i}.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param Nsim The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
#' variance of \eqn{T_k^{inv}}
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#'  
#' @return The simulation estimated variance of Cuzick and Edwards \eqn{T_k^{inv}} test statistic for disease clustering
#'
#' @seealso \code{\link{ceTkinv}} and \code{\link{EV.Tkinv}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' set.seed(123)
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#' n1<-sum(cls==1)
#' k<-2
#'
#' Nmc<-1000
#' varTkinv.sim(Y,k,cls,Nsim=Nmc)
#'
#' set.seed(1)
#' varTrun.sim(Y,cls,Nsim=Nmc)
#' set.seed(1)
#' varTkinv.sim(Y,k=1,cls,Nsim=Nmc)
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' varTkinv.sim(Y,k,fcls,Nsim=Nmc,case.lab="a")
#' 
#' @export 
varTkinv.sim <- function(dat,k,cc.lab,Nsim=1000,case.lab=NULL)
{
  if (length(table(cc.lab))!=2)
  {stop('cc.lab must have two types/labels for the case-control setting')}
  
  if (!is.null(case.lab))
  {lab2<-rep(0,length(cc.lab))
  lab2[cc.lab==case.lab]<-1
  cc.lab<-lab2}
  
  ce.vec<-vector()
  for (i in 1:Nsim)
  {
    slab<-sample(cc.lab) #random labeling of cases and ctrls
    ce<-ceTkinv(dat,k,slab)
    ce.vec<-c(ce.vec,ce$Tk)
  }
  var(ce.vec)
} #end for the function
#'
 
#################################################################

# funsZTkinv
#'
#' @title Z-Test for Cuzick and Edwards \eqn{T_k^{inv}} statistic
#'
#' @description
#' Two functions: \code{ZTkinv} and \code{ZTkinv.sim}, each of which is an object of class \code{"htest"} performing a
#' \eqn{z}-test for Cuzick and Edwards \eqn{T_k^{inv}} test statistic. See \code{\link{ceTkinv}} for a description of 
#' \eqn{T_k^{inv}} test statistic.
#' 
#' The function \code{ZTkinv} performs a \eqn{Z}-test for \eqn{T_k^{inv}} using asymptotic normality with a simulation estimated
#' variance under RL of cases and controls to the given points.
#' And the function \code{ZTkinv.sim} performs test for\eqn{T_k^{inv}} based on MC simulations under the RL hypothesis.
#'  
#' Asymptotic normality for the \eqn{T_k^{inv}} is not established yet, but this seems likely according to 
#' \insertCite{cuzick:1990;textual}{nnspat}. 
#' If asymptotic normality holds, it seems a larger sample size would be needed before this becomes
#' an effective approximation.
#' Hence the simulation-based test \code{ZTkinv.sim} is recommended for use to be safe. 
#' When \code{ZTkinv} is used, this is also highlighted with the warning "asymptotic normality of \eqn{T_k^{inv}} is not yet established, 
#' so simulation-based test is recommended".
#' 
#' All arguments are common for both functions, except for \dots, Nvar.sim which are used in \code{ZTkinv} only,
#' and \code{Nsim}, which is used in \code{ZTkinv.sim} only.
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly.
#' The argument \code{Nvar.sim} represents the number of resamplings (without replacement) in the
#' RL scheme, with default being \code{1000} for estimating the variance of \eqn{T_k^{inv}} statistic in \code{ZTkinv}.
#' The argument \code{Nsim} represents the number of resamplings (without replacement) in the
#' RL scheme, with default being \code{1000} for estimating the \eqn{T_k^{inv}} values in \code{ZTkinv.sim}.
#' 
#' Both functions might take a very long time when data size is large or \code{Nsim} is large.
#' 
#' See also (\insertCite{cuzick:1990;textual}{nnspat}) and the references therein.
#'  
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point, 
#' used in both functions.
#' @param k Integer specifying the number of the closest controls to subject \eqn{i}, used in both functions.
#' @param cc.lab Case-control labels, 1 for case, 0 for control, used in both functions.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"} or \code{"greater"},
#' used in both functions.
#' @param conf.level Level of the upper and lower confidence limits, default is \code{0.95}, 
#' for Cuzick and Edwards \eqn{T_k^{inv}} statistic. Used in both functions.
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}, used in both functions.
#' @param Nvar.sim The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
#' variance of Tkinv, used in \code{ZTkinv} only.
#' @param Nsim The number of simulations, i.e., the number of resamplings under the RL scheme to estimate the 
#' \eqn{T_k^{inv}} values, used in \code{ZTkinv.sim} only.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function. Used in \code{ZTkinv} only.
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The \eqn{Z} test statistic for the Cuzick and Edwards \eqn{T_k^{inv}} test}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding alternative. In \code{ZTkinv}
#' this is computed using the standard normal distribution, while in \code{ZTkinv.sim}, it is based on which percentile
#' the observed \eqn{T_k^{inv}} value is among the generated \eqn{T_k^{inv}} values.}
#' \item{conf.int}{Confidence interval for the Cuzick and Edwards \eqn{T_k^{inv}} value
#' at the given confidence level \code{conf.level} and depends on the type of \code{alternative}.}
#' \eqn{z}-critical values are used in the construction of the confidence interval in \code{ZTkinv}, 
#' while the percentile values are used in the generated sample of \eqn{T_k^{inv}} values in \code{ZTkinv.sim} 
#' \item{estimate}{Estimate of the parameter, i.e., the Cuzick and Edwards \eqn{T_k^{inv}} value.}
#' \item{null.value}{Hypothesized null value for the Cuzick and Edwards \eqn{T_k^{inv}} value
#' which is \eqn{k n_1 (n_1-1)/(n_0+1)} under RL, where the number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0}.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}}
#' 
#' @seealso \code{\link{ceTkinv}} and \code{\link{EV.Tkinv}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsZTkinv
NULL
#'
#' @rdname funsZTkinv
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' set.seed(123)
#' n<-10 #try also 20, 50, 100
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' k<-2
#'
#' ZTkinv(Y,k,cls)
#' ZTkinv(Y,k,cls+1,case.lab = 2,alt="l")
#' ZTkinv(Y,k,cls,conf=.9,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ZTkinv(Y,k,fcls,case.lab="a")
#' 
#' @export
ZTkinv <- function(dat,k,cc.lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,
                   case.lab=NULL,Nvar.sim=1000,...)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  n<-nrow(dat)
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tki<-ceTkinv(dat,k,cc.lab,case.lab,...)$Tk
  
  varTki<-varTkinv.sim(dat,k,cc.lab,Nsim=Nvar.sim,case.lab) #var of Tkinv
  ETki<-EV.Tkinv(n1,n0,k)
  
  stderr<-sqrt(varTki)
  ts<-(Tki-ETki)/stderr
  names(ts) <-paste(c("Test statistic for Cuzick-Edwards Tkinv test for k =", k," with variance estimated via simulation, Z"), collapse=" ")
  
  estimate<-Tki
  names(estimate) <-c("Cuzick-Edwards Tkinv")
  
  method <- "Z-Test for Cuzick-Edwards Tkinv Test for Clustering"
  
  null.val<- ETki
  names(null.val) <-"(expected) value of Cuzick-Edwards Tkinv under the null hypothesis"
  
  alt<- switch(alternative,
         less = { 
           pval <-pnorm(ts)
           cint <-estimate+c(-Inf, qnorm(conf.level))*stderr
         },
         greater = { 
           pval <-pnorm(ts, lower.tail = FALSE)
           cint <-estimate+c(-qnorm(conf.level),Inf)*stderr
         },
         two.sided = { 
           pval <-2 * pnorm(-abs(ts))
           alpha <-1 - conf.level
           cint <-qnorm(1 - alpha/2)
           cint <-estimate+c(-cint, cint)*stderr
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(dat))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  warning("asymptotic normality of Tkinv is not yet established, so simulation-based test is recommended")
  return(rval)
} #end for the function
#'
#' @rdname funsZTkinv
#'
#' @examples
#' n<-10 #try also 20, 50, 100
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)  #or try cls<-rep(0:1,c(10,10))
#' k<-2 # try also 3,5
#'
#' ZTkinv.sim(Y,k,cls)
#' ZTkinv.sim(Y,k,cls+1,case.lab = 2,alt="l")
#' ZTkinv.sim(Y,k,cls,conf=.9,alt="g")
#'
#' #cls as a factor
#' na<-floor(n/2); nb<-n-na
#' fcls<-rep(c("a","b"),c(na,nb))
#' ZTkinv.sim(Y,k,fcls,case.lab="a")
#'
#' #with k=1
#' ZTkinv.sim(Y,k=1,cls)
#' ZTrun(Y,cls)
#' 
#' @export
ZTkinv.sim <- function(dat,k,cc.lab,alternative=c("two.sided", "less", "greater"),conf.level = 0.95,
                       case.lab=NULL,Nsim=1000)
{
  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  
  n<-nrow(dat)
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tk0<-ceTkinv(dat,k,cc.lab,case.lab)$Tk
  
  if (!is.null(case.lab))
  {lab2<-rep(0,length(cc.lab))
  lab2[cc.lab==case.lab]<-1
  cc.lab<-lab2}
  
  ce.vec<-vector()
  for (i in 1:Nsim)
  {
    slab<-sample(cc.lab) #random labeling of cases and ctrls
    ce<-ceTkinv(dat,k,slab)
    ce.vec<-c(ce.vec,ce$Tk)
  }
  
  ETk0<-EV.Tkinv(n1,n0,k)
  ts<-Tk0
  names(ts) <-paste(c("Cuzick-Edwards Tkinv value with k =", k," for the current data"), collapse=" ")
  
  estimate<-Tk0
  names(estimate) <-c("Cuzick-Edwards Tkinv")
  
  method <- "Monte Carlo Test for Cuzick-Edwards Tkinv for Clustering"
  
  null.val<- ETk0
  names(null.val) <-"(expected) value of Cuzick-Edwards Tkinv under the null hypothesis of RL"
  
  alpha <-1 - conf.level
  quants<-quantile(ce.vec,probs=c(alpha/2,alpha,1-alpha,1-alpha/2))
  pval.ls <-(sum(ce.vec<=Tk0)+1)/(Nsim+1); pval.rs <-(sum(ce.vec>=Tk0)+1)/(Nsim+1)
  alt<- switch(alternative,
         less = { 
           pval<-pval.ls
           cint <-c(-Inf, quants[3])
         },
         greater = { 
           pval<-pval.rs
           cint <-c(quants[2],Inf)
         },
         two.sided = { 
           pval <-2 * min(pval.ls,pval.rs)
           alpha <-1 - conf.level
           cint <-c(quants[1],quants[4])
         }
  )
  
  if (is.null(alt)) stop("Alternative must be one of less, greater, or two.sided in quotes")
  
  attr(cint, "conf.level") <-conf.level 
  
  dname <-deparse(substitute(dat))
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.val,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  return(rval)
} #end for the function
#'

###################################################################################
###FUNCTIONS for TANGO's CORRECTION to CUZICK and EDWARDS k-NN TESTS###
###################################################################################

# funsW345values
#'
#' @title \eqn{W_k} values for Tango's \eqn{T} test statistic
#'
#' @description
#' Three functions: \code{W3val}, \code{W4val} and \code{W5val}, each of which is needed to compute \eqn{E[T^3]}
#' (i.e., for the skewness of \eqn{T})
#' where \eqn{T=T(\theta)} which is defined in Equation (2) of \insertCite{tango:2007;textual}{nnspat} as follows:
#' Let \eqn{(z_1,\ldots,z_n )}, \eqn{n = n_0 + n_1}, denote the locations of the points in the combined sample 
#' when the indices have been randomly permuted so that the \eqn{z_i} contain no information about group membership.
#' \deqn{T(\theta)=\sum_{i=1}^{n}\sum_{j=1}^{n}\delta_i \delta_j a_{ij}(\theta)=
#' \boldsymbol \delta^t \boldmath A(\theta)) \boldsymbol \delta} where \eqn{\delta_i=1} if \eqn{z_i} is a case,
#' and 0 if \eqn{z_i} is a control,  \eqn{\boldmath A(\theta) = (a_{ij} (\theta))} could be any matrix of a measure of
#' the closeness between two points \eqn{i} and \eqn{j} with \eqn{a_{ii} = 0} for all \eqn{i = 1,\ldots,n}, and \eqn{\boldsymbol \theta = 
#' (\theta_1,\ldots,\theta_p)^t} denotes the unknown parameter vector related to cluster size and 
#' \eqn{\boldsymbol \delta = (\delta_1,\ldots,\delta_n)^t}. 
#' Here the number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0}  to match the case-control class
#' labeling, which is just the reverse of the labeling in \insertCite{tango:2007;textual}{nnspat}.
#' 
#' If \eqn{\theta=k} in the nearest neighbors model with \eqn{a_{ij}(k) = 1} if \eqn{z_j} is among the \eqn{k}NNs of \eqn{z_i} and 0 
#' otherwise, then the test statistic \eqn{T(\theta) = T_k} is the Cuzick and Edwards \eqn{k}NN test statistic, \eqn{T_k}
#' \insertCite{cuzick:1990;textual}{nnspat}, see also \code{\link{ceTk}}.
#' 
#' \eqn{W_k} values are used for Tango's correction to Cuzick and Edwards \eqn{k}NN test statistic, \eqn{T_k} and
#' \eqn{W_k} here corresponds to \eqn{W_{k-1}} in \insertCite{tango:2007;textual}{nnspat}
#' (defined for consistency with \eqn{p_k}'s and \eqn{alpha_r} having \eqn{r} distinct elements).
#' 
#' The argument of the function is the \eqn{A_{ij}} matrix, \code{a}, which is the output of the function \code{\link{aij.mat}}.
#' However, inside the function we symmetrize the matrix \code{a} as \code{b <- (a+a^t)/2}, to facilitate the formulation.
#' 
#' @param a \eqn{A_{ij}} matrix which is the output of the function \code{\link{aij.mat}}.
#' 
#' @return Each function \code{Wkval} returns the \eqn{W_k} value for \eqn{k=3,4,5}.
#' 
#' @seealso \code{\link{ceTk}}, \code{\link{EV.Tk}}, \code{\link{varTk}}, \code{\link{Xsq.ceTk}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @name funsW345values
NULL
#'
#' @rdname funsW345values
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' k<-sample(1:5,1) # try also 3, 5, sample(1:5,1)
#' k
#' a<-aij.mat(Y,k)
#' W3val(a)
#' W4val(a)
#' W5val(a)
#'
#' a<-aij.mat(Y,k,method="max")
#' W3val(a)
#' W4val(a)
#' W5val(a)
#'
#' @export 
W3val <- function(a)
{ b<-1/2*(a+t(a))
n<-nrow(b)
W3<-0
for (i in 1:n)
{bi<-sum(b[i,])
for (j in 1:n)
{
  Sij<-sum(b[i,]*b[j,]) #short for the below for loop
  #S<-0 #for (m in 1:n)
  #{S<-S+b[i,m]*b[j,m]}
  W3<-W3+24*b[i,j]^2*bi-24*b[i,j]^3+8*b[i,j]*Sij
}
}
W3
} #end for the function
#'
#' @rdname funsW345values
#'
#' @export
W4val <- function(a)
{
  b<-1/2*(a+t(a))
  n<-nrow(b)
  S1<-sum(b)
  W4<-0
  for (i in 1:n)
  {bi<-sum(b[i,])
  for (j in 1:n)
  { bj<-sum(b[j,])
  Sij<-sum(b[i,]*b[j,]) #short for the below for loop
  # S<-0      for (m in 1:n)
  #{S<-S+b[i,m]*b[j,m]}
  W4<-W4+6*S1*b[i,j]^2+52*b[i,j]^3-96*b[i,j]^2*bi+
    8*b[i,j]*bi^2+24*b[i,j]*bi*bj-24*b[i,j]*Sij
  }
  }
  W4
} #end for the function
#'
#' @rdname funsW345values
#'
#' @export 
W5val <- function(a)
{ b<-1/2*(a+t(a))
n<-nrow(b)
S1<-sum(b)
W5<-0
for (i in 1:n)
{bi<-sum(b[i,])
for (j in (1:n)[-i])
{bj<-sum(b[j,])
Sij<-sum(b[i,]*b[j,]) #short for the below for loop
#S<-0     for (m in 1:n)
#{S<-S+b[i,m]*b[j,m]}
W5<-W5-S1*b[i,j]^2-4*b[i,j]^3+10*b[i,j]^2*bi-b[i,j]*bi^2-
  4*b[i,j]*bi*bj+bi^2*bj+2*b[i,j]*Sij
}
}
12*W5
} #end for the function
#'

#################################################################

#' @title Skewness of Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description
#' This function estimates the skewness of Cuzick and Edwards \eqn{T_k} test statistic under the RL hypothesis.
#' Skewness of a random variable \eqn{T} is defined as \eqn{E(T-\mu)^3/(E(T-\mu)^2)^{1.5}} where \eqn{\mu=E T}.
#' 
#' Skewness is used for Tango's correction to Cuzick and Edwards \code{k}NN test statistic, \eqn{T_k}.
#' Tango's correction is a chi-square approximation, and its degrees of freedom is estimated using the skewness
#' estimate (see page 121 of \insertCite{tango:2007;textual}{nnspat}).
#' 
#' The argument, \eqn{n_1}, is the number of cases (denoted as \code{n1} as an argument) 
#' and \code{k} is the number of NNs considered in \eqn{T_k} test statistic.
#' The argument of the function is the \eqn{A_{ij}} matrix, \code{a}, which is the output of the function \code{\link{aij.mat}}.
#' However, inside the function we symmetrize the matrix \code{a} as \code{b <- (a+a^t)/2}, to facilitate the formulation.
#' 
#' The number of cases are denoted as \eqn{n_1} and number of controls as \eqn{n_0} in this function
#' to match the case-control class labeling,
#' which is just the reverse of the labeling in \insertCite{cuzick:1990;textual}{nnspat}.
#' 
#' @param n1 Number of cases
#' @param k Integer specifying the number of NNs (of subject \eqn{i})
#' @param a \eqn{A_{ij}} matrix which is the output of the function \code{\link{aij.mat}}.
#'  
#' @return The skewness of Cuzick and Edwards \eqn{T_k} test statistic for disease clustering
#'
#' @seealso \code{\link{ceTk}}, \code{\link{EV.Tk}}, and \code{\link{varTk}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#' n1<-sum(cls==1)
#'
#' k<-sample(1:5,1) # try also 3, 5, sample(1:5,1)
#' k
#' a<-aij.mat(Y,k)
#'
#' SkewTk(n1,k,a)
#' 
#' @export 
SkewTk <- function(n1,k,a)
{ b<-1/2*(a+t(a))
n<-nrow(b)
n0<-n-n1
S1<-sum(b)
S2<-sum(b^2)
S3<-sum(b^3)

G1<-G2<-G3<-G4<-G5<-0
for (i in 1:n)
{bi.<-sum(b[i,])
for (j in (1:n)[-i])
{ bj.<-sum(b[j,])
Sij<-sum(b[i,]*b[j,]) 
G1<-G1+b[i,j]^2*bi.
G2<-G2+b[i,j]*Sij
G3<-G3+b[i,j]*bi.*bj.
G4<-G4+b[i,j]*(bi.)^2
G5<-G5+(bi.)^2*bj.
}
}

W2<-4*S3
W3<-24*(G1-S3)+8*G2
W4<-6*S1*S2+52*S3-96*G1+8*G4+24*G3-24*G2
W5<-12*(-S1*S2-4*S3+10*G1-G4-4*G3+G5+2*G2)
W6<-S1^3-(W2+W3+W4+W5)

ETk<-EV.Tk(k,n1,n0)
var.Tk<-varTkaij(n1,k,a)

p2<-pk(n,n1,2); p3<-pk(n,n1,3); p4<-pk(n,n1,4)
p5<-pk(n,n1,5); p6<-pk(n,n1,6)

skewness<-((p2*W2+p3*W3+p4*W4+p5*W5+p6*W6)-3*ETk*var.Tk-ETk^3)/((var.Tk)^(1.5))
skewness
} #end for the function
#'

#################################################################

#' @title Chi-square Approximation to Cuzick and Edwards \eqn{T_k} Test statistic
#'
#' @description  
#' An object of class \code{"Chisqtest"} performing  a chi-square approximation for Cuzick and Edwards \eqn{T_k} test statistic
#' based on the number of cases within \code{k}NNs of the cases in the data.
#' 
#' This approximation is suggested by \insertCite{tango:2007;textual}{nnspat} since \eqn{T_k} statistic had high 
#' skewness rendering the normal approximation less efficient. The chi-square approximation is as follows:
#' \eqn{\frac{T_k- ET_k}{\sqrt{Var T_k}} \approx \frac{\chi^2_\nu-\nu}{\sqrt{2 \nu}}} where \eqn{\chi^2_\nu} is a chi-square
#' random variable with \eqn{\nu} df, and \eqn{\nu=8/skewnees(T_k)} (see \code{\link{SkewTk}} for the skewness).
#' 
#' The argument \code{cc.lab} is case-control label, 1 for case, 0 for control, if the argument \code{case.lab} is \code{NULL}, 
#' then \code{cc.lab} should be provided in this fashion, if \code{case.lab} is provided, the labels are converted to 0's 
#' and 1's accordingly.
#' 
#' The logical argument \code{nonzero.mat} (default=\code{FALSE}) is for using the \eqn{A} matrix if \code{FALSE} or just the matrix of nonzero
#' locations in the \eqn{A} matrix (if \code{TRUE}).
#' 
#' The logical argument \code{asy.var} (default=\code{FALSE}) is for using the asymptotic variance or the exact (i.e. finite
#' sample) variance for the variance of \eqn{T_k} in its standardization.
#' If \code{asy.var=TRUE}, the asymptotic variance is used for \eqn{Var[T_k]} (see \code{asyvarTk}), otherwise the exact
#' variance (see \code{varTk}) is used.
#' 
#' See also (\insertCite{tango:2007;textual}{nnspat}) and the references therein.
#'
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param cc.lab Case-control labels, 1 for case, 0 for control
#' @param k Integer specifying the number of NNs (of subject \eqn{i}).
#' @param case.lab The label used for cases in the \code{cc.lab} (if \code{cc.lab} is not provided then the labels are converted
#' such that cases are 1 and controls are 0), default is \code{NULL}.
#' @param nonzero.mat A logical argument (default is \code{TRUE}) to determine whether the \eqn{A} matrix or the matrix of
#' nonzero locations of the \eqn{A} matrix will be used in the computations.
#' If \code{TRUE} the nonzero location matrix is used, otherwise the \eqn{A} matrix itself is used.
#' @param asy.var A logical argument (default is \code{FALSE}) to determine whether the asymptotic variance or 
#' the exact (i.e. finite sample) variance for the variance of \eqn{T_k} in its standardization. 
#' If \code{TRUE}, the asymptotic variance is used for \eqn{Var[T_k]}, otherwise the exact variance is used.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return A \code{list} with the elements
#' \item{statistic}{The chi-squared test statistic for Tango's chi-square approximation to Cuzick & Edwards' \eqn{T_k}
#' test for disease clustering.}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test}
#' \item{df}{Degrees of freedom for the chi-squared test, which is \eqn{8/}skewness where skewness is the output of
#' \code{\link{SkewTk}} function.}
#' \item{estimate}{Estimates, i.e., the observed \eqn{T_k} value.}
#' \item{est.name,est.name2}{Names of the estimates, they are almost identical for this function.}
#' \item{null.value}{Hypothesized null value for Cuzick & Edwards' \eqn{T_k}, which is \eqn{ET_k}.}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set, \code{dat}}
#'  
#' @seealso \code{\link{ceTk}}, \code{\link{ZceTk}} and \code{\link{SkewTk}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' set.seed(123)
#' n<-20
#' Y<-matrix(runif(3*n),ncol=3)
#' cls<-sample(0:1,n,replace = TRUE)
#'
#' k<-sample(1:5,1) # try also 1, 3, 5,
#' k
#'
#' Xsq.ceTk(Y,cls,k)
#' Xsq.ceTk(Y,cls,k,nonzero.mat=FALSE)
#' Xsq.ceTk(Y,cls+1,k,case.lab = 2)
#' Xsq.ceTk(Y,cls,k,method="max")
#'
#' Xsq.ceTk(Y,cls,k,asy.var=TRUE)
#'
#' @export
Xsq.ceTk <- function(dat,cc.lab,k,case.lab=NULL,nonzero.mat=TRUE,asy.var=FALSE,...)
{
  n<-nrow(dat)
  ifelse(is.null(case.lab),
         n1<-sum(cc.lab==1), n1<-sum(cc.lab==case.lab)) #n1 is number of cases
  n0<-n-n1 #n0 is number of controls
  
  Tk<-ceTk(dat,cc.lab,k,case.lab,...)
  ETk<-EV.Tk(k,n1,n0)
  ifelse(asy.var==TRUE,var.Tk<-asyvarTk(dat,n1,k,nonzero.mat,...)$asy.var,
         var.Tk<-varTk(dat,n1,k,nonzero.mat,...)$var)
  stderr <-sqrt(var.Tk)
  zTk<-(Tk-ETk)/stderr #TS.CE
  
  a<-aij.mat(dat,k,...)
  skw<-SkewTk(n1,k,a)
  nu<-8/(skw^2)
  ts<-Xsq<-nu+sqrt(2*nu)*zTk
  
  if (is.na(ts))
  {stop("The test statistic is NaN, so chisquare approximation for Cuzick-Edwards Tk is not defined")}
  
  pval<-pchisq(ts,nu,lower.tail = F)
  
  method <-"Chi-square Approximation for Cuzick-Edwards kNN Test"
  
  dname <-deparse(substitute(dat))
  
  estimate<-Tk
  estimate.name <-paste("Cuzick-Edwards Tk for k=",k,sep="")
  estimate.name2 <-"Cuzick-Edwards Tk"
  
  null.val<- ETk
  names(null.val) <-"(expected) value of Cuzick-Edwards Tk under the null hypothesis"
  
  rval <-list(
    statistic=ts,
    p.value=pval,
    df=nu,
    estimate = estimate,
    est.name = estimate.name,
    est.name2 = estimate.name2,
    null.value = ETk,
    method = method,
    data.name = dname
  )
  
  attr(rval, "class") <-"Chisqtest"
  return(rval)
} #end for the function
#'

#################################################################

#' @title Closeness or Proximity Matrix for Tango's Spatial Clustering Tests
#'
#' @description 
#' This function computes the \eqn{A=a_{ij}(\theta)} matrix useful in calculations for Tango's test \eqn{T(\theta)} 
#' for spatial (disease) clustering (see Eqn (2) of \insertCite{tango:2007;textual}{nnspat}.
#' Here, \eqn{A=a_{ij}(\theta)} is any matrix of a measure of the closeness between two points \eqn{i} and \eqn{j} with \eqn{aii = 0} for all
#' \eqn{i = 1, \ldots,n}, and \eqn{\theta = (\theta_1,\ldots,\theta_p)^t} denotes the unknown parameter vector related 
#' to cluster size and \eqn{\delta = (\delta_1,\ldots,\delta_n)^t}, where \eqn{\delta_i=1} if \eqn{z_i} is a case and 0 
#' otherwise.
#' The test is then
#' \deqn{T(\theta)=\sum_{i=1}^n\sum_{j=1}^n\delta_i \delta_j a_{ij}(\theta)=\delta^t A(\theta) \delta}
#' where \eqn{A=a_{ij}(\theta)}.
#' 
#' \eqn{T(\theta)} becomes Cuzick and Edwards \eqn{T_k} tests statistic (\insertCite{cuzick:1990;textual}{nnspat}),
#' if \eqn{a_{ij}=1} if \eqn{z_j} is among the \code{k}NNs of \eqn{z_i} and 0 otherwise.
#' In this case \eqn{\theta=k} and \code{aij.theta} becomes \code{aij.mat} (more specifically,
#' \code{aij.mat(dat,k)} and \code{aij.theta(dat,k,model="NN")}.
#' 
#' In Tango's exponential clinal model (\insertCite{tango:2000;textual}{nnspat}),
#' \eqn{a_{ij}=\exp\left(-4 \left(\frac{d_{ij}}{\theta}\right)^2\right)} if \eqn{i \ne j}  and 0 otherwise,
#' where \eqn{\theta} is a predetermined scale of cluster such that any pair of cases far apart beyond the distance 
#' \eqn{\theta} cannot be considered as a cluster and \eqn{d_{ij}} denote the Euclidean distance between 
#' two points \eqn{i} and \eqn{j}. 
#' 
#' In the exponential model (\insertCite{tango:2007;textual}{nnspat}),
#' \eqn{a_{ij}=\exp\left(-\frac{d_{ij}}{\theta}\right)} if \eqn{i \ne j}  and 0 otherwise,
#' where \eqn{\theta} and \eqn{d_{ij}} are as above.
#' 
#' In the hot-spot model (\insertCite{tango:2007;textual}{nnspat}),
#' \eqn{a_{ij}=1} if \eqn{d_{ij} \le \theta} and \eqn{i \ne j}  and 0 otherwise,
#' where \eqn{\theta} and \eqn{d_{ij}} are as above.
#' 
#' The argument \code{model} has four options, \code{NN}, \code{exp.clinal}, \code{exponential}, and 
#' \code{hot.spot}, with \code{exp.clinal} being the default.
#' And the \code{theta} argument specifies the scale of clustering or the clustering parameter in the particular
#' spatial disease clustering model.
#' 
#' See also (\insertCite{tango:2007;textual}{nnspat}) and the references therein.
#'
#' @param dat The data set in one or higher dimensions, each row corresponds to a data point.
#' @param theta A predetermined cluster scale so that any pair of cases farther apart then the distance 
#' \eqn{\theta} is unlikely to be cluster.
#' @param model Type of Tango's spatial clustering model with four options: 
#' \code{NN}, \code{exp.clinal} (default), \code{exponential}, and \code{hot.spot}.
#' @param \dots are for further arguments, such as \code{method} and \code{p}, passed to the \code{\link[stats]{dist}} function.
#' 
#' @return The \eqn{A=a_{ij}(\theta)} matrix useful in calculations for Tango's test \eqn{T(\theta)}.
#'  
#' @seealso \code{\link{aij.mat}}, \code{\link{aij.nonzero}} and \code{\link{ceTk}}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-20  #or try sample(1:20,1)
#' Y<-matrix(runif(3*n),ncol=3)
#' k<-3#1 #try also 2,3
#'
#' #aij for CE's Tk
#' Aij<-aij.theta(Y,k,model = "NN")
#' Aij2<-aij.mat(Y,k)
#' sum(abs(Aij-Aij2)) #check equivalence of aij.theta and aij.mat with model="NN"
#'
#' Aij<-aij.theta(Y,k,method="max")
#' Aij2<-aij.mat(Y,k)
#' range(Aij-Aij2)
#'
#' theta=.2
#' aij.theta(Y,theta,model = "exp.clinal")
#' aij.theta(Y,theta,model = "exponential")
#' aij.theta(Y,theta,model = "hot.spot")
#' 
#' @export
aij.theta <- function(dat,theta,model="exp.clinal",...)
{ 
  ipd<-ipd.mat(dat,...)
  n<-nrow(dat)
  a<-matrix(0,n,n)
  
  if (n<=1)
  { return(a)}
 a<- switch(model,
         NN = { a<-aij.mat(dat,theta,...) },
         exp.clinal = {   a<-matrix(0,n,n)
         for (i in 1:(n-1))
         {
           for (j in (i+1):n)
           {
             a[i,j]<-exp(-4*(ipd[i,j]/theta)^2)
             a[j,i]<-a[i,j]
           }
         }
         a },
         exponential = {   a<-matrix(0,n,n)
         for (i in 1:(n-1))
         {
           for (j in (i+1):n)
           {
             a[i,j]<-exp(-(ipd[i,j]/theta))
             a[j,i]<-a[i,j]
           }
         }
         a },
         hot.spot = {   a<-matrix(0,n,n)
         for (i in 1:(n-1))
         {
           for (j in (i+1):n)
           {
             a[i,j]<-sum(ipd[i,j] <= theta)
             a[j,i]<-a[i,j]
           }
         }
         a}
  )
 
 if (is.null(a)) stop("Enter numbers 1-4 or NN, exp.clinal, exponential, hotspot in quotes for model")
  a
} #end for the function
#'