# AuxRFuncs4NNCTClasses.R
###############################################
#Auxiliary functions for class cellhtest
###############################################
#'
#' @title Print a summary of a \code{cellhtest} object
#'
#' @description 
#' Printing objects of class "\code{cellhtest}" by simple \code{\link[base]{print}} methods.
#'
#' @param x	 object of class "\code{summary.cellhtest}"
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link[base]{strwrap}} for displaying the method component 
#' of the \code{classhtest} object.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @export
print.cellhtest <- function(x, digits = getOption("digits"), prefix = "\t",...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  
  if (!inherits(x, "cellhtest"))
    stop("object must be of class \"cellhtest\"")
  
  if (is.null(x$data.name)) 
  {cat("contingency table:",x$ct.name,"\n")
  } else
  {cat("data:",x$data.name,"\n")}
  cat(x$stat.names,"=\n",sep = "")
  print(x$statistic)
  
  cat("\n")
  
  cat("p-values=\n",sep = "")
  print(x$p.value)
  
  cat("\n")
  
  #this is for the alternative part
  if(!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if(!is.null(x$null.value)) {
      
      alt.char <-
        switch(x$alternative,
               two.sided = "not equal to",
               less = "less than",
               greater = "greater than")
      
      if (all(x$null.value==0))
      {
        cat("true (expected value of each of the) ", x$null.name, " is ", alt.char," ", 0, "\n",sep = "")
      } else
      {
        cat("true (expected value of the) ", x$null.name, " is ", alt.char, "\n",sep = "")
        print(x$null.value, digits=digits,...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  cat("\n")
  #for the confidence intervals (LCL and UCL as matrices)
  if(!is.null(x$LCL)) {
    cat(100 * x$cnf.lvl, "% lower confidence limits of ",x$est.name2,":\n",sep = "")
    print(x$LCL)
    cat("\n")
  }
  
  if(!is.null(x$UCL)) {
    cat(100 * x$cnf.lvl, "% upper confidence limits of ",x$est.name2,":\n",sep = "")
    print(x$UCL)
    cat("\n")
  }
  
  #for the confidence intervals
  if(!is.null(x$conf.int)) {
    cat(100 * x$cnf.lvl, "% confidence intervals for ",x$est.name2,":\n",sep = "")
    print(x$conf.int)
    cat("\n")
  }
  
  if(!is.null(x$estimate)) {
    cat("sample estimates of ",x$est.name2,":\n",sep = "")
    print(x$estimate,...)
  }
} #end of the function
#'
###############################################
#Auxiliary functions for class classhtest
###############################################
#'
#' @title Print a summary of a \code{classhtest} object
#'
#' @description 
#' Printing objects of class "\code{classhtest}" by simple \code{\link[base]{print}} methods.
#'
#' @param x	 object of class "\code{summary.classhtest}"
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link[base]{strwrap}} for displaying the method component 
#' of the \code{classhtest} object.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @export
print.classhtest <- function(x, digits = getOption("digits"), prefix = "\t",...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  
  if (!inherits(x, "classhtest"))
    stop("object must be of class \"classhtest\"")
  
  est<-x$estimate
  clnames<-rownames(est) 
  k<-nrow(est)
  d.f<-x$df
  
  stat<-x$statistic
  p.val<-x$p.value
  
  if (is.null(x$data.name)) 
  {cat("contingency table:",x$ct.name,"\n")
  } else
  {cat("data:",x$data.name,"\n")}
  cat(x$stat.names," =\n",sep = "")
  for (i in 1:k)
    cat(x$type,"class",clnames[i],": X-squared =",stat[i],"  p-value =",p.val[i],"  df =",d.f,"\n")
  
  cat("\n")
  
  #this is for the alternative part
  cat("alternative hypothesis: true expected values of NNCT entries for \n")
  for (i in 1:k)
    cat(x$type,"class",clnames[i],"are different from",x$null.value[i,], "\n")
  
  cat("\n")
  
  if(!is.null(x$estimate)) {
    cat("sample estimates of NNCT entries for \n")
    for (i in 1:k)
      cat(x$type,"class",clnames[i],"are",x$estimate[i,], "\n")
  }
} #end of the function
#'
#'
###############################################
#Auxiliary functions for class Chisqtest
###############################################
#'
#' @title Print a summary of a \code{Chisqtest} object
#'
#' @description 
#' Printing objects of class "\code{Chisqtest}" by simple \code{\link[base]{print}} methods.
#'
#' @param x	 object of class "\code{summary.Chisqtest}"
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link[base]{strwrap}} for displaying the method component 
#' of the \code{classhtest} object.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @export
print.Chisqtest <- function(x, digits = getOption("digits"), prefix = "\t",...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  
  if (!inherits(x, "Chisqtest"))
    stop("object must be of class \"Chisqtest\"")
  
  est<-x$estimate
  clnames<-rownames(est) 
  k<-nrow(est)
  d.f<-x$df
  
  stat<-x$statistic
  p.val<-x$p.value
  
  if (is.null(x$data.name)) 
  {cat("contingency table:",x$ct.name,"\n")
  } else
  {cat("data:",x$data.name,"\n")}
  
  cat("Chi-squared =",stat,"  p-value =",p.val,"  df =",d.f,"\n")
  
  cat("\n")
  
  #this is for the alternative part
  if (all(x$null.value==0))
  {
    cat("alternative hypothesis: true expected value of each of ",x$est.name," is different from 0 \n",sep = "")    
  } else
  { 
    if (length(x$estimate)>1)
    {cat("alternative hypothesis: true expected values of ",x$est.name," are different from \n",sep = "")}
    else {cat("alternative hypothesis: true expected value of ",x$est.name," are different from \n",sep = "")}
    print(x$null.value, digits=digits,...)
  }
  
  cat("\n")
  
  if(!is.null(x$estimate)) {
    if (length(x$estimate)>1)
    {cat("sample estimates of ",x$est.name2,":\n",sep = "")}
    else {cat("sample estimate of ",x$est.name2,":\n",sep = "")}
    print(est,...)
  }
} #end of the function
#'
#'
###############################################
#Auxiliary functions for class refhtest
###############################################
#'
#' @title Print a summary of a \code{refhtest} object
#'
#' @description 
#' Printing objects of class "\code{refhtest}" by simple \code{\link[base]{print}} methods.
#'
#' @param x	 object of class "\code{summary.refhtest}"
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link[base]{strwrap}} for displaying the method component 
#' of the \code{classhtest} object.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @export
print.refhtest <- function(x, digits = getOption("digits"), prefix = "\t",...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  
  if (!inherits(x, "refhtest"))
    stop("object must be of class \"refhtest\"")
  
  est<-x$estimate
  clnames<-rownames(est) 
  k<-nrow(est)
  d.f<-x$df
  
  stat<-x$statistic
  p.val<-x$p.value
  
  if (is.null(x$data.name)) 
  {cat("contingency table:",x$ct.name,"\n")
  } else
  {cat("data:",x$data.name,"\n")}
  
  cat(x$stat.names," \n",sep = "")
  cat("self-reflexivity: Z = ",stat[1],"  p-value = ",p.val[1],"\n")
  cat("mixed-non-reflexivity: Z = ",stat[2],"  p-value = ",p.val[2],"\n")
  
  cat("\n")
  
  #this is for the alternative part
  if(!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if(!is.null(x$null.value)) {
      
      alt.char <-
        switch(x$alternative,
               two.sided = "not equal to ",
               less = "less than ",
               greater = "greater than ")
      
      cat("true (expected value of the) ", x$null.name, " are ", alt.char, "\n",sep = "")
      print(x$null.value, digits=digits,...)
    }
    else cat(x$alternative, "\n", sep = "")
  }
  cat("\n")
  #for the confidence intervals
  if(!is.null(x$conf.int)) {
    cat(100 * x$cnf.lvl, " % confidence intervals for diagonal rct entries for \n",sep = "")
    cat("self-reflexivity =",x$conf.int[1,], "\n",sep = " ")
    cat("mixed-non-reflexivity =",x$conf.int[2,], "\n",sep = " ")
  }
  cat("\n")
  
  if(!is.null(x$estimate)) {
    cat("sample estimates of diagonal rct entries are \n")
    print(x$estimate, digits=digits,...)
  }
} #end of the function
#'
#'
###############################################
#Auxiliary functions for class SpatPatterns
###############################################
#'
#' @title Print a \code{SpatPatterns} object
#'
#' @description Prints the \code{call} of the object of class '\code{SpatPatterns}'
#' and also the \code{type} (or description) of the pattern).
#'
#' @param x A \code{SpatPatterns} object.
#' @param \dots Additional arguments for the S3 method '\code{print}'.
#'
#' @return
#' The \code{call} of the object of class '\code{SpatPatterns}'
#' and also the \code{type} (or description) of the pattern).
#'
#' @seealso \code{\link{summary.SpatPatterns}}, \code{\link{print.summary.SpatPatterns}}, and \code{\link{plot.SpatPatterns}}
#'
#' @examples
#' #TBF (to be filled)
#'
#' @export
print.SpatPatterns <- function(x,...)
{
  if (!inherits(x, "SpatPatterns"))
    stop("x must be of class \"SpatPatterns\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
} #end of the function
#'
########################
#'
#' @title Return a summary of a \code{SpatPatterns} object
#'
#' @description Returns the below information about the \code{object}:
#'
#' \code{call} of the function defining the \code{object}, the \code{type} of the pattern, \code{parameters} of the pattern,
#' study window, some sample points from the generated pattern, reference points (if any for the bivariate pattern),
#' and number of points for each class
#'
#' @param object Object of class \code{SpatPatterns}.
#' @param \dots Additional parameters for \code{summary}.
#'
#' @return
#' The \code{call} of the object of class '\code{SpatPatterns}', the \code{type} of the pattern, \code{parameters} of the pattern,
#' study window, some sample points from the generated pattern, reference points (if any for the bivariate pattern),
#' and number of points for each class
#'
# #' @seealso \code{\link{print.SpatPatterns}}, \code{\link{print.summary.SpatPatterns}}, and \code{\link{plot.SpatPatterns}}
#'
#' @examples
#' #TBF
#' @export
summary.SpatPatterns <- function(object,...)
{
  if (!inherits(object, "SpatPatterns"))
    stop("object must be of class \"SpatPatterns\"")
  
  dtlab <-object$lab
  
  typ <- object$type
  ref.pts<-object$ref.points #reference points
  gen.pts<-object$gen.points #generated points
  dat.pts<-object$dat.points #data points
  
  switch(object$pat.type, 
         "1c" = {
           ifelse(is.null(gen.pts),xv<-dat.pts,xv<-gen.pts)
           yv<-NULL
         }, 
         "2c" ={
           dtset<- gen.pts
           xv<-as.matrix(dtset[dtlab==1,])
           yv<-as.matrix(dtset[dtlab==2,])
         }, 
         cc = {
           dtset<- dat.pts
           xv<-as.matrix(dtset[dtlab==1,])
           yv<-as.matrix(dtset[dtlab==0,])
         },
         ref.gen={
           xv<-gen.pts
           yv<-ref.pts
         }
  )
  
  nx<-min(5,nrow(xv))
  ny<-min(5,nrow(yv))
  Npts<-object$num.points
  
  if (nx==0) {xm<-NULL} else {xm<-xv[1:nx,]}
  if (ny==0) {ym<-NULL} else {ym<-yv[1:ny,]}
  
  res <- list(pat.type=object$pat.type,
              pat.desc=object$desc.pat,
              call=object$call,
              xmat=xm,
              ymat=ym,
              param=object$parameters,
              type=typ,
              num.pts=Npts,
              Xlim=object$xlimit,
              Ylim=object$ylimit 
  )
  
  class(res) <- "summary.SpatPatterns"
  res
} #end of the function
#'
########################
#'
#' @title Print a summary of a \code{SpatPatterns} object
#'
#' @description Prints some information about the \code{object}.
#'
#' @param x	 object of class "\code{summary.SpatPatterns}", generated by \code{summary.SpatPatterns}.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @seealso \code{\link{print.SpatPatterns}}, \code{\link{summary.SpatPatterns}}, and \code{\link{plot.SpatPatterns}}
#'
#' @export
print.summary.SpatPatterns <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the Pattern:\n")
  print(x$type)
  
  cat("\nParameters of the Pattern:\n")
  print(x$param)
  
  cat("\nStudy Window:\n")
  cat("range in x-coordinate =", x$Xlim,"\n")
  cat("range in y-coordinate =", x$Ylim,"\n")
  
  switch(x$pat.type, 
         "1c" = {cat("\n Generated Points from",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$xmat)
           cat("\nNumber of Points:\n ")
           print(x$num.pts)
         }, 
         "2c" ={cat("\n Class 1 points generated according to",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$xmat)
           
           cat("\n Class 2 points generated according to",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$ymat)
           
           cat("\nNumber of Points:\n n1 = number of class 1 points\n n2 = number of class 2 points\n")
           print(x$num.pts)
         }, 
         cc = {
           cat("\n Data points labeled as cases according to",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$xmat)
           
           cat("\n Data points labeled as controls according to",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$ymat)
           
           cat("\nNumber of Points:\n n1 = number of cases\n n0 = number of controls\n")
           print(x$num.pts)
         },
         ref.gen={
           cat("\n Generated Points from",x$pat.desc ," \n (first 5 or fewer are printed) \n")
           print(x$xmat)
           
           cat("\n Reference Points in the Region 
      (first 5 or fewer are printed) \n")
           print(x$ymat)
           
           cat("\nNumber of Points:\n n2 = number of generated points according to the pattern\n n1 = number of reference points\n")
           print(x$num.pts)
         }
  )
} #end of the function
#'
########################
#'
#' @title Plot a \code{SpatPatterns} object
#'
#' @description Plots the points generated from the pattern (color coded for each class) together with the
#' study window
#'
#' @param x Object of class \code{SpatPatterns}.
#' @param asp A numeric value, giving the aspect ratio for y axis to x-axis y/x (default is \code{NA}),
#' see the official help for \code{asp} by typing "? asp".
#' @param xlab,ylab Titles for the x and y axes, respectively (default is \code{xlab}="x" and \code{ylab}="y").
#' @param \dots Additional parameters for \code{plot}.
#'
#' @return
#' None
#'
# #' @seealso \code{\link{print.SpatPatterns}}, \code{\link{summary.SpatPatterns}}, and \code{\link{print.summary.SpatPatterns}}
#'
#' @examples
#' #TBF
#'
#' @export
plot.SpatPatterns<-function (x, asp=NA,xlab="x",ylab="y",...)
{
  dtlab <-x$lab
  
  ref.pts<-x$ref.points #reference points
  gen.pts<-x$gen.points #generated points
  dat.pts<-x$dat.points
  
  Xlim<-x$xlimit
  Ylim<-x$ylimit
  
  xf<-.01*(Xlim[2]-Xlim[1])
  yf<-.01*(Ylim[2]-Ylim[1])
  
  switch(x$pat.type, 
         "1c" = {
           X1<-x$gen.points
           plot(rbind(X1),asp=asp,pch=16,col=1,lwd=2, xlab=xlab,ylab=ylab,main=x$mtitle,
                xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf),...)
         }, 
         "2c" ={
           X1<-gen.pts[dtlab==1,]
           X2<-gen.pts[dtlab==2,]
           plot(rbind(X1),asp=asp,pch=16,col=2,lwd=2, xlab=xlab,ylab=ylab,main=x$mtitle,
                xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf),...)
           points(rbind(X2),pch=1,lwd=2)
           points(rbind(x$init.cases),pch=4,col=c(2,1),lwd=2)
           
           myl <- levels(as.factor(dtlab))
           legend(Xlim[2]-9*xf,Ylim[2]-yf,legend=myl,col=c(2,1),pch=c(16,1),cex=1.2)
         }, 
         cc = {
           X1<-dat.pts[dtlab==1,]
           X0<-dat.pts[dtlab==0,]
           plot(rbind(X1),asp=asp,pch=16,col=2,lwd=2, xlab=xlab,ylab=ylab,main=x$mtitle,
                xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf),...)
           points(rbind(X0),col=1,pch=1,lwd=2)
           points(rbind(x$init.cases),pch=4,col=2,lwd=2)
           points(rbind(x$cont.cases),pch=1,col="blue",lwd=2)
           
           pch.vec<-c(1,16)
           ifelse(length(levels(as.factor(dtlab)))>1,{col.cc<-c(1,2);pch.cc=pch.vec},
                  {cc.ind<-dtlab[1]+1;col.cc<-cc.ind;pch.cc<-pch.vec[cc.ind]})
           # {
           myl <- levels(as.factor(dtlab))
           legend(Xlim[2]-9*xf,Ylim[2]-yf,legend=myl,col=col.cc,pch=pch.cc,cex=1.2)
           # }
           
         },
         ref.gen={
           X1<-gen.pts
           X2<-ref.pts
           plot(rbind(X1),asp=asp,pch=16,col=2,lwd=2, xlab=xlab,ylab=ylab,main=x$mtitle,
                xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf),...)
           points(rbind(X2),pch=4,lwd=2)
         }
  )
  
} #end of the function
#'
#'
###############################################
#Auxiliary functions for class Clusters
###############################################
#'
#' @title Print a \code{Clusters} object
#'
#' @description Prints the \code{call} of the object of class '\code{Clusters}'
#' and also the \code{type} (or description) of the pattern).
#'
#' @param x A \code{Clusters} object.
#' @param \dots Additional arguments for the S3 method '\code{print}'.
#'
#' @return
#' The \code{call} of the object of class '\code{Clusters}'
#' and also the \code{type} (or description) of the pattern).
#'
#' @seealso \code{\link{summary.Clusters}}, \code{\link{print.summary.Clusters}}, and \code{\link{plot.Clusters}}
#'
#' @examples
#' #TBF (to be filled)
#'
#' @export
print.Clusters <- function(x,...)
{
  if (!inherits(x, "Clusters"))
    stop("x must be of class \"Clusters\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
} #end of the function
#'
########################
#'
#' @title Return a summary of a \code{Clusters} object
#'
#' @description Returns the below information about the \code{object}:
#'
#' \code{call} of the function defining the \code{object}, the \code{type} of the pattern, \code{parameters} of the pattern,
#' study window, some sample points from the generated pattern, reference points (if any for the bivariate pattern),
#' and number of points for each class
#'
#' @param object Object of class \code{Clusters}.
#' @param \dots Additional parameters for \code{summary}.
#'
#' @return
#' The \code{call} of the object of class '\code{Clusters}', the \code{type} of the pattern, \code{parameters} of the pattern,
#' study window, some sample points from the generated pattern, reference points (if any for the bivariate pattern),
#' and number of points for each class
#'
# #' @seealso \code{\link{print.Clusters}}, \code{\link{print.summary.Clusters}}, and \code{\link{plot.Clusters}}
#'
#' @examples
#' #TBF
#' @export
summary.Clusters <- function(object,...)
{
  if (!inherits(object, "Clusters"))
    stop("object must be of class \"Clusters\"")
  
  typ <- object$type
  gen.pts<-object$gen.points #generated points
  
  nx<-min(5,nrow(gen.pts))
  Npts<-object$num.points
  
  res <- list(pat.desc=object$desc.pat,
              call=object$call,
              xmat=gen.pts[1:nx,],
              param=object$parameters,
              type=typ,
              num.pts=Npts,
              Xlim=object$xlimit,
              Ylim=object$ylimit,
              gen.pts=gen.pts
  )
  
  class(res) <- "summary.Clusters"
  res
} #end of the function
#'
########################
#'
#' @title Print a summary of a \code{Clusters} object
#'
#' @description Prints some information about the \code{object}.
#'
#' @param x	 object of class "\code{summary.Clusters}", generated by \code{summary.Clusters}.
#' @param \dots Additional parameters for \code{print}.
#'
#' @return
#' None
#'
#' @seealso \code{\link{print.Clusters}}, \code{\link{summary.Clusters}}, and \code{\link{plot.Clusters}}
#'
#' @export
print.summary.Clusters <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the Pattern:\n")
  print(x$type)
  
  if (!is.null(x$param))
  { cat("\nParameters of the Pattern:\n")
    print(x$param)
  }
  
  cat("\nStudy Window:\n")
  cat("range in x-coordinate =", x$Xlim,"\n")
  cat("range in y-coordinate =", x$Ylim,"\n")
  
  cat("\n Generated points from the pattern of",x$pat.desc ," \n (first 5 or fewer are printed) \n")
  print(x$xmat)
  
  cat("\nNumber of generated points according to the clustering pattern\n")
  print(x$num.pts)
} #end of the function
#'
########################
#'
#' @title Plot a \code{Clusters} object
#'
#' @description Plots the points generated from the pattern (color coded for each class) together with the
#' study window
#'
#' @param x Object of class \code{Clusters}.
#' @param asp A numeric value, giving the aspect ratio for y axis to x-axis y/x (default is \code{NA}),
#' see the official help for \code{asp} by typing "? asp".
#' @param xlab,ylab Titles for the x and y axes, respectively (default is \code{xlab}="x" and \code{ylab}="y").
#' @param \dots Additional parameters for \code{plot}.
#'
#' @return
#' None
#'
# #' @seealso \code{\link{print.Clusters}}, \code{\link{summary.Clusters}}, and \code{\link{print.summary.Clusters}}
#'
#' @examples
#' #TBF
#'
#' @export
plot.Clusters<-function (x, asp=NA,xlab="x",ylab="y",...)
{
  gen.pts<-x$gen.points #generated points
  
  Xlim<-x$xlimit
  Ylim<-x$ylimit
  
  xf<-.01*(Xlim[2]-Xlim[1])
  yf<-.01*(Ylim[2]-Ylim[1])
  
  plot(gen.pts,asp=asp,pch=16,col=1,lwd=2, xlab=xlab,ylab=ylab,main=x$mtitle,
       xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf),...)
} #end of the function
#'