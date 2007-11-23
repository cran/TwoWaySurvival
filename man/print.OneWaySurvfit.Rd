\name{print.OneWaySurvfit}
\alias{print.OneWaySurvfit}
\title{Prints an Object of class 'OneWaySurvfit'}
\description{Prints estimates of fixed parameters of smoothing varying coefficients
as well as penalty values for their random parts.
}
\usage{
\method{print}{OneWaySurvfit}(x,...)
}
\arguments{
  \item{x}{  object of class 'OneWaySurvfit'}
  \item{...}{  additional parameters for print}
}
\references{Kauerman G. (2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{summary.OneWaySurvfit}}, \code{\link{plot.OneWaySurvfit}}}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.OneWaySurv.object<-OneWaySurv(surv.time=surv.time,status=status)
x<-OneWaySurvfitCreate(my.OneWaySurv.object)
print.OneWaySurvfit(x)
}
\keyword{print}

