\name{summary.OneWaySurvfit}
\alias{summary.OneWaySurvfit}
\title{Objects Summaries}
\description{prints detailed results from the fitting of the one-way hazards modelling.
}
\usage{
\method{summary}{OneWaySurvfit}(object,...)
}
\arguments{
  \item{object}{an object of class 'OneWaySurvfit'}
  \item{...}{additional arguments}
}
\value{the printed values are estimates of fixed and random components of varying coefficients with
their estimated variances as well as penalty values and marginal log-likelihoods at start point and at optimum.
}
\references{Kauerman G. (2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{print.OneWaySurvfit}}, \code{\link{plot.OneWaySurvfit}}}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.OneWaySurv.object<-OneWaySurv(surv.time=surv.time,status=status)
x<-OneWaySurvfitCreate(my.OneWaySurv.object)
summary.OneWaySurvfit(x)
}
\keyword{print}

