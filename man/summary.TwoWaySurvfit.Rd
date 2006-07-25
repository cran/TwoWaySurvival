\name{summary.TwoWaySurvfit}
\alias{summary.TwoWaySurvfit}
\title{Objects Summaries}
\description{prints detailed results from the fitting of the two-way hazards modelling.
}
\usage{summary.TwoWaySurvfit(object,...)
}
\arguments{
  \item{object}{an object of class 'TwoWaySurvfit'}
  \item{...}{additional arguments}
}
\value{the printed values are estimates of fixed and random components of varying coefficients with
their estimated variances as well as penalty values and marginal log-likelihoods at start point and at optimum.
}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{print.TwoWaySurvfit}}, \code{\link{plot.TwoWaySurvfit}}}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
birth.time<-c(1,11,23,4,23,5,7,12,5,12,23,12,15,21,4,6,1,4,13,11)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.TwoWaySurv.object<-TwoWaySurv(surv.time=surv.time,birth.time=birth.time,status=status)
x<-TwoWaySurvfitCreate(my.TwoWaySurv.object)
summary.TwoWaySurvfit(x)
}
\keyword{print}

