\name{print.TwoWaySurvfit}
\alias{print.TwoWaySurvfit}
\title{Prints an Object of class 'TwoWaySurvfit'}
\description{Prints estimates of fixed parameters of smoothing varying coefficients
as well as penalty values for their random parts.
}
\usage{print.TwoWaySurvfit(x,...)
}
\arguments{
  \item{x}{  object of class 'TwoWaySurvfit'}
  \item{...}{  additional parameters for print}
}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{summary.TwoWaySurvfit}}, \code{\link{plot.TwoWaySurvfit}}}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
birth.time<-c(1,11,23,4,23,5,7,12,5,12,23,12,15,21,4,6,1,4,13,11)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.TwoWaySurv.object<-TwoWaySurv(surv.time=surv.time,birth.time=birth.time,status=status)
x<-TwoWaySurvfitCreate(my.TwoWaySurv.object)
print.TwoWaySurvfit(x)
}
\keyword{print}