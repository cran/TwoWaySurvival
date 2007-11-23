\name{OneWaySurv}
\alias{OneWaySurv}
\title{Creates a One-Way Survival Object
}
\description{Creates an object of class 'OneWaySurv' for using as a response in the model formula
in the call to the 'OneWaySurvitCreate' function.
}
\usage{
OneWaySurv(surv.time, status)
}
\arguments{
  \item{surv.time}{   the time to follow up for right censored data.}
  \item{status}{   the (0,1)-indicator, 1 for event and 0 for censored observation.}
}
\details{The function handls only right censored data. All arguments must be 
numeric vectors. NA's can be supplied but will not be supported  by  the function 'OneWaySurvfitCreate' in the further analysis.
}
\value{An object of class 'OneWaySurv'.}
\references{Kauermann G.(2005): Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.
}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{OneWaySurvfitCreate}}}
\examples{
surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.OneWaySurv.object<-OneWaySurv(surv.time=surv.time,status=status)
}
\keyword{misc}

