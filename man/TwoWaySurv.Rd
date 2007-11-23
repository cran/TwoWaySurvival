\name{TwoWaySurv}
\alias{TwoWaySurv}
\title{Creates a Two-Way Survival Object
}
\description{Creates an object of class 'TwoWaySurv' for using as a response in the model formula
in the call to the 'TwoWaySurvitCreate' function.
}
\usage{
TwoWaySurv(surv.time, birth.time, status)
}
\arguments{
  \item{surv.time}{   the time to follow up for right censored data.}
  \item{birth.time}{ the (calendar) time to enter in the study.}
  \item{status}{   the (0,1)-indicator, 1 for event and 0 for censored observation.}
}
\details{The function handls only right censored data. All arguments must be 
numeric vectors. NA's can be supplied but will not be supported  by  the function 'TwoWaySurvfitCreate' in the further analysis.
}
\value{An object of class 'TwoWaySurv'.}
\references{Kauermann G. and Khomski P. (2006): Additiv two-way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, 51, 1944-1956.
}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{TwoWaySurvfitCreate}}}
\examples{
surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
birth.time<-c(1,11,23,4,23,5,7,12,5,12,23,12,15,21,4,6,1,4,13,11)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.TwoWaySurv.object<-TwoWaySurv(surv.time=surv.time,birth.time=birth.time,status=status)
}
\keyword{misc}

