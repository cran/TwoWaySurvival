\name{plot.TwoWaySurvfit}
\alias{plot.TwoWaySurvfit}
\title{Plots the Smoothed Varying Coefficients}
\description{'plot.TwoWaySurvfit' is a proposed function for plotting of smooth components of the model,
one plot for each baseline and covariate (or factor level) effect, for survival time and birth time accordingly.
}
\usage{plot.TwoWaySurvfit(x,...)
}
\arguments{
 \item{x}{object of class 'TwoWaySurvfit'}
 \item{...}{additional plot parameters}
}
\details{All plots will be made in the device, which is specified by the user. 
One can use the object components for producing his own plots. 
}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{print.TwoWaySurvfit}}, \code{\link{summary.TwoWaySurvfit}}}
\examples{pdf(file="~/my.plot.pdf",h=12,w=12)
par(mfrow=c(2,2))
surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
birth.time<-c(1,11,23,4,23,5,7,12,5,12,23,12,15,21,4,6,1,4,13,11)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.TwoWaySurv.object<-TwoWaySurv(surv.time=surv.time,birth.time=birth.time,status=status)
x<-TwoWaySurvfitCreate(my.TwoWaySurv.object~1)
plot.TwoWaySurvfit(x)
dev.off()
}
\keyword{hplot}