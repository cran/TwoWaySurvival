\name{plot.OneWaySurvfit}
\alias{plot.OneWaySurvfit}
\title{Plots the Smoothed Varying Coefficients}
\description{'plot.OneWaySurvfit' is a proposed function for plotting of smooth components of the one-dimensional model,
one plot for each baseline and covariate (or factor level) effect, for survival time.
}
\usage{
\method{plot}{OneWaySurvfit}(x,...)
}
\arguments{
 \item{x}{object of class 'OneWaySurvfit'}
 \item{...}{additional plot parameters}
}
\details{All plots will be made in the device, which is specified by the user. 
One can use the object components for producing his own plots. 
}
\references{Kauermann G. (2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis},  49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{print.OneWaySurvfit}}, \code{\link{summary.OneWaySurvfit}}}
\examples{pdf(file="myplot.pdf",h=12,w=12)
par(mfrow=c(1,2))
surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.OneWaySurv.object<-OneWaySurv(surv.time=surv.time,status=status)
x<-OneWaySurvfitCreate(my.OneWaySurv.object~1)
plot.OneWaySurvfit(x)
dev.off()
}
\keyword{hplot}

