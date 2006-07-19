\name{TwoWaySurvfitCreate}
\alias{TwoWaySurvfitCreate}
\title{Additiv Two-Way Hazards Modelling with Varying Coefficients
}
\description{This function is for the fitting of models with hazard depending on both survival time and
calendar (here birth or entry in the study) time. The effects of covariates are modelled as varying coefficients, additiv in both times. 
For flexible smoothing of the effects truncated polynomial splines are employed.
}
\usage{TwoWaySurvfitCreate(formula=formula, data=data, na.action=na.fail, control=control)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of
model specification are given below.}
  \item{data}{an optional data frame containing variables in the model. If not supplied the handling of missing values will be accordingly 
        to the value of 'na.action'.}
  \item{na.action}{an optional parameter for handling with missing data}
  \item{control}{an optional list argument with control parameters needed for optimization routine. If not specified,
  default values from 'TwoWaySurvfitControl' will be taken; otherwise the values, which are defined by user, will
  overwrite the defaults.}
}
\details{
A formula argument can be specified in differentt ways. The response or left part of the formula can be either an object of class 'TwoWaySurv' or it 
can be a call to the 'TwoWaySurv' function with its arguments specified in it. The right part may be left unspecified (without '\eqn{\sim}{~}')
or it can be of '\eqn{\sim}{~}1' type. In this case no covariates are supplied to the model, hence only baseline hazard will be adopted. 
In the case with covariates these must be supplied explicitly through their names in a standard way. In all cases the 'formula' argument 
must be supplied.

If factor variables are to be considered, their reference categories should be specified bevor applied for optimization. 

A 'control' argument can be absent, in this case the defaults from the 'TwoWaySurvfitControl'
will be taken, or it can be either a list with some or all components of it specified by their names with
according values, or it can be a call to the 'TwoWaySurvfitControl' function with parameters specified in it.

The most of computation and optimization is implemented in the "internal" function 'TwoWayVarCoefBackfit', which is 
in the body of 'TwoWaySurvfitCreate'. Its code is well commented and can be clearly understood following the supplemented paper.

The working spline bases are ones consisted from truncated polynomials. Other alternativ could be B-splines, but has not
been implemented yet and might be considered in the following update of the package.
}
\value{An object of class 'TwoWaySurvit'. It has methods for 'print', 'summary' and 'plot'.
}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{TwoWaySurv}}, \code{\link{TwoWaySurvfitControl}}, \code{\link{TwoWaySurvfitObject}}
}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
birth.time<-c(1,11,23,4,23,5,7,12,5,12,23,12,15,21,4,6,1,4,13,11)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.TwoWaySurv.object<-TwoWaySurv(surv.time=surv.time,birth.time=birth.time,status=status)
my.TwoWaySurvit.object<-TwoWaySurvfitCreate(formula=my.TwoWaySurv.object, na.action=na.fail, 
                  control = list(niter.epoch=50,print.epoch=FALSE,method="fix"))
x1<-rnorm(20)
x2<-rt(20,2)
my.TwoWaySurvit.object<-TwoWaySurvfitCreate(formula = my.TwoWaySurv.object~x1+x2,data=data.frame(x1,x2),
                  na.action = na.fail,control = TwoWaySurvfitControl(niter.epoch=50,print.epoch=FALSE,method="fix"))
}
\keyword{survival}
\keyword{nonparametric}
\keyword{smooth}