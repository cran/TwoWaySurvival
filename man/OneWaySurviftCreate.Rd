\name{OneWaySurvfitCreate}
\alias{OneWaySurvfitCreate}
\title{Additiv One-Way Hazards Modelling with Varying Coefficients
}
\description{This function is for the fitting of models with hazard depending on a survival time only. Other time scales are not considered
(compare for this with 'TwoWaySurvfitCreate' function). The effects of covariates are modelled as varying coefficients.
For flexible smoothing of the effects whether a basis of truncated polynomial splines or the B-splines can be chosen.
}
\usage{OneWaySurvfitCreate(formula=formula, data=data, na.action=na.fail, control=control, basis=basis)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of
model specification are given below.}
  \item{data}{an optional data frame containing variables in the model. If not supplied the handling of missing values will be accordingly 
        to the value of 'na.action'.}
  \item{na.action}{an optional parameter for handling with missing data.}
  \item{control}{an optional list argument with control parameters needed for optimization routine. If not specified,
  default values from 'OneWaySurvfitControl' will be taken; otherwise the values, which are defined by user, will
  overwrite the defaults.}
  \item{basis}{character string needed for specifying the function basis, which will be chosen for estimating the varying coefficients. 
  Can only be on of 'trunc' or 'bspline'}
}
\details{
A formula argument can be specified in different ways. The response or left part of the formula can be either an object of class 'OneWaySurv' or it 
can be a call to the 'OneWaySurv' function with its arguments specified in it. The right part may be left unspecified (without '\eqn{\sim}{~}')
or it can be of '\eqn{\sim}{~}1' type. In this case no covariates are supplied to the model, hence only baseline hazard will be adopted. 
In the case with covariates these must be supplied explicitly through their names in a standard way. In all cases the 'formula' argument 
must be supplied.

If factor variables are to be considered, their reference categories should be specified bevor applied for optimization. 

A 'control' argument can be absent, in this case the defaults from the 'OneWaySurvfitControl'
will be taken, or it can be either a list with some or all components of it specified by their names with
according values, or it can be a call to the 'OneWaySurvfitControl' function with parameters specified in it. 

The most of computation and optimization is implemented in the "internal" functions 'TruncatedOneway' or 'BsplineOneway', which are called from
the body of 'OneWaySurvfitCreate' in depending on the bases, which are chosen for the fitting routine.
The code for both routines is well commented and can be clearly understood following the supplemented paper.

The computation procedure can take considerable time, depending on data supplied. All evaluations are written in R-code, so no external
program code like C or Fortran has been used. We advice initialy to set the number of iterations not to be too large, and eventually increase it
if needed.
}
\value{An object of class 'OneWaySurvit'. It has methods for 'print', 'summary' and 'plot'.
}
\references{Kauermann G.(2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{OneWaySurv}}, \code{\link{OneWaySurvfitControl}}, \code{\link{OneWaySurvfitObject}}
}
\examples{surv.time<-c(1,2,4,3,5,6,1,6,2,3,4,7,5,1,6,7,10,2,10,5)
status<-c(1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1)
my.OneWaySurv.object<-OneWaySurv(surv.time=surv.time,status=status)
my.OneWaySurvfit.object<-OneWaySurvfitCreate(formula=my.OneWaySurv.object, na.action=na.fail, 
                         control = list(niter.epoch=50,print.epoch=FALSE,method="fix"),basis="bspline")
x1<-rnorm(20)
x2<-rt(20,2)
my.OneWaySurvfit.object<-OneWaySurvfitCreate(formula = my.OneWaySurv.object~x1+x2,data=data.frame(x1,x2),
                         na.action = na.fail,control = OneWaySurvfitControl(niter.epoch=50,print.epoch=FALSE,method="fix"))
}
\keyword{survival}
\keyword{nonparametric}
\keyword{smooth}

