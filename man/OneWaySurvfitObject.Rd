\name{OneWaySurvfitObject}
\alias{OneWaySurvfitObject}
\title{Fitted one-dimensional Object}
\description{An object of class 'OneWaySurvfit' returned as a list by the 'OneWaySurvfitCreate' function
and representing a fitted one-way hazard model with varying coefficients. It has methods for generic functions print, plot and summary.
}
\value{
  \item{fix.coef}{named vector of not penalized parameters of the model. It has the value 'NULL' if B-spline bases are chosen for modelling.}
  \item{random.coef}{named vector of penalized parameters of the model.}
  \item{penalty}{named vector of penalty values for the random parts of the model.}
  \item{var.fix}{estimated variances of not penalized parameters of the model. It has the value 'NULL' if B-spline bases are chosen for modelling.}
  \item{var.random}{estimated variances of the predicted values of (penalized) random components of the model.}
  \item{log.lik.margin.start}{initial value of the marginal log-likelihood.}
  \item{log.lik.margin}{the value of the marginal log-likelihood at optimum.}
  \item{df}{estimated degrees of freedom for the parts of the model.}
  \item{df.total}{total or summed estimated degrees of freedom.}
  \item{niter.epoch}{number of outer loops in optimization routine.}
  \item{varying.frame}{data frame with named columns, according to the names of the varying coefficients.
                       It is needed for plotting of the smooth components.}
  \item{deviation.frame}{data frame with named columns, according to the names of the varying coefficients.
                         It is needed for plotting of the confidence bands of the smooth components.}
  \item{grid.frame}{data frame with just one column, for survival time, needed for plotting issues.}
  \item{p}{number of covariates; for factors including their categories (excluding reference category).}
  \item{factor.names}{covariate names; for factors the names of categories (excluding reference category).}
  } 
\references{Kauerman G. (2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{OneWaySurv}}, \code{\link{OneWaySurvfitCreate}}}
\keyword{misc}  

