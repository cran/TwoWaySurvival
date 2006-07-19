\name{TwoWaySurvfitObject}
\alias{TwoWaySurvfitObject}
\title{Fitted Object}
\description{An object of class 'TwoWaySurvfit' returned as a list by the 'TwoWaySurvfitCreate' function
and representing a fitted additiv two-way hazard model with varying coefficients. It has methods for generic functions print, plot and summary
}
\value{
  \item{fix.coef}{named vector of fixed parameters, beta, of the model}
  \item{random.coef}{named vector of random parameters, u, of the model}
  \item{penalty}{named vector of penalty values for the random part of the model}
  \item{var.fix}{estimated variances of the fixed components, beta, of the model}
  \item{var.random}{estimated variances of the random components, u, of the model}
  \item{log.lik.margin.start}{initial value of the marginal log-likelihoods for t- and b- direction}
  \item{log.lik.margin}{the value of the marginal log-likelihoods for t- and b- direction at optimum}
  \item{df}{estimated degrees of freedom for the parts of the model}
  \item{df.total}{total or summed estimated degrees of freedom}
  \item{niter.epoch}{number of outer loops in optimization routine}
  \item{varying.frame}{data frame with named columns, according to the names of the varying coefficients, 
  for each t- and b- direction. It is needed for plotting of the smooth components.}
  \item{deviation.frame}{data frame with named columns, according to the names of the varying coefficients, 
  for each t- and b- direction. It is needed for plotting of the confidence bands of the smooth components.}
  \item{grid.frame}{data frame with two columns, the first for survival time and the second for the birth time, needed
  for plotting issues.}
  \item{p}{number of covariates; for factors including their categories (excluding reference category)}
  \item{factor.names}{covariate names; for factors names of the categories (excluding reference category)}
  } 
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press
}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{TwoWaySurv}}, \code{\link{TwoWaySurvfitCreate}}}
\keyword{misc}  