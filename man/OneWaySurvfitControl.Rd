\name{OneWaySurvfitControl}
\alias{OneWaySurvfitControl}
\title{Control Values for Fitting of two-way Hazards}
\description{The values supplied in the call of 'OneWaySurvfitCreate' function will replace the defaults, and
a list with all possible arguments is returned. The returned list
is used as the 'control' argument to the 'OneWaySurvfitCreate' function.
}
\usage{
OneWaySurvfitControl(niter.epoch = 100, niter.penalty = 2, tol.epoch.theta = 1e-08, tol.epoch.variance = 1e-08,
                     tol.penalty = 1e-08, print.epoch = FALSE, print.penalty=FALSE, print.log.lik=FALSE, 
		     print.estimates=FALSE, method="fix", number.int=60, ...)
}
\arguments{
  \item{niter.epoch}{maximum number of iterations of the outer loop in the backfitting optimization algorithm.}
  \item{niter.penalty}{maximum number of iterations of the inner loop in the backfitting optimization algorithm.}
  \item{tol.epoch.theta}{tolerance for the convergence criterion for the fixed and random parameters.}
  \item{tol.epoch.variance}{tolerance for the convergence criterion for the penalty values w.r.t. outer loop.}
  \item{tol.penalty}{tolerance for the convergence criterion for the penalty values w.r.t. inner loop.}
  \item{print.epoch}{logical value for printing results from each outer loop.}
  \item{print.penalty}{logical value for printing results from each inner loop.}
  \item{print.log.lik}{logical value for printing the marginal log-likelihood in each iteration.}
  \item{print.estimates}{logical value for printing estimates of the fixed parameters and penalty from each iteration.}
  \item{method}{character string specifying optimization method for smoothing parameter; currently can be only "fix" menaning optimization via fix point 
                iteration.}
  \item{number.int}{a specified number of the integration points for approximation of the cumulative hazard function}		 
  \item{\dots}{other parameters which can only be 'num.knots.t' for the number of spline knots for 
  survival time, as well as starting value 'start.penalty.t' for penalized terms 
  for the random parts of the model.
  }
}
\details{The defaults or user specified values are applied as the 'control' argument in the
call of the 'OneWaySurvfitCreate' function. It can be an (empty) list object or a call to the 'OneWaySurvfitControl' function itself,
whether or not with supplied arguments to be changed from their default values. If not supplied,
the value for 'num.knots.t' will be determined internally. The value of 'number.int' should be chosen moderately; this because of
generation of artificial poisson distributed data, for each observation respectively. (See suplemented paper describing the generation mechanism.)
}
\value{a list with components for each of the possible arguments.
}
\references{Kauerman G. (2005). Penalised Spline Fitting in Multivariable Survival Models with Varying Coefficients
\emph{Computational Statistics and Data Analysis}, 49, 169-186.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{OneWaySurvfitCreate}}}
\keyword{list}

