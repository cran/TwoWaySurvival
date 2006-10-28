\name{TwoWaySurvfitControl}
\alias{TwoWaySurvfitControl}
\title{Control Values for Fitting of two-way Hazards}
\description{The values supplied in the call of 'TwoWaySurvfitCreate' function will replace the defaults, and
a list with all possible arguments is returned. The returned list
is used as the 'control' argument to the 'TwoWaySurvfitCreate' function.
}
\usage{
TwoWaySurvfitControl(niter.epoch = 100, niter.penalty = 2, tol.epoch.theta = 1e-08, tol.epoch.variance = 1e-08,
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
  \item{print.log.lik}{logical value for printing the marginal log-likelihoods in each iteration.}
  \item{print.estimates}{logical value for printing estimates of the fixed parameters and penalty from each iteration.}
  \item{method}{specified optimization method for smoothing parameter; can be either "fix" menaning optimization via fix point 
                iteration, or "NR" meaning Newton-Raphson optimization. We advice to set method="fix" because
		 evaluations then seem to go faster.}
  \item{number.int}{a specified number of the integration points for approximation of the cumulative hazard function}		 
  \item{\dots}{other parameters which can only be 'num.knots.t' and 'num.knots.b' for the number of spline knots for 
  survival time and birth time accordingly, as well as starting values 'start.penalty.t' and 'start.penalty.b' for penalized terms 
  for the random parts of the model. The latter should be applied together, otherwise both of them will be initialy set to 1 internaly.
  }
}
\details{The defaults or user specified values are applied as the 'control' argument in the
call of the 'TwoWaySurvfitCreate' function. It can be an (empty) list object or a call to the 'TwoWaySurvfitControl' function itself,
whether or not with supplied arguments to be changed from their default values. If not supplied,
the values for 'num.knots.t' or 'num.knots.b' will be determined internally. The value of 'number.int' should be chosen moderately; this because of
generation of artificial poisson distributed data, for each observation respectively. (See suplemented paper describing the generation mechanism.)
}
\value{a list with components for each of the possible arguments.
}
\references{Kauerman G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients.
\emph{Computational Statistics and Data Analysis}, in press.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{TwoWaySurvfitCreate}}}
\keyword{list}

