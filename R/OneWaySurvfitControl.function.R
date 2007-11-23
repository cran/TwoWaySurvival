"OneWaySurvfitControl" <-
function(niter.epoch=100,niter.penalty=2,tol.epoch.theta=1e-08,tol.epoch.variance=1e-08,tol.penalty=1e-08,print.epoch=FALSE,print.penalty=FALSE,print.log.lik=FALSE,print.estimates=FALSE,method="fix",number.int=60,...)
{
list(niter.epoch=100,niter.penalty=2,tol.epoch.theta=1e-08,tol.epoch.variance=1e-08,tol.penalty=1e-08,print.epoch=FALSE,print.penalty=FALSE,print.log.lik=FALSE,print.estimates=FALSE,method="fix",number.int=60,...) 
}

#... can be some or all of 
#num.knots.t=num.knots.t (specified number of knots for splines in surv.time-direction)
#start.penalty.t=start.penalty.t  (specified initialization vector of length=number.of.factors +1 for penalties in t-direction; default is vector of 1's)