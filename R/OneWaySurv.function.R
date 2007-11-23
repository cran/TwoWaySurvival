"OneWaySurv" <-
function(surv.time, status)
{
 
if (!is.numeric(surv.time))
    stop("not numeric time argument")
if (!is.numeric(status))
    stop("censoring must be binary variable: 1 for death, 0 for censored")
  
    Surv.object <- cbind(surv.time, status)
    dimnames(Surv.object) <- list(NULL, c("time","status"))
    class(Surv.object) <- c("OneWaySurv","matrix")
    Surv.object
}

