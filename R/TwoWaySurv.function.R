"TwoWaySurv" <-
function(surv.time, birth.time, status)
{
 
if (!is.numeric(surv.time) || !is.numeric(birth.time))
    stop("not numeric times arguments")
if (!is.numeric(status))
    stop("censoring must be binary variable: 1 for death, 0 for censored")
  
    Surv.object <- cbind(surv.time, birth.time, status)
    dimnames(Surv.object) <- list(NULL, c("time", "birth", "status"))
    class(Surv.object) <- c("TwoWaySurv","matrix")
    Surv.object
}

