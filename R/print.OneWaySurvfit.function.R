"print.OneWaySurvfit" <-
function(x,...)
{
cat("fixed coefficients","\n\n")
print(x$fix.coef,...)
cat("\n","Penalties for random components","\n\n")
print(x$penalty,...)
}

