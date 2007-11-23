\name{unemployed}
\docType{data}
\alias{unemployed}
\title{Sample of german unemployed}
\description{This data give the time spent in the unemployment spell and the calendar time of entry in that spell for a sample of german
unemployed individuals. The full data have been resived from the \emph{Deutsches Institut fuer Wirtschaftsforschung} 
and were analysed in the larger study. The results can be found in the supplemented paper.
}
\usage{data(unemployed)
}
\format{A data frame containing 300 sampled observations. Description of columns:
\itemize{
\item{birth}{Calendar time of entry in the unemployment spell, measured in months, ranging from the year 1983 till 2000.}
\item{surv}{Time spent in the unemployment spell, measured in months. The spells which lasted more than 36 months have been
 truncated on 36 months and denoted as sensored.}
\item{status}{Censoring indicator in the usual notation for right censored data.}
\item{nation}{Nationality of the unemployed, german vs. others.}
\item{gender}{Gender of the person.}
\item{age}{Age has been categorized in young (till 25), middleage and old (over 50).}
\item{training}{Whether an individual received a professional training.}
\item{university}{Whether an individual has a university degree.}
}
}

\source{Socio Economic Panel (SOEP), s. \url{www.diw.de/deutsch/sop}.}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients, in press.
}
\keyword{datasets}

