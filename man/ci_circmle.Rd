% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_circmle.r
\name{ci_circmle}
\alias{ci_circmle}
\title{ci_circmle}
\usage{
ci_circmle(circmle, model)
}
\arguments{
\item{circmle}{A list consisting of the output from function 'circ_mle'}

\item{model}{character string indicating the model to be used to estimate parameter uncertainty.
Must be one of c("M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B",
"M5A", "M5B"). Default = the $bestmodel from the circmle object.}
}
\value{
A data frame with a column for the parameter name, the maximum likelihood estimate (MLE),
standard error (SE), and 95\% confidence interval (CI) for each estimated model parameter.
}
\description{
Calculate the 95\% confidence interval for estimated model parameters
}
\examples{
testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
fit <- circ_mle(testdata)
ci_circmle(fit)
}
\keyword{circ_mle}
