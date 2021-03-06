% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M4B.r
\name{M4B}
\alias{M4B}
\title{Model M4B MLE function}
\usage{
M4B(data, BadStart, nchains, method, niter, lambda.min)
}
\arguments{
\item{data}{A vector of class 'circular'}

\item{BadStart}{An integer to replace the log likelihood when starting
parameters of the optimizer fall outside the preset bounds.  This is
usually set to a large number. Can also be set to Inf, but will
result in an error if a method other than "Nelder-Mead" is chosen. (default = 10^9)}

\item{nchains}{A positive integer indicating the number of chains to run.
Only the chain with the lowest log likelihood is returned (default = 5)}

\item{method}{A character string indicating the optimizing algorithm
to use.  Either "BFGS" or "Nelder-Mead" are recommended. See ?optim
for more details. (default = "BFGS")}

\item{niter}{The maximum number of iterations for the optimizing
algorithm.  Equivalent to the 'maxit' control parameter of the
optim() function.  See ?optim for more details. (default = 5000)}

\item{lambda.min}{The minimum proportional size of the first distribution.
Must be between 0 and 1. (default = 0.25)}
}
\value{
A list with elements (same as for function optim()):

$par:  Vector with the optimized mean angle (mu1),
concentration parameter (kappa1), concentration parameter (kappa2),
and proportional size (lambda) of the first distribution.

$lik:  The negative log likelihood

$counts:  A two-element integer vector giving the number of calls to
‘fn’ and ‘gr’ respectively. See ?optim() for details.

$convergence:  An integer code. ‘0’ indicates successful completion
(which is always the case for ‘"SANN"’ and ‘"Brent"’). Possible error codes are:

‘1’ indicates that the iteration limit ‘maxit’ had been reached.

‘10’ indicates degeneracy of the Nelder-Mead simplex.

‘51’ indicates a warning from the ‘"L-BFGS-B"’ method; see
component ‘message’ for further details.

‘52’ indicates an error from the ‘"L-BFGS-B"’ method; see
component ‘message’ for further details.

$message:  A character string giving any additional information returned
by the optimizer, or ‘NULL’.
}
\description{
Run Maximum likelihood estimation for model M4B.
}
\examples{
testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
M4B(testdata)
}
\keyword{M4B}
