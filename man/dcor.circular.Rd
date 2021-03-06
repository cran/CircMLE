% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dcor.circular.r
\name{dcor.circular}
\alias{dcor.circular}
\title{circular distance correlation function}
\usage{
dcor.circular(x, y, method = "chord", type = "c-c", ...)
}
\arguments{
\item{x}{A vector of class 'circular', or numeric vector of angles measured in radians}

\item{y}{A vector of class 'circular', numeric vector of angles measured in radians, or numeric vector}

\item{method}{the distance measure to be used. This must be one of the following functions:
‘"angularseparation"’, ‘"chord"’, '"geodesic"’, or '"circ.range"' (default = "chord").
see ?dist.circular for additional details.}

\item{type}{if ‘type == "c-c"’ then perform a circular-circular distance
corellation, else if ‘type == "c-l"’ then perform a circular-linear distance
corellation (default = "c-c").}

\item{...}{additional parameters passed to the dcor.test function}
}
\value{
Same as from the \code{\link[energy]{dcor.test}} function:
a list with class ‘htest’containing

method: description of test

statistic: observed value of the test statistic

estimate: dCov(x,y) or dCor(x,y)

estimates: a vector: [dCov(x,y), dCor(x,y), dVar(x), dVar(y)]

replicates: replicates of the test statistic

p.value: approximate p-value of the test

n: sample size

data.name: description of data
}
\description{
Perform a distance correlation between circular datasets or between circular and linear datasets.
}
\examples{
# Circular-circular distance corellation
x <- circular::rvonmises(n = 50, mu = circular::circular(0), kappa = 3)
y <- x + circular::rvonmises(n = 50, mu = circular::circular(pi), kappa = 10)
dcor.circular(x, y)

# Run permutation test with 9999 iterations
dcor.circular(x, y, R = 9999)

# Circular-linear distance corellation
x <- circular::rvonmises(n = 50, mu = circular::circular(0), kappa = 3)
y <- as.numeric(x) + rnorm(50, mean = 5, sd = 2)
dcor.circular(x, y, type = "c-l", R = 9999)
}
\seealso{
\code{\link[energy]{dcor}}
\code{\link[energy]{dcov}}
\code{\link[energy]{DCOR}}
\code{\link[energy]{dcor.test}}
\code{\link[circular]{dist.circular}}
}
\keyword{correlation}
\keyword{distance}
