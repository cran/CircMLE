#' circular distance correlation function
#'
#' Perform a distance correlation between circular datasets or between circular and linear datasets.
#' @param x A vector of class 'circular', or numeric vector of angles measured in radians
#' @param y A vector of class 'circular', numeric vector of angles measured in radians, or numeric vector
#' @param method the distance measure to be used. This must be one of the following functions:
#' ‘"angularseparation"’, ‘"chord"’, '"geodesic"’, or '"circ.range"' (default = "chord").
#' see ?dist.circular for additional details.
#' @param type if ‘type == "c-c"’ then perform a circular-circular distance
#' corellation, else if ‘type == "c-l"’ then perform a circular-linear distance
#' corellation (default = "c-c").
#' @param ... additional parameters passed to the dcor.test function
#' @keywords distance correlation
#' @return Same as from the \code{\link[energy]{dcor.test}} function:
#' a list with class ‘htest’containing
#'
#' method: description of test
#'
#' statistic: observed value of the test statistic
#'
#' estimate: dCov(x,y) or dCor(x,y)
#'
#' estimates: a vector: [dCov(x,y), dCor(x,y), dVar(x), dVar(y)]
#'
#' replicates: replicates of the test statistic
#'
#' p.value: approximate p-value of the test
#'
#' n: sample size
#'
#' data.name: description of data
#' @import circular
#' @import energy
#' @export
#' @seealso
#' \code{\link[energy]{dcor}}
#' \code{\link[energy]{dcov}}
#' \code{\link[energy]{DCOR}}
#' \code{\link[energy]{dcor.test}}
#' \code{\link[circular]{dist.circular}}
#' @examples
#' # Circular-circular distance corellation
#' x <- circular::rvonmises(n = 50, mu = circular::circular(0), kappa = 3)
#' y <- x + circular::rvonmises(n = 50, mu = circular::circular(pi), kappa = 10)
#' dcor.circular(x, y)
#'
#' # Run permutation test with 9999 iterations
#' dcor.circular(x, y, R = 9999)
#'
#' # Circular-linear distance corellation
#' x <- circular::rvonmises(n = 50, mu = circular::circular(0), kappa = 3)
#' y <- as.numeric(x) + rnorm(50, mean = 5, sd = 2)
#' dcor.circular(x, y, type = "c-l", R = 9999)

dcor.circular <- function(x, y, method = "chord", type = "c-c", ...){

	#Check input parameters
	if (is.null(x) | missing(x)) stop("Please provide input data vectors")
	if (is.null(y) | missing(y)) stop("Please provide input data vectors")
	if (missing(method)) method = "chord" else method = method
	if (missing(type)) type = "c-c" else type = type
	if (length(method) != 1 | length (type) != 1 | !is.character(method) | !is.character(type)) stop("Please set the input parameters properly.")
	funcs <- list("angularseparation", "chord", "geodesic", "circ.range")
	if (! method %in% funcs) stop("\"method\" not specified properly.")
	if (! type %in% c("c-c", "c-l")) stop("\"type\" not specified properly.")
	
	# Coerce into correct format
	x <- check_data(x)
	if (type == "c-c") y <- check_data(y)
	
	# Get all circular attributes from data
    paramsx <- circularp(x)
    if (type == "c-c") paramsy <- circularp(y)
	
	# Make sure lengths are equal!
	if(length(x) != length(y)) stop("Error: x and y must have the same length.")
	
	# Get the method to use	
	f <- get(funcs[[which(funcs %in% method)]])
	
	# Get x pairwise distances
	xdist <- customDist(x, f)
	
	# Get y pairwise distances
	if (type == "c-c") ydist <- customDist(y, f)
	if (type == "c-l") ydist <- stats::dist(y)	

	# Calculate distance corellation and perform permutation test
	out <- dcor.test(xdist, ydist, ...)
	return(out)
}

#############################################################
############ Circular distance functions  ###################
#############################################################
customDist <- function(x, metric, diag = FALSE, upper = FALSE) {
  xlen <- length(x)
  dm <- matrix(nrow = xlen, ncol = xlen)
  for (j in 1:(xlen-1)) for (i in j:(xlen-1)) dm[i+1,j] <- metric(x[[i+1]], x[[j]])
  stats::as.dist(dm, upper = upper, diag = diag)
}

# Angular separation
angularseparation <- function(x1, x2){
	return(1 - cos(x1 - x2))
}

# Chord distance
chord <- function(x1, x2){
	return(sqrt(2 * (1 - cos(x1 - x2))))
}
  
# Geodesic distance
geodesic <- function(x1, x2){
	return(pi - abs(pi - abs(x1 - x2)))
}

# Circular range distance
circ.range <- function(x1, x2){
	d <- as.circular(c(x1, x2), modulo = "2pi", type = "angles", units = "radians", zero = 0, template = "none", rotation = "counter")
	return(as.numeric(range.circular(d)))
}