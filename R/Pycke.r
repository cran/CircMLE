#' Pycke test function
#'
#' Perform the Pycke test.
#' @param data A vector of class 'circular', or numeric vector of 
#' angles measured in radians
#' @param iter The number of bootstrap replicates to perform in order 
#' to estimate the p-value of the test. (default = 9999)
#' @keywords Pycke test
#' @return A numeric vector of the test statistic (T) and associated p-value
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(20, mu = circular::circular(pi), kappa = 3)
#' pycke_test(testdata, iter = 999)

pycke_test <- function(data, iter = 9999){

	#Check input parameters
	if (is.null(data) | missing(data)) stop("Please provide input data vector")
	if (missing(iter)) iter = 9999 else iter = iter
	if (length(iter) != 1 | !is.numeric(iter)) stop("Please set the number of iterations correctly for p-value estimation")
	if (iter <= 99) warning("WARNING: You set an awfully small number of iterations for p-value estimation, it is recommended to run at least 999 for accurate estimation.")
	
	# Coerce into correct format
	data = check_data(data)
	
	# Get all circular attributes from data
    params = circularp(data)
	
	# Begin Pycke test	
	n <- length(data)
	testset <- rep(0, iter)
	for (f in 1:iter){
		data1 <- matrix(rcircularuniform(n, control.circular = params))
		testset[f] <- PyckeT(data1)
	}
	Tsample <- PyckeT(data)
	counter <- 1
	for(j in 1:iter){
		if(testset[j] >= Tsample){
			counter <- counter + 1
		}
	}
	p <- counter / (iter + 1)
	out <- c(as.numeric(Tsample), p)
	names(out) <- c("Test statistic (T)", "p-value")
	return(out)
}

# Pycke test formula, as reported by Landler et al. 2019
PyckeT <- function(data){
	n <- length(data)
	total <- 0
	for (i in 2:n){
		for (j in 1:(i - 1)){
			numerator <- cos(data[i] - data[j]) - sqrt(0.5)
			denominator <- 1.5 - (2 * sqrt(0.5) * cos(data[i] - data[j]))
			total <- total + (numerator / denominator)
		}
	}
	T <- (2 * total) / (n - 1)
	return(T)
}