#' Hermans-Rasson test function
#'
#' Perform variants of the Hermans-Rasson test.
#' @param data A vector of class 'circular', or numeric vector of 
#' angles measured in radians
#' @param original A logical of whether or not to run the original 
#' version of the Hermans-Rasson test or the newer version described
#' in Landler et al. (2019) doi: 10.1186/s12898-019-0246-8 (default = F)
#' @param iter The number of bootstrap replicates to perform in order 
#' to estimate the p-value of the test. (default = 9999)
#' @keywords Hermans-Rasson test
#' @return A numeric vector of the test statistic (T) and associated p-value
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(20, mu = circular::circular(pi), kappa = 3)
#' HR_test(testdata, iter = 999)

HR_test <- function(data, original = F, iter = 9999){

	#Check input parameters
	if (is.null(data) | missing(data)) stop("Please provide input data vector")
	if (missing(original)) original = F else original = original
	if (!is.logical(original)) stop ("Please set \"original\" to TRUE (T) or FALSE (F)")
	if (missing(iter)) iter = 9999 else iter = iter
	if (length(iter) != 1 | !is.numeric(iter)) stop("Please set the number of iterations correctly for p-value estimation")
	if (iter <= 99) warning("WARNING: You set an awfully small number of iterations for p-value estimation, it is recommended to run at least 999 for accurate estimation.")
	
	# Coerce into correct format
	data = check_data(data)
	
	# Get all circular attributes from data
    params = circularp(data)
	
	# Begin Hermans-Rasson test
	n <- length(data)
	# sample <- circular(data, units = "radians")  # not needed because of check_data() above
	testset<- rep(0, iter)
	for (f in 1:iter){
		data1 <- matrix(rcircularuniform(n, control.circular = params))
		if(original) testset[f] <- HermansRassonT(data1) else testset[f] <- HermansRasson2T(data1)
	}
	
	# Run appropriate method for T
	if(original) Tsample <- HermansRassonT(data) else Tsample <- HermansRasson2T(data)
	counter <- 1
	for(j in 1:iter){
		if(testset[j] > Tsample){
			counter <- counter + 1
		}
	}
	p <- counter/(iter + 1)
	out <- c(as.numeric(Tsample), p)
	names(out) <- c("Test statistic (T)", "p-value")
	return(out)
}

# Original Hermans-Rasson test formula, referred to as HR-infinity by Landler et al. 2019
HermansRassonT <- function(data){
	n <- length(data)
	total <- 0
	for (i in 1:n){
		for (j in 1:n){
			total <- total + abs(sin(data[i] - data[j]))
		}
	}
	T <- abs((n / pi) - (total / (2 * n)))
	return(T)
}

# Updated Hermans-Rasson test formula, referred to as HR by Landler et al. 2019
HermansRasson2T <- function(data){
	n <- length(data)
	total <- 0
	for (i in 1:n){
		for (j in 1:n){
			total <- total + abs(abs(data[i] - data[j]) - pi) - (pi / 2)
			total <- total - (2.895 * (abs(sin(data[i] - data[j])) - (2 / pi)))
		}
	}
	T <- total / n
	return(T)
}