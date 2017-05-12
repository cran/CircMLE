#' circ_mle
#'
#' Run all 10 maximum likelihood models of circular orientation
#' @param data A vector of class 'circular'
#' @param criterion chose from either "AIC", "AICc", or "BIC" for the
#' model comparison information criterion. (default = "AIC")
#' @param BadStart An integer to replace the log likelihood when starting
#' parameters of the optimizer fall outside the preset bounds.  This is
#' usually set to a large number. Can also be set to Inf, but will
#' result in an error if a method other than "Nelder-Mead" is chosen. (default = 10^9)
#' @param nchains A positive integer indicating the number of chains to run.
#' Only the chain with the lowest log likelihood is returned (default = 5)
#' @param method A character string indicating the optimizing algorithm
#' to use.  Either "BFGS" or "Nelder-Mead" are recommended. See ?optim
#' for more details. (default = "BFGS")
#' @param niter The maximum number of iterations for the optimizing
#' algorithm.  Equivalent to the 'maxit' control parameter of the
#' optim() function.  See ?optim for more details. (default = 5000)
#' @keywords circ_mle
#' @return A list with 3 elements:
#' @return $results: A data frame consisting of a row for each model (rownames) with the columns:
#' 1 = number of free parameters, 2 = mu1, 3 = kappa1, 4 = lamda, 5 = mu2, 6 = kappa2,
#' 7 = negative log likelihood, 8 = Counts.function, 9 = Counts.Gradient,
#' 10 = convergence, 11 = message. See ?optim for details on columns 8-10.
#' @return $bestmodel:  The best model according to the criterion chosen
#' @return $rt:  A two-element vector giving the test statistic and p-value for the Rayleigh Test
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' circ_mle(testdata)

circ_mle = function(data, criterion = "AIC", nchains = 5, BadStart = 10^9, niter = 5000, method = "BFGS"){
    
    # Check all the input for errors
    if (is.null(criterion)) stop("Please provide a model of either AIC, AICc, or BIC")
    if (missing(BadStart)) BadStart = 10^9 else BadStart = BadStart
    if (BadStart < 0) stop("The value for starting parameters outside the preset limits must be >0")
    if (missing(nchains)) nchains = 5 else nchains = nchains
    if (nchains < 1) stop("Must set the number of chains to an integer >=1")
    if (missing(niter)) niter = 5000 else niter = niter
    if (niter < 1000) warning("At least 1000 iterations is recommended but not required. Check ?optim for details.")
    if (missing(method)) method = "BFGS" else method = method
    if (method != "Nelder-Mead" & BadStart == Inf) stop("Except for Nelder-Mead, all other optimization algorithms require finite starting parameters")

# First check that data is in the right format
    data=check_data(data)

# Get Rayleigh test results
    rt = unlist(unlist(rayleigh.test(data))[1:2])

# Get sample size
    ss=length(data)

# Run all 10 models
    m1.out = M1(data)
    m2a.out = M2A(data, BadStart, nchains, method, niter)
    m2b.out = M2B(data, BadStart, nchains, method, niter)
    m2c.out = M2C(data, BadStart, nchains, method, niter)
    m3a.out = M3A(data, BadStart, nchains, method, niter)
    m3b.out = M3B(data, BadStart, nchains, method, niter)
    m4a.out = M4A(data, BadStart, nchains, method, niter)
    m4b.out = M4B(data, BadStart, nchains, method, niter)
    m5a.out = M5A(data, BadStart, nchains, method, niter)
    m5b.out = M5B(data, BadStart, nchains, method, niter)

# Build results table
    results=as.data.frame(matrix(ncol=11,nrow=10))
    rownames(results) = c("M1", "M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B", "M5A", "M5B")
    colnames(results) = c("params", "q1", "k1", "lamda", "q2", "k2", "Likelihood", "Counts.Function", "Counts.Gradient", "Convergence", "Message")
    results[1,] = c("0", "pi", "0", "1", "0", "0", m1.out, "NA", "NA", 0, "NA")
    results[2,] = c(2, unlist(m2a.out)[c(1:2)], "1", "0", "0", unlist(m2a.out)[c(3:7)])
    results[3,] = c(2, unlist(m2b.out)[c(1:2)], "0.5", "0", "0", unlist(m2b.out)[c(3:7)])
    results[4,] = c(3, unlist(m2c.out)[c(1:3)], "0", "0", unlist(m2c.out)[c(4:8)])
    results[5,] = c(2, unlist(m3a.out)[c(1:2)], "0.5", as.numeric(unlist(m3a.out)[1])+pi, unlist(m3a.out)[2], unlist(m3a.out)[c(3:7)])
    results[6,] = c(3, unlist(m3b.out)[c(1:2)], "0.5", as.numeric(unlist(m3b.out)[1])+pi, unlist(m3b.out)[3], unlist(m3b.out)[c(4:8)])
    results[7,] = c(3, unlist(m4a.out)[c(1:3)], as.numeric(unlist(m4a.out)[1])+pi, unlist(m4a.out)[2], unlist(m4a.out)[c(4:8)])
    results[8,] = c(4, unlist(m4b.out)[c(1, 2, 4)], as.numeric(unlist(m4b.out)[1])+pi, unlist(m4b.out)[3], unlist(m4b.out)[c(5:9)])
    results[9,] = c(4, unlist(m5a.out)[c(1, 2, 4, 3, 2)], unlist(m5a.out)[c(5:9)])
    results[10,] = c(5, unlist(m5b.out)[c(1, 2, 5, 3, 4)], unlist(m5b.out)[c(6:10)])

# Add AIC, BIC, other model comparison statistics
    AIC = 2 * (as.numeric(results$Likelihood) + as.numeric(results$params))
    AICc = 2 * (as.numeric(results$Likelihood)) + 2 * (as.numeric(results$params)) * (ss / (ss - (as.numeric(results$params)) - 1))
    BIC = 2 * as.numeric(results$Likelihood) + (log(ss) * as.numeric(results$params))
    results = cbind(results, AIC, AICc, BIC)

# Get deltaAIC
    results = results[order(results$AIC),]
    deltaAIC = results$AIC - results$AIC[1]
    results = cbind(results, deltaAIC)

# Get deltaAICc
    results = results[order(results$AICc),]
    deltaAICc = results$AICc - results$AICc[1]
    results = cbind(results, deltaAICc)

# Get deltaBIC
    results = results[order(results$BIC),]
    deltaBIC = results$BIC - results$BIC[1]
    results = cbind(results, deltaBIC)

# Check convergence diagnostics
    if(any(results$Convergence != 0)) warning("Maximum likelihood search has not converged for one or more models. Check results table")

# Output results as a list
if (criterion == "AIC"){
    results = results[order(results$AIC),]
    bestmodel = rownames(results)[1]
    output = list(results, rt, bestmodel)
    names(output) = c("results", "rt", "bestmodel")
    return(output)
}
if (criterion == "AICc"){
    results = results[order(results$AICc),]
    bestmodel = rownames(results)[1]
    output = list(results, rt, bestmodel)
    names(output) = c("results", "rt", "bestmodel")
    return(output)
}
if (criterion == "BIC"){
    results = results[order(results$BIC),]
    bestmodel = rownames(results)[1]
    output = list(results, rt, bestmodel)
    names(output) = c("results", "rt", "bestmodel")
    return(output)
}
}
