#' Model M3B MLE function
#'
#' Run Maximum likelihood estimation for model M3B.
#' @param data A vector of class 'circular'
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
#' @keywords M3B
#' @return A list with elements (same as for function optim()):
#' @return $par:  Vector with the optimized mean angle (mu1),
#' concentration parameter (kappa1), and concentration parameter (kappa2).
#' @return $lik:  The negative log likelihood
#' @return $counts:  A two-element integer vector giving the number of calls to
#' ‘fn’ and ‘gr’ respectively. See ?optim() for details.
#' @return $convergence:  An integer code. ‘0’ indicates successful completion
#' (which is always the case for ‘"SANN"’ and ‘"Brent"’). Possible error codes are:
#' @return ‘1’ indicates that the iteration limit ‘maxit’ had been reached.
#' @return ‘10’ indicates degeneracy of the Nelder-Mead simplex.
#' @return ‘51’ indicates a warning from the ‘"L-BFGS-B"’ method; see
#' component ‘message’ for further details.
#' @return ‘52’ indicates an error from the ‘"L-BFGS-B"’ method; see
#' component ‘message’ for further details.
#' @return $message:  A character string giving any additional information returned
#' by the optimizer, or ‘NULL’.
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' M3B(testdata)

M3B = function(data, BadStart, nchains, method, niter){
    
    if (missing(BadStart)) BadStart = 10^9 else BadStart = BadStart
    if (BadStart < 0) stop("The value for starting parameters outside the preset limits must be >0")
    if (missing(nchains)) nchains = 5 else nchains = nchains
    if (nchains < 1) stop("Must set the number of chains to an integer >=1")
    if (missing(niter)) niter = 5000 else niter = niter
    if (niter < 1000) warning("At least 1000 iterations is recommended but not required. Check ?optim for details.")
    if (missing(method)) method = "BFGS" else method = method
    if (method != "Nelder-Mead" & BadStart == Inf) stop("Except for Nelder-Mead, all other optimization algorithms require finite starting parameters")

    m3b = function(params){
        if (params[1] < 0 | params[1] > 2*pi | params[2] <= 0 | params[2] > 227 | params[3] <= 0 | params[3] > 227){
            R = BadStart
            return(R)
        }
        else {
        	P = circularp(data)
            R = dmixedvonmises(data, mu1 = as.circular(params[1], control.circular = P), mu2 = as.circular(params[1]+pi, control.circular = P), kappa1 = params[2], kappa2 = params[3], prop = 0.5)
            R = -sum(log(R))
            return(R)
        }
    }
    
    # Randomize starting parameters for the optimization
    q1 = as.numeric(rcircularuniform(nchains, control.circular = list(modulo = "2pi")))
    k1 = as.numeric(sample(1:5, nchains, replace = T))
    k2 = as.numeric(sample(1:5, nchains, replace = T))
    
    # Run optimization
    m3b.out = list()
    for (i in 1:nchains){
        chain.out = suppressWarnings(stats::optim(c(q1[i], k1[i], k2[i]), fn = m3b, method = method, control = list(maxit = niter), hessian = T))
        names(chain.out)[2] = "lik"
        m3b.out[[i]] = chain.out
    }
    min = which.min(sapply(m3b.out,function(x) x[2]))
    return(m3b.out[[min]])
}
