#' Model M1 MLE function
#'
#' Run Maximum likelihood estimation for model M1.
#' @param data A vector of class 'circular'
#' @keywords M1
#' @return A list with the elements:
#' @return $lik: The negative log likelihood of data for model M1
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' M1(testdata)

M1 = function(data){
    
    # Since there are no free parameters, just calculate the negative log likelihood
    P = circularp(data)
    R = -sum(log(dmixedvonmises(data, mu1 = as.circular(pi, control.circular = P), mu2 = as.circular(0, control.circular = P), kappa1 = 0, kappa2 = 0, prop = 1)))
    R = list(R)
    names(R) = "lik"
    return(R)
}
