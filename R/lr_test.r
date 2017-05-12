#' Likelihood ratio test for nested models
#'
#' Statistically test nested models for the rejection of the null model
#' in favor of the alternative model.
#' @param data A vector of class 'circular'
#' @param null_model character string indicating the null model to be used.
#' Must be one of c("M1", "M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B",
#' "M5A", "M5B").
#' @param alt_model character string indicating the alternative model to be used.
#' Must be one of c("M1", "M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B",
#' "M5A", "M5B").
#' @keywords likelihood ratio test
#' @keywords chi-square distribution
#' @return A list with elements:
#' @return $null_model: Vector with the model name and number of free parameters
#' @return $alt_model: Vector with the model name and number of free parameters
#' @return $lr: likelihood ratio test statistic 2 * (-log(L1) - -log(L2))
#' @return $df: degrees of freedom (difference in number of parameters between models)
#' @return $p_val: probability of rejecting null model in favor of alternative
#' model due to chance (e.g, if <0.05, data favors alternative model).
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' lr_test(testdata, "M1", "M2A")

lr_test = function(data, null_model, alt_model){
    
    # Data checks
    if (missing(data)) stop("No data specified")
    if (circularp(data)$units == "degrees") {
        warning("Converting from degrees to radians...")
        data <- conversion.circular(data, units = "radians")
    }
    if (circularp(data)$modulo != "2pi") {
        warning("Converting to modulo 2pi...")
        data <- conversion.circular(data, modulo = "2pi")
    }
    if (!is.circular(data)) {
        warning("Data is being coerced into class 'circular'")
        data = check_data(data)
    }
    if (missing(null_model)) stop("The null model must be specified")
    if (missing(alt_model)) stop("The alternative model must be specified")
    if (alt_model == null_model) stop("The null model and alternative model must be different!")

    models = c("M1", "M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B", "M5A", "M5B")
    np = c(0, 2, 2, 3, 2, 3, 3, 4, 4, 5)
    
    if(null_model %in% models == FALSE) stop("The null model was not found or not specified properly. See ?lr_test for details.")
    if(alt_model %in% models == FALSE) stop("The alternative model was not found or not specified properly. See ?lr_test for details.")

    L1 = do.call(null_model, list(data))$lik
    L2 = do.call(alt_model, list(data))$lik
    np1 = np[which(models == null_model)]
    np2 = np[which(models == alt_model)]
    df = np2 - np1
    if (df < 0) warning("Cannot have negative degrees of freedom.  Please use an alternative model with more free parameters than the null model")
    lr = 2 * (L1 - L2)
    p_val = 1 - (stats::pchisq(lr, df))
    result = list(null_model, alt_model, lr, df, p_val)
    return(result)
}
