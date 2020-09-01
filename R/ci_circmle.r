#' ci_circmle
#'
#' Calculate the 95\% confidence interval for estimated model parameters
#' @param circmle A list consisting of the output from function 'circ_mle'
#' @param model character string indicating the model to be used to estimate parameter uncertainty.
#' Must be one of c("M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B",
#' "M5A", "M5B"). Default = the $bestmodel from the circmle object.
#' @keywords circ_mle
#' @return A data frame with a column for the parameter name, the maximum likelihood estimate (MLE),
#' standard error (SE), and 95\% confidence interval (CI) for each estimated model parameter.
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' fit <- circ_mle(testdata)
#' ci_circmle(fit)

ci_circmle = function(circmle, model){
	# Check the the input circmle object is specified properly
	if (is.null(circmle) | missing(circmle)) stop("Please provide the output list from the \"circ_mle\" function")
	if (!identical(names(circmle),c("results", "rt", "bestmodel", "hessians"))) stop("The output list from the \"circ_mle\" function was not specified properly")
    
    # If no model specified, choose the best model
    if (missing(model)) model = circmle$bestmodel
    
    # Check that the model is specified properly
    if (length(model) != 1) stop("Only 1 model can be specified.")
    if (model == "M1") stop("No confidence interval can be calculated for model M1 since it has no parameters.")
    if (!any(c("M2A", "M2B", "M2C", "M3A", "M3B", "M4A", "M4B", "M5A", "M5B") == model)) stop("Model not specified correctly.")
    
    # Get the correct Hessian matrix for the specified model
    H = circmle$hessians[which(names(circmle$hessians) %in% model)][[1]]
    
    # Get the standard error
    se = sqrt(diag(solve(H)))
    
    # Initialize output table
	out = data.frame("Parameter" = vector(), "MLE" = vector(), "SE" = vector(), "CI" = vector())
	r = which(rownames(circmle$results) == model)
	
    # Fill the output table
    if (model == "M2A"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
    }
    if (model == "M2B"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
    }
    if (model == "M2C"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	l = circmle$results[r, 4]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
     	out[3,] <- c("lambda", l, round(se[3], 3), paste0(round(l - 1.96*se[3], 3), " - ", round(l + 1.96*se[3], 3)))
    }
    if (model == "M3A"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
    }
    if (model == "M3B"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	k2 = circmle$results[r, 6]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
     	out[3,] <- c("k2", k2, round(se[3], 3), paste0(round(k2 - 1.96*se[3], 3), " - ", round(k2 + 1.96*se[3], 3)))
    }
    if (model == "M4A"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	l = circmle$results[r, 4]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
     	out[3,] <- c("lambda", l, round(se[3], 3), paste0(round(l - 1.96*se[3], 3), " - ", round(l + 1.96*se[3], 3)))
    }
    if (model == "M4B"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	k2 = circmle$results[r, 6]
    	l = circmle$results[r, 4]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
     	out[3,] <- c("k2", k2, round(se[3], 3), paste0(round(k2 - 1.96*se[3], 3), " - ", round(k2 + 1.96*se[3], 3)))
     	out[4,] <- c("lambda", l, round(se[4], 3), paste0(round(l - 1.96*se[4], 3), " - ", round(l + 1.96*se[4], 3)))
    }
    if (model == "M5A"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	q2 = circmle$results[r, 5]
    	l = circmle$results[r, 4]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
    	out[3,] <- c("q2", q2, round(se[3], 3), paste0(round(q2 - 1.96*se[3], 3), " - ", round(q2 + 1.96*se[3], 3)))
     	out[4,] <- c("lambda", l, round(se[4], 3), paste0(round(l - 1.96*se[4], 3), " - ", round(l + 1.96*se[4], 3)))
    }
    if (model == "M5B"){
    	q1 = circmle$results[r, 2]
    	k1 = circmle$results[r, 3]
    	q2 = circmle$results[r, 5]
    	k2 = circmle$results[r, 6]
    	l = circmle$results[r, 4]
    	out[1,] <- c("q1", q1, round(se[1], 3), paste0(round(q1 - 1.96*se[1], 3), " - ", round(q1 + 1.96*se[1], 3)))
     	out[2,] <- c("k1", k1, round(se[2], 3), paste0(round(k1 - 1.96*se[2], 3), " - ", round(k1 + 1.96*se[2], 3)))
    	out[3,] <- c("q2", q2, round(se[3], 3), paste0(round(q2 - 1.96*se[3], 3), " - ", round(q2 + 1.96*se[3], 3)))
     	out[4,] <- c("k2", k2, round(se[4], 3), paste0(round(k2 - 1.96*se[4], 3), " - ", round(k2 + 1.96*se[4], 3)))
     	out[5,] <- c("lambda", l, round(se[5], 3), paste0(round(l - 1.96*se[5], 3), " - ", round(l + 1.96*se[5], 3)))
    }
    return(out)
}
