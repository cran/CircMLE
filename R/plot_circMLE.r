#' plot_circMLE
#'
#' Plot the observed data and specific model from circ_mle output
#' @param data A vector of class 'circular'
#' @param table A list containing the output from the circ_mle function
#' @param model The name of the model to be plotted. (default = first model of "table")
#' @param bins The number of bins for the circular histogram. (defualt = 18)
#' @param shrink The value by which to shrink the size of the plotted circle.
#' Larger numbers shrink the circle, smaller numbers expand the circle.
#' (default = 1.5)
#' @param col Vector of colors used for plotting, up to four colors can be specified.
#' The order is: histogram color, mean vector color, model density color, predicted
#' mean direction(s) color(s). (default = c("grey", "red", "black", "black"))
#' @param lwd Vector of line weights used for plotting, up to 3 weights can be
#' specified. The order is: mean vector weight, model density weight, predicted
#' mean direction(s) weight(s). (default = c(2, 2, 2))
#' @param lty Vector of line weights used for plotting, up to 3 weights can be
#' specified. The order is: mean vector weight, model density weight, predicted
#' mean direction(s) weight(s). (default = c("solid", "dashed", "dashed"))
#' @keywords plot_circMLE
#' @return A plot is returned to the default image device.
#' @note In the resulting plot, the observed vector length is equal to rho
#' (vector strength). However, the predicted direction(s) from the model are
#' drawn with vector length of 1. Also, if the model "M1" is selected, by
#' definition no predicted direction is plotted.
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3)
#' out = circ_mle(testdata)
#' plot_circMLE(testdata, out)
#' plot_circMLE(testdata, out, model = "M4A")

plot_circMLE = function(data, table, model, bins, shrink, col, lwd, lty){
    if (is.null(data) | missing(data)) stop("Please provide input data vector")
    if (is.null(table) | missing(table)) stop("Please provide the output list from the \"circ_mle\" function")
    
    # If no model specified, choose the first model
    if (missing(model)) model = rownames(table$results)[1]
    if (length(model) != 1) stop("Only 1 model can be specified")
    if (!all(model %in% rownames(table$results))) stop("Model not specified correctly.")
 
    # Check all plotting options
    if (missing(bins)) bins = 18 else bins = bins
    if (length(bins) != 1 | !is.numeric(bins)) stop("Must specify the number of bins as a numeric vector of length 1")
    if (missing(shrink)) shrink = 1.5 else shrink = shrink
    if (length(shrink) != 1 | !is.numeric(shrink)) stop("Must specify \"shrink\" as a numeric vector of length 1")
    
    if (missing(col)) col = c("grey", "red", "black", "black") else col = col
    if (missing(lwd)) lwd = c(2, 2, 2) else lwd = lwd
    if (missing(lty)) lty = c("solid", "dashed", "dashed") else lty = lty
            
    # Coerce into correct format
    data = check_data(data)
    
    # Get all circular attributes from data
    params = circularp(data)
    
    # Make observed histogram (rose plot)
    rose.diag(data, bins = bins, col = col[1], shrink = shrink)
    
    # Plot observed data mean angle
    arrows.circular(mean.circular(data, na.rm = T), y = rho.circular(data, na.rm = T), col = col[2], lwd = lwd[1], lty = lty[1])
    
    # Get the 5 parameters for the model selected
    model.vector = table$results[which(rownames(table$results) == model), 2:6]
    q1 = as.circular(model.vector[1], control.circular = params)
    q2 = as.circular(model.vector[4], control.circular = params)
    k1 = as.numeric(model.vector[2])
    k2 = as.numeric(model.vector[5])
    l = as.numeric(model.vector[3])
    
    # Plot density of best model
     plot.function.circular(function(x) dmixedvonmises(x, q1, q2, k1, k2, l), add = T, lwd = lwd[2], shrink = shrink, lty = lty[2], col = col[3], zero = params$zero, rotation=params$rotation)

    # Plot arrows for the mean angle(s) from the model specified
    if (any(c("M2A", "M2B", "M2C") == model)){
        arrows.circular(q1, col = col[4], lwd = lwd[3], lty = lty[3])
    }
    if (any(c("M3A", "M3B", "M4A", "M4B", "M5A", "M5B") == model)){
        arrows.circular(circular(q1, zero = 0, rotation = "counter"), col = col[4], lwd = lwd[3], lty = lty[3])
        arrows.circular(circular(q2, zero = 0, rotation = "counter"), col = col[4], lwd = lwd[3], lty = lty[3])
    }
}