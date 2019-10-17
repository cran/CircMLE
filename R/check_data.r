#' Data Checking Function
#'
#' Make sure data is in the right format.  Datasets are coerced into class 'circular'
#' of type = angles, units = radians, and modulo = 2pi.  It is recommended to set
#' these atributes ahead of time.
#' @param data A vector, class 'circular' is recommended but not required
#' @keywords data checking
#' @import circular
#' @export
#' @examples
#' testdata = circular::rvonmises(100, mu = circular::circular(pi), kappa = 3,
#'    control.circular = list(units = "degrees"))
#' check_data(testdata)

check_data = function(data){
    data <- data[!is.na(data)]
    if (length(data) == 0) {
        stop("Dataset is empty after removing 'NA' values")
        return(NULL)
    }
    
    if (!is.circular(data)) {
    	message("Converting data to an object of class \"circular\", assumed to be in units of radians.\n")
        data = circular(data, type = "angles", units = "radians", modulo = "2pi")
    }
    
    if (circularp(data)$units=="degrees") {
        warning("Converting from degrees to radians...")
        data <- conversion.circular(data, units="radians")
    }
    if (circularp(data)$modulo != "2pi") {
        warning("Converting to modulo 2pi...")
        data <- conversion.circular(data, modulo="2pi")
    }
    return(data)
}
