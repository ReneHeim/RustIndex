#' Export relevant results from multiple feature selection runs  
#'
#' @param results Object selected from a list containing multiple VSURF result 
#' outputs
#' @param namesfrom A data frame including specified colums where to get the 
#' names for the return of this function 
#' @return A vector of selected features named according to the waveband that 
#' was selected
#' @examples
#' export.vsurf(feature.set[[i]]$varselect.pred, data[, 2:202])

export_vsurf <- function(results, namesfrom){
        
        hepp <- as.numeric(results)
        gg.names <- names(namesfrom[hepp])
        features <- as.numeric(gsub("X", "", paste(gg.names)))
		features
}

