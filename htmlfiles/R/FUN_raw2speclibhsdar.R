#' Create spectral library using hsdar pkg
#'
#' @param dat A data frame containing spectral data
#' @return A spectral library ready to use in the hsdar pkg
#' @examples
#' DropClass(ori.data, ori.data$Type, "Healthy")
#' Note: It is important the column names are numbers for each waveband only and 
#' we must prevent R checking the names when loading a spectral dataset. R would 
#' add an X in front of each number as this is what R thinks is good name for a 
#' column. Also, the input df must have a first column containing the response
#' variables.

raw2speclib <- function(dat){
    
    require(hsdar)
    
    data.new <- t(dat[,-1])
    colnames(data.new) <- dat[,1]
    rownames(data.new) <-  c()
    
    Wavelength <- colnames(dat[,2:length(dat)])
    Wavelength <- as.numeric(Wavelength)
    
    # NB -- have investigated warning here and it is fine
    spectra <- suppressWarnings(speclib(data.new,Wavelength))
    
    
    mat <- as.matrix(dat[,1])
    colnames(mat) <- colnames(dat[1])
    
    return(spectra)
    
    
}