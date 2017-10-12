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
    
    attribute(spectra) <- mat
    
    return(spectra)
    
    
}