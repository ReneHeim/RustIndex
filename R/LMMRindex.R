LMMR.index <- function(Band,model){ 
    # Band contains the four important wavelength
    
    res =(exp(coef(model)[1]+(coef(model)[2]*log(Band[1]))+(coef(model)[3]*log(Band[2]))+(coef(model)[4]*log(Band[3]))+(coef(model)[5]*log(Band[4]))))/(1+(exp(coef(model)[1]+(coef(model)[2]*log(Band[1]))+(coef(model)[3]*log(Band[2]))+(coef(model)[4]*log(Band[3]))+(coef(model)[5]*log(Band[4])))))
    
    return(res)
}