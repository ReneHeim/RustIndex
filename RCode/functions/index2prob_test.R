Index2Prob <- function(dat, specdat, index){
    
    SVI <- vegindex(specdat, index)
    
    dat <- cbind(dat, SVI)
    
    SVI.logit <- glm(dat[,1]~dat$SVI, family=binomial(link=probit))
    SVI.prob <- SVI.logit$fitted.values
    
    
    return(SVI.prob)
}


