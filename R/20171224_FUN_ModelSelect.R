GLM.featureselect <- function(dat, pred){
    
    require(glmulti)
    
    #pred must be numeric variables
    #dat must have response in first col
    
    res <- glmulti(y=names(dat)[1],xr=paste0('X',pred),data.log,maxsize=4,level = 1, family=binomial)
    
    model <- glm(as.formula(summary(res)$bestmodel),dat,family=binomial)
    
    return(model)
} 
