
index4col <- function(df.best.bands, model){

  result <- vector('list', length = length(1:nrow(df.best.bands)))
        #Initiate a list for the loop to fill in the results
  
  for( i in 1:nrow(df.best.bands)){
 
    result[[i]] <- LMMR.index(df.best.bands[i,best.bands],model)
  }# Loops through each row of the input df and calc index. Then output into result list

result <- do.call(rbind, result)# takes the list apart to create a df again
names(result)[1] <- 'index_val'

return(result)
}
