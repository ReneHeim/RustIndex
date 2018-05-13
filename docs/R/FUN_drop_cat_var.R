#' Dropping factor levels
#'
#' @param data A data frame
#' @param column Column containing factors
#' @param classtodrop Factor level to drop
#' @return A data frame without the specified factor and according levels
#' @examples
#' DropClass(ori.data, ori.data$Type, "Healthy")
#' DropClass(iris, iris$Species, "setosa")

drop_class <-  function(data, column, classtodrop){
        
        require('gdata')
        
        data2 <- subset(data, column != classtodrop)
        
        data2 <- drop.levels(data2)
        
        return(data2)
}

