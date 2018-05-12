####
# Introduction
#
# This function is taking a df containing a categorical response variable called "Type". The 
# second argument in the function is the response variable of "Type" to drop.
####

DropClass <-  function(data, column, classtodrop){
        
        require('gdata')
        
        data2 <- subset(data, column != classtodrop)
        
        #Check if levels are structured correctly and adapt if necessary (gdata) #
        
        data2 <- drop.levels(data2)
        
        return(data2)
}

