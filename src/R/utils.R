getParsedArgs = function(){
    args = commandArgs(TRUE)
    ## Default setting when no arguments passed
    if(length(args) < 1) {
         args = c("--help")
    }
 
    ## Help section
    if("--help" %in% args) {
          cat("Help message.")
                  q(save="no")
    }

    ## parse expecte arguments
    parseArgs = function(x) strsplit( sub("^--", "", x), "=")
    argsDF = as.data.frame( do.call(rbind,parseArgs(args)) )
    argsList = as.list(as.character(argsDF$V2))        
    names(argsList) = argsDF$V1

    return(argsList)
}