####
# This function exports an object from VSURF results (varselect.pred/thresh/inter)
####

export.VSURF <- function(results, namesfrom){
        
        hepp <- as.numeric(results)
        gg.names <- names(namesfrom[hepp])
        features <- as.numeric(gsub("X", "", paste(gg.names)))
		features
}

