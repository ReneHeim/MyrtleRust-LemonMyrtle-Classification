####
# This function is taking a wide dataset from a spectral data assessment and transforms it
# into an object that can be plotted by ggplot2. It needs an X for all the predictor names
# and a response variable "Type"
####

prep.gg <- function(datawide){
        aggdata<-aggregate(.~Type, data=datawide, mean)
        datamelt<-melt(aggdata, id=c("Type"))
        names(datamelt)[2:length(names(datamelt))] <-  c('Wavelength', 'Reflectance')
        datamelt$Wavelength<-gsub("X", "", paste(datamelt$Wavelength))
        datamelt$Wavelength <- as.numeric(datamelt$Wavelength)
        
        datamelt
}


