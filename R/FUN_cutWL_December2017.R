cutWL <- function(raw,startWL,endWL){

start <- match(startWL, names(raw)) #Find start wavelength
end <- match(endWL, names(raw)) #Find end wavelength

data.rmv.noise <-
        raw[, start:end] #to find correct bands according to col name use match('500',names(x))

Type <- raw$Type #Get response variables back

data.wo.noise <- cbind(Type, data.rmv.noise)

stopifnot(colnames(data.wo.noise[2]) == startWL &
                  colnames(rev(data.wo.noise)[1]) == endWL) # Test if selection correct
return(data.wo.noise)
}
