####
# 0. Script Description
####

# This script contains the full analysis for the first chapter of my PhD: Myrtle Rust (Puccinia
# psidii) Classification on Lemon Myrtle (Backhousia citriodora). All following steps will be 
# commented.

####
# 0. Preparation Steps
####


####
# 1. Loading Packages and Original Data
####

library(dplyr)
library(hsdar)
library(ggplot2)
library(funHDDC)
library(fda) # For Outlier detection
library(fda.usc) # For Outlier detection
library(prospectr) # For spectral binning
library(gdata) 
library(HDclassif)
library(caret)
library(reshape2)
library(gridExtra)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(VSURF)
library(colourpicker)
library(fitdistrplus)
library(png)
library(grid)
library(rlist)
library(qpcR)

data.original <- read.csv('Input_for_C1_AllSpectraABGPlantation_LeafClip.csv', as.is = TRUE, check.names = FALSE)

####
# 2. Cleaning the Original Data (Removing outlier, removing noise and spectral binning)
####

# 2.1 Removing noisy ends from electromagnetic spectrum

names(data.original[,c(3,2153)])  #Check start and end before rmv
data.rmv.noise <- data.original[,153:2153]
names(data.rmv.noise[,c(1,2001)]) #Check bands after rmv

Type <- data.original$Type #Get response variable back 

data.wo.noise <- cbind(Type,data.rmv.noise)

results$data.wo.noise <- data.wo.noise

# 2.2 Remove outlier using fda package

data.wo.noise.mat <- as.matrix(data.wo.noise[,2:2002]) #As matrix to be able to transform the object

#Create fdata object

labnames <- list(main="With", xlab="Wavelength [nm]", ylab="Reflectance [%]") #labnames to have plot information within fdata object
myfdata <- fdata(data.wo.noise.mat, argvals = as.integer(names(data.wo.noise[,2:2002])), names = labnames) #Why as integer??

#Removing outliers

outlier.mat <- outliers.depth.trim(myfdata, dfunc = depth.FM, nb = 10, smo = 0.1, trim = 0.06)

saveRDS(outlier.mat, "outlier.object.rds")
outlier.mat <- readRDS('outlier.object.rds')

#Create fdata object and substract outlier

labnames.out <- list(main="Out", xlab="Wavelength [nm]", ylab="Reflectance [%]")
myfdata.out <- fdata(data.wo.noise.mat,argvals=as.integer(names(data.wo.noise[,2:2002])), names=labnames.out)
myfdata.out <- myfdata[-as.integer(outlier.mat$outliers), ] 

#pdf("Compare including and without Outlier.pdf", width = 16, height = 9)

par(mfrow=c(1,2))

plot(myfdata)
plot(myfdata.out)

#dev.off()

#Create a dataframe and remove outlier for further pre-processing

outlier.vector <- as.numeric(outlier.mat$outliers)
data.wo.noise.for.resampling <- data.wo.noise[-outlier.vector, ]

# Spectral resampling --------------------------------------------------------
Type2 <- data.wo.noise.for.resampling[,1]
data.resamp<-data.wo.noise.for.resampling[,2:2002] #removing Type column
data.bin <- binning(data.resamp, bin.size=10)
data.after.bin <- cbind(Type2, as.data.frame(data.bin))
data.after.bin <- dplyr::rename(data.after.bin, Type = Type2)

# Writing final file without outlier for classification -------------------

write.csv(data.after.bin, 'data.wo.out.binned.cut.csv', row.names=FALSE)

####
# Creating 1st derivative data --------------------------------------------
####

#Load transposed data.wo.out.binned.cut.csv for hsdar package
data4Deri <- read.csv("data.wo.out.binned.cut.transposed.csv", sep=",") #I transposed this file outside of R using xcl. Could not create the correct format using R.

#Create data matrix
myspec <- t(data4Deri[,2:681])
myspec.t <- t(myspec)
spectra <- speclib(myspec.t,data4Deri$Wavelength)

#Splitting rownames stored in spectra@ID to create attributes
attributes <- attribute(spectra)
mySamples<-t(as.data.frame(strsplit(spectra@ID, ".", fixed="T")))
attributes <- data.frame(sampleID=mySamples[,2], Type=mySamples[,1])
attribute(spectra) <- attributes

#Test if spectral object is functional
str(spectra)
He <- subset(spectra, Type == "Healthy")
Tr <- subset(spectra, Type == "Treated")
Un <- subset(spectra, Type == "Untreated")

#Calculate and export 1st derivative data table
d1 <- derivative.speclib(spectra)

He.D <- subset(d1, Type == "Healthy")
Tr.D <- subset(d1, Type == "Treated")
Un.D <- subset(d1, Type == "Untreated")

plot(d1)#to test if functional

data.derivative <- as.data.frame(d1@spectra@spectra_ma)
names(data.derivative) <- d1@wavelength
Type <- attributes$Type
derivative.data <- cbind(Type, data.derivative)
write.csv(derivative.data, '1st.Derivative.Spectra.cleaned.csv', row.names=FALSE)

################################################################################################

####
#
####

data4RF.clean.normal <- read.csv('Input/data.wo.out.binned.cut.csv')
data4RF.clean.1stDeri <- read.csv('Input/1st.Derivative.Spectra.cleaned.csv')

