####
# Script Description
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

library(hsdar)
library(fda) # For Outlier detection
library(fda.usc) # For Outlier detection
library(prospectr) # For spectral binning
library(gdata) #To drop levels of a df
library(reshape2)
library(cowplot)
library(VSURF)
library(colourpicker)
library(tidyverse)
library(caret)
library(randomForest)
dir.create("output", FALSE, FALSE)

data.original <- read.csv('data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv', as.is = TRUE, check.names = FALSE)

####
# 2. Cleaning the Original Data (Removing outlier, removing noise and spectral binning)
####

# 2.1 Removing noisy ends from electromagnetic spectrum (500:2500)
start <- match('500', names(data.original)) #Find start wavelength
end <- match('2500', names(data.original)) #Find end wavelength

data.rmv.noise <- data.original[,start:end] #to find correct bands according to col name use match('500',names(x))

Type <- data.original$Type #Get response variables back 

data.wo.noise <- cbind(Type,data.rmv.noise)

stopifnot(colnames(data.wo.noise[2])=='500' & colnames(rev(data.wo.noise)[1])=='2500')

# 2.2 Remove outlier using fda package

data.wo.noise.mat <- as.matrix(data.wo.noise[,2:2002]) #As matrix to be able to transform the object

#Create fdata object

labnames <- list(main="With", xlab="Wavelength [nm]", ylab="Reflectance [%]") #labnames to have plot information within fdata object
myfdata <- fdata(data.wo.noise.mat, argvals = as.integer(names(data.wo.noise[,2:2002])), names = labnames)

#Removing outliers

outlier.mat <- outliers.depth.trim(myfdata, dfunc = depth.FM, nb = 10, smo = 0.1, trim = 0.06)

saveRDS(outlier.mat, "output/outlier.object.rds")
outlier.mat <- readRDS('output/outlier.object.rds')

#Create fdata object and remove outlier

labnames.out <- list(main="Out", xlab="Wavelength [nm]", ylab="Reflectance [%]")
myfdata.out <- fdata(data.wo.noise.mat,argvals=as.integer(names(data.wo.noise[,2:2002])), names=labnames.out)
myfdata.out <- myfdata[-as.integer(outlier.mat$outliers), ] 

pdf("output/Compare including and without Outlier.pdf", width = 16, height = 9)

par(mfrow=c(1,2))

plot(myfdata)
plot(myfdata.out)

dev.off()

#Create a dataframe and remove outlier for further pre-processing

outlier.vector <- as.numeric(outlier.mat$outliers)
data.wo.noise.for.resampling <- data.wo.noise[-outlier.vector, ]

# 2.4 Spectral resampling

Type <- data.wo.noise.for.resampling[,1]
data.resamp<-data.wo.noise.for.resampling[,2:2002] #removing Type column
data.bin <- binning(data.resamp, bin.size=10)
data.after.bin <- cbind(Type, as.data.frame(data.bin))

# 2.5 Writing final file without outlier for classification

write.csv(data.after.bin, 'output/data.wo.out.binned.cut.csv', row.names=FALSE)

####
# 3. Transforming spectra into 1st order derivatives and store dataframe
####

# 3.1 Create data for spectral library

data.new = t(data.after.bin[,-1])
colnames(data.new) <- data.after.bin[,1]
rownames(data.new) <-  c()

# 3.2 Create wavelength object for spectral library

Wavelength <- colnames(data.after.bin[,2:202])
Wavelength <- as.numeric(Wavelength)

# 3.3 Create spectral library "spectra"
# NB -- have investigated warning here and it is fine
spectra <- suppressWarnings(speclib(data.new,Wavelength))

# 3.4 Create attribute object for spectra
mat <- as.matrix(data.after.bin$Type)
colnames(mat) <- c("Type")

attribute(spectra) <- mat

# 3.5 Test if spectral object is functional

str(spectra)
He <- subset(spectra, Type == "Healthy")
Tr <- subset(spectra, Type == "Treated")
Un <- subset(spectra, Type == "Untreated")

# 3.6 Calculate and export 1st derivative data table
d1 <- derivative.speclib(spectra)

He.D <- subset(d1, Type == "Healthy")
Tr.D <- subset(d1, Type == "Treated")
Un.D <- subset(d1, Type == "Untreated")

plot(d1 )#to test if functional

data.derivative <- as.data.frame(d1@spectra@spectra_ma)
names(data.derivative) <- d1@wavelength

Type <- attribute(spectra)$Type
derivative.data <- cbind(Type, data.derivative)
write.csv(derivative.data, 'output/1st.Derivative.Spectra.cleaned.csv', row.names=FALSE)

################################################################################################

####
# 4. Random Forest Classification of Primary and 1st order derivative spectra
####


data4RF.clean.normal <- read.csv('output/data.wo.out.binned.cut.csv')
data4RF.clean.1stDeri <- read.csv('output/1st.Derivative.Spectra.cleaned.csv')

source('R/RandomForest_May2017.r') #data split 80:20 as default
source('R/DropCatVar_Type_RF_May2017.R')


# 4.1 Apply RF on primary spectra


HvsT <- RFsubset(data4RF.clean.normal,"Untreated")
HvsU <- RFsubset(data4RF.clean.normal,"Treated")
TvsU <- RFsubset(data4RF.clean.normal,"Healthy")  

PrimarySpectra_Results <- list()

PrimarySpectra_Results[[1]] <- RFapply(data4RF.clean.normal,repeats=1,trees=1,seq(1,70,5))
PrimarySpectra_Results[[2]] <- RFapply(HvsT,repeats=1,trees=1,mtry=seq(1,70,5))
PrimarySpectra_Results[[3]] <- RFapply(HvsU,repeats=1,trees=1,mtry=seq(1,70,5))
PrimarySpectra_Results[[4]] <- RFapply(TvsU,repeats=1,trees=1,mtry=seq(1,70,5))

tmp <- lapply(PrimarySpectra_Results, function(i){  capture.output( print(i) , file="output/20170528_Results_PrimarySpectra.txt", append=TRUE)})



# 4.2 Apply RF on first order derivative spectra


HvsT.D <- RFsubset(data4RF.clean.1stDeri,"Untreated")
HvsU.D <- RFsubset(data4RF.clean.1stDeri,"Treated")
TvsU.D <- RFsubset(data4RF.clean.1stDeri,"Healthy")  

DerivativeSpectra_Results <- list()

DerivativeSpectra_Results[[1]] <- RFapply(data4RF.clean.1stDeri,repeats=1,trees=1,seq(1,70,5))
DerivativeSpectra_Results[[2]] <- RFapply(HvsT.D,repeats=1,trees=1,mtry=seq(1,70,5))
DerivativeSpectra_Results[[3]] <- RFapply(HvsU.D,repeats=1,trees=1,mtry=seq(1,70,5))
DerivativeSpectra_Results[[4]] <- RFapply(TvsU.D,repeats=1,trees=1,mtry=seq(1,70,5))

tmp <- lapply(DerivativeSpectra_Results, function(i){  capture.output( print(i) , file="output/20170528_Results_DerivativeSpectra.txt", append=TRUE)})

######################################################################################

####
# 5. VSURF (Variable Selection using Random Forest)
#
# The feature selection will be applied on both, primary and first-order derivaive
# spectra, and on all three classes (Healthy, Treated and Infected)
# and also only on classes present on the plantation (Treated and Infected).
####

source('R/DropCatVar_Type_RF_May2017.R')

# 5.1 Take the subsetted data from 4.1 and 4.2 to get four datasets

Full.Deri <- read.csv('output/1st.Derivative.Spectra.cleaned.csv')
Full.Prim <- read.csv('output/data.wo.out.binned.cut.csv')
        
Planta.Deri <- RFsubset(Full.Deri,"Healthy")
Planta.Prim <- RFsubset(Full.Prim,"Healthy")        

# 5.2 Start feature selection
VSURF.Results <- list()

VSURF.Results[[1]] <- VSURF(Full.Deri[,2:202], Full.Deri[,1], clusterType = "FORK", ntree = 2000, mtry = 50)
VSURF.Results[[2]] <- VSURF(Full.Prim[,2:202], Full.Prim[,1], clusterType = "FORK", ntree = 2000, mtry = 50)
VSURF.Results[[3]] <- VSURF(Planta.Deri[,2:202], Planta.Deri[,1], clusterType = "FORK", ntree = 2000, mtry = 50)
VSURF.Results[[4]] <- VSURF(Planta.Prim[,2:202], Planta.Deri[,1], clusterType = "FORK", ntree = 2000, mtry = 50)

# 5.3 Save results 
tmp <- lapply(VSURF.Results, function(i){  capture.output( print(i$varselect.pred) , file="output/20170603_Results_VSURF.txt", append=TRUE)})

saveRDS(VSURF.Results, "output/VSURF.Results.2.rds")
VSURF.Results <- readRDS('output/VSURF.Results.2.rds')

# 5.4 Prepare ggplot (Figure X)

source('R/prepgg_June2017.R')

Full.Deri.gg <- prep.gg(Full.Deri)# needs wide dataset 
Full.Prim.gg <- prep.gg(Full.Prim)
Planta.Deri.gg <- prep.gg(Planta.Deri)
Planta.Prim.gg <- prep.gg(Planta.Prim)

source('R/export_VSURF_June2017.R')

features.FullDeri <- export.VSURF(VSURF.Results[[1]]$varselect.pred, Full.Deri[,2:202])
features.FullPrim <- export.VSURF(VSURF.Results[[2]]$varselect.pred, Full.Deri[,2:202])
features.PlantaDeri <- export.VSURF(VSURF.Results[[3]]$varselect.pred, Full.Deri[,2:202])
features.PlantaPrim <- export.VSURF(VSURF.Results[[4]]$varselect.pred, Full.Deri[,2:202])


levels(Full.Deri.gg$Type)[levels(Full.Deri.gg$Type)=="Treated"] <- "Uninfected"
levels(Full.Deri.gg$Type)[levels(Full.Deri.gg$Type)=="Untreated"] <- "Infected"

levels(Full.Prim.gg$Type)[levels(Full.Prim.gg$Type)=="Treated"] <- "Uninfected"
levels(Full.Prim.gg$Type)[levels(Full.Prim.gg$Type)=="Unreated"] <- "Infected"

levels(Planta.Deri.gg$Type)[levels(Planta.Deri.gg$Type)=="Treated"] <- "Uninfected"
levels(Planta.Deri.gg$Type)[levels(Planta.Deri.gg$Type)=="Untreated"] <- "Infected"

levels(Planta.Prim.gg$Type)[levels(Planta.Prim.gg$Type)=="Treated"] <- "Uninfected"
levels(Planta.Prim.gg$Type)[levels(Planta.Prim.gg$Type)=="Untreated"] <- "Infected"

# 5.5 Prepare vlines from features

p1 <-ggplot(Full.Deri.gg, aes(Wavelength, Reflectance, colour = Type))+
        annotate("rect", xmin = 500, xmax = 570, ymin = -Inf, ymax = Inf, alpha = .2, fill='green')+
        annotate("rect", xmin = 570, xmax = 590, ymin = -Inf, ymax = Inf, alpha = .2, fill='yellow')+
        annotate("rect", xmin = 590, xmax = 610, ymin = -Inf, ymax = Inf, alpha = .2, fill='orange')+
        annotate("rect", xmin = 610, xmax = 700, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("red"))+
        annotate("rect", xmin = 700, xmax = 1300, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("lightgrey"))+
        annotate("rect", xmin = 1300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = .2, fill='white')+
        annotate("text", x= 600, y= 5, label="VIS", fontface="bold", size=5)+
        annotate("text", x= 1000, y= 5, label="NIR", fontface="bold", size=5)+
        annotate("text", x= 1900, y= 5, label="SWIR", fontface="bold", size=5)+
        geom_vline(xintercept = features.FullDeri, col = "black", linetype = "twodash", size = 1, alpha=.5)+
        geom_line(size=.2)+
        theme_bw()+
        scale_x_continuous(breaks=seq(500,2500,150),expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0))+
        labs(x="Wavelength [nm]", y="1st Derivative Reflectance")+
        scale_color_manual(values=c("#30D17E", "#030003", "#FF00FF"))+
        theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
        theme(legend.title = element_text(colour="black", size=12, face="bold"))+
        theme(axis.text = element_text(colour = "black",size=12))+  
        theme(axis.title.y = element_text(size = 12, angle = 90,face="bold"))+
        theme(axis.title.x = element_text(size = 12,face="bold"))+
        theme(legend.position="bottom")+
        theme(plot.title = element_text(lineheight=.8, face="bold", size = 12)) 
        

p1

p2 <- ggplot(Full.Prim.gg, aes(Wavelength, Reflectance, colour = Type))+
        annotate("rect", xmin = 500, xmax = 570, ymin = -Inf, ymax = Inf, alpha = .2, fill='green')+
        annotate("rect", xmin = 570, xmax = 590, ymin = -Inf, ymax = Inf, alpha = .2, fill='yellow')+
        annotate("rect", xmin = 590, xmax = 610, ymin = -Inf, ymax = Inf, alpha = .2, fill='orange')+
        annotate("rect", xmin = 610, xmax = 700, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("red"))+
        annotate("rect", xmin = 700, xmax = 1300, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("lightgrey"))+
        annotate("rect", xmin = 1300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = .2, fill='white')+
        annotate("text", x= 600, y= 30, label="VIS", fontface="bold", size=5)+
        annotate("text", x= 1000, y= 30, label="NIR", fontface="bold", size=5)+
        annotate("text", x= 1900, y= 30, label="SWIR", fontface="bold", size=5)+
        geom_vline(xintercept = features.FullPrim, col = "black", linetype = "twodash", size = 1, alpha=.5)+
        geom_line(size=.2)+
        theme_bw()+
        scale_x_continuous(breaks=seq(500,2500,150),expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0))+
        labs(x="Wavelength [nm]", y="Reflectance [%]")+
        scale_color_manual(values=c("#30D17E", "#030003", "#FF00FF"))+
        theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
        theme(legend.title = element_text(colour="black", size=12, face="bold"))+
        theme(axis.text = element_text(colour = "black",size=12))+  
        theme(axis.title.y = element_text(size = 12, angle = 90,face="bold"))+
        theme(axis.title.x = element_text(size = 12,face="bold"))+
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position="none")+
        #ggtitle("Plantation + Botanical Garden") + 
        theme(plot.title = element_text(lineheight=.8, face="bold", size = 12))+ coord_cartesian(xlim = c(650, 750))

p2

p3 <- ggplot(Planta.Deri.gg, aes(Wavelength, Reflectance, colour = Type))+
        annotate("rect", xmin = 500, xmax = 570, ymin = -Inf, ymax = Inf, alpha = .2, fill='green')+
        annotate("rect", xmin = 570, xmax = 590, ymin = -Inf, ymax = Inf, alpha = .2, fill='yellow')+
        annotate("rect", xmin = 590, xmax = 610, ymin = -Inf, ymax = Inf, alpha = .2, fill='orange')+
        annotate("rect", xmin = 610, xmax = 700, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("red"))+
        annotate("rect", xmin = 700, xmax = 1300, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("lightgrey"))+
        annotate("rect", xmin = 1300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = .2, fill='white')+
        annotate("text", x= 600, y= 5, label="VIS", fontface="bold", size=5)+
        annotate("text", x= 1000, y= 5, label="NIR", fontface="bold", size=5)+
        annotate("text", x= 1900, y= 5, label="SWIR", fontface="bold", size=5)+
        geom_vline(xintercept = features.PlantaDeri, col = "black", linetype = "twodash", size = 1, alpha=.5)+
        geom_line(size=.2)+
        theme_bw()+
        scale_x_continuous(breaks=seq(500,2500,150),expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0))+
        labs(x="Wavelength [nm]", y="1st Derivative Reflectance")+
        scale_color_manual(values=c("#030003", "#FF00FF"))+
        theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
        theme(legend.title = element_text(colour="black", size=12, face="bold"))+
        theme(axis.text = element_text(colour = "black",size=12))+  
        theme(axis.title.y = element_text(size = 12, angle = 90,face="bold"))+
        theme(axis.title.x = element_text(size = 12,face="bold"))+
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position="bottom")+
        #ggtitle("The Channon+ABG MA") + 
        theme(plot.title = element_text(lineheight=.8, face="bold", size = 12))

p3

p4 <- ggplot(Planta.Prim.gg, aes(Wavelength, Reflectance, colour = Type))+
        annotate("rect", xmin = 500, xmax = 570, ymin = -Inf, ymax = Inf, alpha = .2, fill='green')+
        annotate("rect", xmin = 570, xmax = 590, ymin = -Inf, ymax = Inf, alpha = .2, fill='yellow')+
        annotate("rect", xmin = 590, xmax = 610, ymin = -Inf, ymax = Inf, alpha = .2, fill='orange')+
        annotate("rect", xmin = 610, xmax = 700, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("red"))+
        annotate("rect", xmin = 700, xmax = 1300, ymin = -Inf, ymax = Inf, alpha = .2, fill=c("lightgrey"))+
        annotate("rect", xmin = 1300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = .2, fill='white')+
        annotate("text", x= 600, y= 30, label="VIS", fontface="bold", size=5)+
        annotate("text", x= 1000, y= 30, label="NIR", fontface="bold", size=5)+
        annotate("text", x= 1900, y= 30, label="SWIR", fontface="bold", size=5)+
        geom_vline(xintercept = features.PlantaPrim, col = "black", linetype = "twodash", size = 1, alpha=.5)+
        geom_line(size=.2)+
        theme_bw()+
        scale_x_continuous(breaks=seq(500,2500,150),expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0))+
        labs(x="Wavelength [nm]", y="Reflectance [%]")+
        scale_color_manual(values=c("#030003", "#FF00FF"))+
        theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
        theme(legend.title = element_text(colour="black", size=12, face="bold"))+
        theme(axis.text = element_text(colour = "black",size=12))+  
        theme(axis.title.y = element_text(size = 12, angle = 90,face="bold"))+
        theme(axis.title.x = element_text(size = 12,face="bold"))+
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position="none")+
        #ggtitle("Plantation Only") + 
        theme(plot.title = element_text(lineheight=.8, face="bold", size = 12))


plot.res <- plot_grid(p2, p4, p1, p3, labels=c("A", "B", 'C', 'D'), ncol = 2, nrow = 2)
plot.res
ggsave("output/FeatureSelectionPlot_July2017_Vers2.pdf", plot=plot.res, width = 40, height = 20, units = "cm", dpi = 400)

