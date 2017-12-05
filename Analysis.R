####
# Script Description
####

# This script contains the full analysis for the first chapter of my PhD: Myrtle Rust
#(Austrouccinia psidii) Classification on Lemon Myrtle (Backhousia citriodora).
# All following steps will be commented.

# 1. Loading packages and raw data -------------------------------------------

library(hsdar) # Creating spectral libraries
library(fda) # Outlier detection
library(fda.usc) # Outlier detection
library(prospectr) # Spectral resampling
library(gdata) # Drop factor levels
library(reshape2) # Data manipulation
library(cowplot) # Multiplot grids
library(VSURF) # Classification and feature selection
library(colourpicker) # Plotting
library(tidyverse) # Plotting and data manipulation
library(caret) # Classification
library(randomForest) # Classification

source("R/DropCatVar_Type_RF_May2017.R")
source("R/export_VSURF_June2017.R")
source("R/FUN_cutWL_December2017.R")
source("R/prepgg_June2017.R")
source("R/RandomForest_May2017.R") # data split 80:20 as default

dir.create("output", FALSE, FALSE) # Creates folder for generated results
dir.create("data", FALSE, FALSE)
dir.create("docs", FALSE, FALSE)
dir.create("R", FALSE, FALSE)

data.original <-
        read.csv(
                'data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv',
                as.is = TRUE,
                check.names = FALSE
        ) # Data must contain a 1st col named "Type" (Response Var), 2nd col Wavelength (Spectra origin) and all the following labelled according to wavelength (500,501,502...)

# 2. Cleaning raw data (removing outlier, removing noise)  --------


# 2.1 Removing noisy ends from electromagnetic spectrum (500:2500)

data.wo.noise <- cutWL(data.original,'500','2500') # Input raw data and choose start and end wavelength

# 2.2 Remove outlier using fda package

data.wo.noise.mat <-
        as.matrix(data.wo.noise[, 2:2002]) # As matrix to be able to transform the object

labnames <-
        list(main = "With", xlab = "Wavelength [nm]", ylab = "Reflectance [%]") # Labnames to have plot information within fdata object
myfdata <-
        fdata(data.wo.noise.mat,
              argvals = as.integer(names(data.wo.noise[, 2:2002])),
              names = labnames) # fda pkg needs fdata object to run, created here.

outlier.mat <-
        outliers.depth.trim(
                myfdata,
                dfunc = depth.FM,
                nb = 10,
                smo = 0.1,
                trim = 0.06
        ) # Use fdata object and run outlier detection

saveRDS(outlier.mat, "output/outlier.object.rds") # Created object containing outliers
outlier.mat <- readRDS('output/outlier.object.rds')

labnames.out <-
        list(main = "Out", xlab = "Wavelength [nm]", ylab = "Reflectance [%]")

myfdata.out <-
        fdata(data.wo.noise.mat,
              argvals = as.integer(names(data.wo.noise[, 2:2002])),
              names = labnames.out)
myfdata.out <- myfdata[-as.integer(outlier.mat$outliers),] # Created a second fdata after outlier detection for visualization

pdf(
        "output/Compare including and without Outlier.pdf",
        width = 16,
        height = 9
)

par(mfrow = c(1, 2))

plot(myfdata)
plot(myfdata.out)

dev.off() # Visualize outlier detection before and after to check if functional

outlier.vector <- as.numeric(outlier.mat$outliers)
data.wo.noise.for.resampling <- data.wo.noise[-outlier.vector,] # Final df created that no longer contains functional outliers

# 2.4 Spectral resampling

Type <- data.wo.noise.for.resampling[, 1]
data.resamp <-
        data.wo.noise.for.resampling[, 2:2002] #removing "Type" column to resample
data.bin <- binning(data.resamp, bin.size = 10)
data.after.bin <- cbind(Type, as.data.frame(data.bin))

# 2.5 Writing final file without outlier for classification

write.csv(data.after.bin,
          'output/data.wo.out.binned.cut.csv',
          row.names = FALSE)


# 3. Transforming spectra into 1st order derivatives ----------------------

# 3.1 Create data for spectral library

data.after.bin <- read.csv("output/data.wo.out.binned.cut.csv", check.names = FALSE)
Type <- data.after.bin[, 1] 

data.new = t(data.after.bin[, -1])
colnames(data.new) <- data.after.bin[, 1]
rownames(data.new) <-  c()

# 3.2 Create wavelength object for spectral library

Wavelength <- colnames(data.after.bin[, 2:202])
Wavelength <- as.numeric(Wavelength)

# 3.3 Create spectral library "spectra"
# NB -- have investigated warning here and it is fine
spectra <- suppressWarnings(speclib(data.new, Wavelength))

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

plot(d1)#to test if functional

data.derivative <- as.data.frame(d1@spectra@spectra_ma)
names(data.derivative) <- d1@wavelength

derivative.data <- cbind(Type, data.derivative)
write.csv(derivative.data,
          'output/1st.Derivative.Spectra.cleaned.csv',
          row.names = FALSE)

################################################################################################

# 4. Random Forest Classification of Primary and 1st order derivat --------

data4RF.clean.normal <-
        read.csv('output/data.wo.out.binned.cut.csv')
data4RF.clean.1stDeri <-
        read.csv('output/1st.Derivative.Spectra.cleaned.csv')


# 4.1 Apply RF on primary spectra


HvsT <- RFsubset(data4RF.clean.normal, "Untreated")
HvsU <- RFsubset(data4RF.clean.normal, "Treated")
TvsU <- RFsubset(data4RF.clean.normal, "Healthy")

PrimarySpectra_Results <- list()

PrimarySpectra_Results[[1]] <-
        RFapply(data4RF.clean.normal,
                repeats = 1, # Can be set manually
                trees = 1, # Can be set manually
                seq(1, 70, 5)) # Can be set manually
PrimarySpectra_Results[[2]] <-
        RFapply(HvsT,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually
PrimarySpectra_Results[[3]] <-
        RFapply(HvsU,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually
PrimarySpectra_Results[[4]] <-
        RFapply(TvsU,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually

tmp <-
        lapply(PrimarySpectra_Results, function(i) {
                capture.output(print(i) , file = "output/20171205_Results_PrimarySpectra.txt", append =
                                       TRUE)
        })



# 4.2 Apply RF on first order derivative spectra


HvsT.D <- RFsubset(data4RF.clean.1stDeri, "Untreated")
HvsU.D <- RFsubset(data4RF.clean.1stDeri, "Treated")
TvsU.D <- RFsubset(data4RF.clean.1stDeri, "Healthy")

DerivativeSpectra_Results <- list()

DerivativeSpectra_Results[[1]] <-
        RFapply(data4RF.clean.1stDeri,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                seq(1, 70, 5))# Can be set manually
DerivativeSpectra_Results[[2]] <-
        RFapply(HvsT.D,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually
DerivativeSpectra_Results[[3]] <-
        RFapply(HvsU.D,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually
DerivativeSpectra_Results[[4]] <-
        RFapply(TvsU.D,
                repeats = 1,# Can be set manually
                trees = 1,# Can be set manually
                mtry = seq(1, 70, 5))# Can be set manually

tmp <-
        lapply(DerivativeSpectra_Results, function(i) {
                capture.output(print(i) , file = "output/20171205_Results_DerivativeSpectra.txt", append =
                                       TRUE)
        })

# Table 1 A and B has been created based on the output of 4.1 and 4.2
######################################################################################

# 5. VSURF (Variable Selection using Random Forest) -----------------------

####
# The feature selection will be applied on both, primary and first-order derivaive
# spectra, and on all three classes (Healthy, Treated and Infected)
# and also only on classes present on the plantation (Treated and Infected).
####

# 5.1 Take the subsetted data from 4.1 and 4.2 to get four datasets

Full.Deri <- read.csv('output/1st.Derivative.Spectra.cleaned.csv')
Full.Prim <- read.csv('output/data.wo.out.binned.cut.csv')

Planta.Deri <- RFsubset(Full.Deri, "Healthy")
Planta.Prim <- RFsubset(Full.Prim, "Healthy")

# 5.2 Start feature selection

VSURF.Results <- list()

VSURF.Results[[1]] <-
        VSURF(
                Full.Deri[, 2:202],
                Full.Deri[, 1],
                clusterType = "FORK",
                ntree = 2000,
                mtry = 50
        )
VSURF.Results[[2]] <-
        VSURF(
                Full.Prim[, 2:202],
                Full.Prim[, 1],
                clusterType = "FORK",
                ntree = 2000,
                mtry = 50
        )
VSURF.Results[[3]] <-
        VSURF(
                Planta.Deri[, 2:202],
                Planta.Deri[, 1],
                clusterType = "FORK",
                ntree = 2000,
                mtry = 50
        )
VSURF.Results[[4]] <-
        VSURF(
                Planta.Prim[, 2:202],
                Planta.Deri[, 1],
                clusterType = "FORK",
                ntree = 2000,
                mtry = 50
        )

# 5.3 Save results

tmp <-
        lapply(VSURF.Results, function(i) {
                capture.output(print(i$varselect.pred) ,
                               file = "output/20170603_Results_VSURF.txt",
                               append = TRUE)
        })

saveRDS(VSURF.Results, "output/VSURF.Results.rds") 
VSURF.Results <- readRDS('output/VSURF.Results.rds')

# Table 2 A and B has been created based on the output of 5.2

# 6. Prepare Figure 2 -----------------------------------------------------

# 6.1 Manipulate plot data

Full.Deri.gg <- prep.gg(Full.Deri) # Uses sourced function "prep.gg" to transform table lond to wide
Full.Prim.gg <- prep.gg(Full.Prim)
Planta.Deri.gg <- prep.gg(Planta.Deri)
Planta.Prim.gg <- prep.gg(Planta.Prim)


features.FullDeri <-
        export.VSURF(VSURF.Results[[1]]$varselect.pred, Full.Deri[, 2:202]) # Takes VSURF results object extract selected wavebands
features.FullPrim <-
        export.VSURF(VSURF.Results[[2]]$varselect.pred, Full.Deri[, 2:202])
features.PlantaDeri <-
        export.VSURF(VSURF.Results[[3]]$varselect.pred, Full.Deri[, 2:202])
features.PlantaPrim <-
        export.VSURF(VSURF.Results[[4]]$varselect.pred, Full.Deri[, 2:202])

tmp.bands <-
        capture.output(sort(features.FullDeri), file = 'output/bands.txt')
tmp.bands.2 <-
        capture.output(sort(features.FullPrim), file = 'output/bands.txt', append = TRUE)
tmp.bands.3 <-
        capture.output(sort(features.PlantaDeri),
                       file = 'output/bands.txt',
                       append = TRUE)
tmp.bands.4 <-
        capture.output(sort(features.PlantaPrim),
                       file = 'output/bands.txt',
                       append = TRUE)


levels(Full.Deri.gg$Type)[levels(Full.Deri.gg$Type) == "Healthy"] <-
        "Naïve" #Renames factor levels


levels(Full.Prim.gg$Type)[levels(Full.Prim.gg$Type) == "Healthy"] <-
        "Naïve"


levels(Planta.Deri.gg$Type)[levels(Planta.Deri.gg$Type) == "Healthy"] <-
        "Naïve"


levels(Planta.Prim.gg$Type)[levels(Planta.Prim.gg$Type) == "Healthy"] <-
        "Naïve"


# 6.2 Prepare Figure 2 a (p2),b (p4),c (p1) and d (p3)

p1 <-
        ggplot(Full.Deri.gg, aes(Wavelength, Reflectance, colour = Type)) +
        annotate(
                "rect",
                xmin = 500,
                xmax = 570,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'green'
        ) +
        annotate(
                "rect",
                xmin = 570,
                xmax = 590,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'yellow'
        ) +
        annotate(
                "rect",
                xmin = 590,
                xmax = 610,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'orange'
        ) +
        annotate(
                "rect",
                xmin = 610,
                xmax = 700,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("red")
        ) +
        annotate(
                "rect",
                xmin = 700,
                xmax = 1300,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("lightgrey")
        ) +
        annotate(
                "rect",
                xmin = 1300,
                xmax = 2500,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'white'
        ) +
        annotate(
                "text",
                x = 600,
                y = 5,
                label = "VIS",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1000,
                y = 5,
                label = "NIR",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1900,
                y = 5,
                label = "SWIR",
                fontface = "bold",
                size = 5
        ) +
        geom_vline(
                xintercept = features.FullDeri,
                col = "black",
                linetype = "twodash",
                size = 1,
                alpha = .5
        ) +
        geom_line(size = .2) +
        theme_bw() +
        scale_x_continuous(breaks = seq(500, 2500, 150), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Wavelength [nm]", y = "1st Derivative Reflectance") +
        scale_color_manual(values = c("#30D17E", "#030003", "#FF00FF")) +
        theme(legend.text = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(legend.title = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(axis.text = element_text(colour = "black", size = 12)) +
        theme(axis.title.y = element_text(size = 12, angle = 90, face =
                                                  "bold")) +
        theme(axis.title.x = element_text(size = 12, face = "bold")) +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(
                lineheight = .8,
                face = "bold",
                size = 12
        ))


p1

p2 <-
        ggplot(Full.Prim.gg, aes(Wavelength, Reflectance, colour = Type)) +
        annotate(
                "rect",
                xmin = 500,
                xmax = 570,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'green'
        ) +
        annotate(
                "rect",
                xmin = 570,
                xmax = 590,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'yellow'
        ) +
        annotate(
                "rect",
                xmin = 590,
                xmax = 610,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'orange'
        ) +
        annotate(
                "rect",
                xmin = 610,
                xmax = 700,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("red")
        ) +
        annotate(
                "rect",
                xmin = 700,
                xmax = 1300,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("lightgrey")
        ) +
        annotate(
                "rect",
                xmin = 1300,
                xmax = 2500,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'white'
        ) +
        annotate(
                "text",
                x = 600,
                y = 30,
                label = "VIS",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1000,
                y = 30,
                label = "NIR",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1900,
                y = 30,
                label = "SWIR",
                fontface = "bold",
                size = 5
        ) +
        geom_vline(
                xintercept = features.FullPrim,
                col = "black",
                linetype = "twodash",
                size = 1,
                alpha = .5
        ) +
        geom_line(size = .2) +
        theme_bw() +
        scale_x_continuous(breaks = seq(500, 2500, 150), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Wavelength [nm]", y = "Reflectance [%]") +
        scale_color_manual(values = c("#30D17E", "#030003", "#FF00FF")) +
        theme(legend.text = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(legend.title = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(axis.text = element_text(colour = "black", size = 12)) +
        theme(axis.title.y = element_text(size = 12, angle = 90, face =
                                                  "bold")) +
        theme(axis.title.x = element_text(size = 12, face = "bold")) +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position = "none") +
        #ggtitle("Plantation + Botanical Garden") +
        theme(plot.title = element_text(
                lineheight = .8,
                face = "bold",
                size = 12
        ))

p2

p3 <-
        ggplot(Planta.Deri.gg, aes(Wavelength, Reflectance, colour = Type)) +
        annotate(
                "rect",
                xmin = 500,
                xmax = 570,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'green'
        ) +
        annotate(
                "rect",
                xmin = 570,
                xmax = 590,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'yellow'
        ) +
        annotate(
                "rect",
                xmin = 590,
                xmax = 610,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'orange'
        ) +
        annotate(
                "rect",
                xmin = 610,
                xmax = 700,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("red")
        ) +
        annotate(
                "rect",
                xmin = 700,
                xmax = 1300,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("lightgrey")
        ) +
        annotate(
                "rect",
                xmin = 1300,
                xmax = 2500,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'white'
        ) +
        annotate(
                "text",
                x = 600,
                y = 5,
                label = "VIS",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1000,
                y = 5,
                label = "NIR",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1900,
                y = 5,
                label = "SWIR",
                fontface = "bold",
                size = 5
        ) +
        geom_vline(
                xintercept = features.PlantaDeri,
                col = "black",
                linetype = "twodash",
                size = 1,
                alpha = .5
        ) +
        geom_line(size = .2) +
        theme_bw() +
        scale_x_continuous(breaks = seq(500, 2500, 150), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Wavelength [nm]", y = "1st Derivative Reflectance") +
        scale_color_manual(values = c("#030003", "#FF00FF")) +
        theme(legend.text = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(legend.title = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(axis.text = element_text(colour = "black", size = 12)) +
        theme(axis.title.y = element_text(size = 12, angle = 90, face =
                                                  "bold")) +
        theme(axis.title.x = element_text(size = 12, face = "bold")) +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position = "bottom") +
        #ggtitle("The Channon+ABG MA") +
        theme(plot.title = element_text(
                lineheight = .8,
                face = "bold",
                size = 12
        ))

p3

p4 <-
        ggplot(Planta.Prim.gg, aes(Wavelength, Reflectance, colour = Type)) +
        annotate(
                "rect",
                xmin = 500,
                xmax = 570,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'green'
        ) +
        annotate(
                "rect",
                xmin = 570,
                xmax = 590,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'yellow'
        ) +
        annotate(
                "rect",
                xmin = 590,
                xmax = 610,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'orange'
        ) +
        annotate(
                "rect",
                xmin = 610,
                xmax = 700,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("red")
        ) +
        annotate(
                "rect",
                xmin = 700,
                xmax = 1300,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = c("lightgrey")
        ) +
        annotate(
                "rect",
                xmin = 1300,
                xmax = 2500,
                ymin = -Inf,
                ymax = Inf,
                alpha = .2,
                fill = 'white'
        ) +
        annotate(
                "text",
                x = 600,
                y = 30,
                label = "VIS",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1000,
                y = 30,
                label = "NIR",
                fontface = "bold",
                size = 5
        ) +
        annotate(
                "text",
                x = 1900,
                y = 30,
                label = "SWIR",
                fontface = "bold",
                size = 5
        ) +
        geom_vline(
                xintercept = features.PlantaPrim,
                col = "black",
                linetype = "twodash",
                size = 1,
                alpha = .5
        ) +
        geom_line(size = .2) +
        theme_bw() +
        scale_x_continuous(breaks = seq(500, 2500, 150), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Wavelength [nm]", y = "Reflectance [%]") +
        scale_color_manual(values = c("#030003", "#FF00FF")) +
        theme(legend.text = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(legend.title = element_text(
                colour = "black",
                size = 12,
                face = "bold"
        )) +
        theme(axis.text = element_text(colour = "black", size = 12)) +
        theme(axis.title.y = element_text(size = 12, angle = 90, face =
                                                  "bold")) +
        theme(axis.title.x = element_text(size = 12, face = "bold")) +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
        theme(legend.position = "none") +
        #ggtitle("Plantation Only") +
        theme(plot.title = element_text(
                lineheight = .8,
                face = "bold",
                size = 12
        ))


plot.res <-
        plot_grid(
                p2,
                p4,
                p1,
                p3,
                labels = c("a", "b", 'c', 'd'),
                ncol = 2,
                nrow = 2
        )
plot.res
ggsave(
        "output/FeatureSelectionPlot_August2017.pdf",
        plot = plot.res,
        width = 40,
        height = 20,
        units = "cm",
        dpi = 400
)

