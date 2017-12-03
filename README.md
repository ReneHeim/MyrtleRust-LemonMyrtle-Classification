## This repository contains code needed to reproduce the article:

**Heim RHJ, Wright IJ, Chang HC, Carnegie A, Pegg G, Lancaster EK, Falster DS and Oldeland J** Detecting Myrtle Rust (*Austropuccinia psidii*) on lemon myrtle trees using spectral signatures and machine learning. *In review.*

**Supporting Information:**
  
JOURNAL
Date, Volume X, Issue X, pp xxx-xxx

See the article at: LINK
***
    
1. For access to the source code, see *Analysis.R* file or [source code](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/Analysis.R)  
2. For access to the raw data, see *data* folder or [raw data](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv) 
	+ This file contains the following columns:
		- Type: Categorical variables referring to spectral class to be classified.
		- Wavelength: Contains class IDs -> ABG = Australian Botanical Garden Mount Annan, Treated = Fungicide treated plants plantation, Untreated = Untreated plants plantation
		- 350,351,....,2500: Each column contains values of spectral reflectance at the specified wavelength in [%].
3. For access to the used packages and versions, see *packrat* folder or [packages](https://github.com/paternogbc/2015_Rohr_et_al_JAEcol/tree/master/packrat)  
    
***
When using the [raw data](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv), please cite the original publication.

Contact rene.heim@hdr.mq.edu.au for any further information.  

+ This repository follows the principles of reproducible research ([Peng, 2011](http://www.sciencemag.org/content/334/6060/1226)).

# Instructions

All analyses were done in `R`. Prior to running the script you may need to install the appropriate packages, all available on the CRAN repository. To install, open R and execute the following line of code:

```
install.packages(c("hsdar", "fda", "fda.usc", "prospectr", "gdata", "reshape2", "cowplot", "VSURF", "colourpicker", "tidyverse", "caret", "randomForest"))
```

To recreate the results, run the commands in the file `analysis.R`. To achieve this [download this repository](https://github.com/reneheim/myrtlerust-lemonmyrtle-classification/archive/master.zip), and then open an R session with working directory set to the root of the project.


Note: This is the first analysis I have made fully reproducible and contains what I would consider now some antiquated coding. But it works!
