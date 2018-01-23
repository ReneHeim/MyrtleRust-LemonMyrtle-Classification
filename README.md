## This repository contains code needed to reproduce the article:

Heim, R. H. J., Wright, I. J., Chang, H.-C., Carnegie, A. J., Pegg, G. S., Lancaster, E. K., Falster, D. S. and Oldeland, J. (2018) 
‘Detecting myrtle rust ( Austropuccinia psidii ) on lemon myrtle trees using spectral signatures and machine learning’, Plant Pathology, 12(10), pp. 3218–3221. doi: 10.1111/ppa.12830.

Find the article [HERE](http://onlinelibrary.wiley.com/doi/10.1111/ppa.12830/full)
***
    
1. For access to the source code, see `Analysis.R` file or [source code](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/Analysis.R)  
2. For access to the raw data, see *data* folder or [raw data](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv) 
	+ This file contains the following columns:
		- Type: Categorical variables referring to spectral class to be classified.
		- Wavelength: Contains class IDs -> ABG = Australian Botanical Garden Mount Annan, Treated = Fungicide treated plants plantation, Untreated = Untreated plants plantation
		- 350,351,....,2500: Each column contains values of spectral reflectance at the specified wavelength in [%].
    
***
When using the [raw data](https://github.com/ReneHeim/MyrtleRust-LemonMyrtle-Classification/blob/master/data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv), please cite the original publication.

Contact rene.heim@hdr.mq.edu.au for any further information.  

+ This repository follows the principles of reproducible research ([Peng, 2011](http://www.sciencemag.org/content/334/6060/1226)).

# Instructions

All analyses were done in `R`. Prior to running the script you may need to install the appropriate packages, all available on the CRAN repository. To install, open `R` and execute the following line of code:

```
install.packages(c("hsdar", "fda", "fda.usc", "prospectr", "gdata", "reshape2", "cowplot", "VSURF", "colourpicker", "tidyverse", "caret", "randomForest", "e1071"))
```

To recreate the results, run the commands in the file `Analysis.R`. To achieve this [download this repository](https://github.com/reneheim/myrtlerust-lemonmyrtle-classification/archive/master.zip), and then open an `R` session with working directory set to the root of the project.


Note: 

+ This is the first analysis I have made fully reproducible and contains what I would consider now some antiquated coding. But it works!
+ Table 1 A and B has been created based on the analysis output of 4.1 and 4.2
+ Table 2 A and B has been created based on the analysis output of 5.2
+ All classification (4.1, 4.2) and feature selection (5.2) steps can be set to minimum repetition values to run the code swiftly for testing. To reproduce the paper's results the mentioned steps must be set to the values given in the article. (Ten-fold cross-validation has been repeated 100 times,  random forest was run using 2000 trees and sequential mtry up until 70 included predictor variables starting at 1 predictor variable and increasing in steps of 5. VSURF was run using 2000 trees and mtry set to 50.

