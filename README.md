# Studies on MyrtleRust-LemonMyrtle-Classification in Northern NSW

This repository contains code needed to reproduce the article:

**Heim RHJ, IJ Wright, HC Chang, A Carnegie, G Pegg, E Lancaster, DS Falster and J Oldeland** Detecting Myrtle Rust (*Austropuccinia psidii*) on lemon myrtle trees using spectral signatures and machine learning. *In review.*

## Data description

Raw data used in the analysis are provided in the file `data/Input_for_C1_AllSpectraABGPlantation_LeafClip.csv`. This file contains the following columns:

- Type: XXX
- Wavelength: XXXX
- 350,351,....,2500: Each column contains values of .... at the specified wavelength in (units).

## Instructions

All analyses were done in `R`. Prior to running the script you may need to install the appropriate packages, all available on the CRAN repository. To install, open R

```
install.packages(c("hsdar", "fda", "fda.usc", "prospectr", "gdata", "reshape2", "cowplot", "VSURF", "colourpicker", "tidyverse", "caret", "randomForest"))
```

To recreate the results, install the appropriate packages and run the commands in the file `analysis.R`. To achieve this [download this repository](https://github.com/reneheim/myrtlerust-lemonmyrtle-classification/archive/master.zip), and then open an R session with working directory set to the root of the project.


Note: This is the first analysis I have made fully reproducible and contains what I would consider now some antiquated coding. But it works!
