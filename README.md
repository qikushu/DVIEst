# DVIEst

A package to conduct development index (DVI) modeling, a.k.a developmental rate 
(DVR) model. DVI models are phenological models that simulate reproductive 
initiation and flowering of crops, especially rice. Parameters of a DVI model
can be estimated based on given climate information, sowing and heading dates 
of a variety. Prediction of heading dates is then achieved using the DVI model 
with the estimated parameters and sowing dates.

## Introduction


## Prerequisit
DVIEst requires the following packages.
```
install.packages("GA")
install.packages("ggplot2")
```
Building the package vignette may also requires the following packages.
```
install.packages("knitr")
install.packages("rmarkdown")
```

## Installation
You can install `DVIEst` from the GitHub repository.
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("tomoyukif/DVIEst", build_vignettes = TRUE)
```

For more information see vignette or run the following code on a R console.
```
browseVignettes(package = "DVIEst")
```

## Citations
In preparation for publication...
