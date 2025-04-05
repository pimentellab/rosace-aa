# rosaceAA

## Overview

`rosaceAA` is an R package for analyzing growth-based deep mutational scanning screen data. It is built on  `rosace`. 

It offers the following functionalities to convert raw counts to `rosaceAA` scores. Any functionality can be used independently or swapped with custom implementation.
- Raw count filtering (on minimum count and/or max na ratio)
- Missing data imputation (various methods)
- Data normalization (wrt wildtype or total count)
- Score inference from normalized count

The first three functionalities are model-independent and can be used for general purposes.

For details, check out the manuscript on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.09.632281v1).

## Installation

`rosaceAA` uses [cmdstanr](https://mc-stan.org/cmdstanr/) to run inference. Please ensure that `cmdstanr` is properly installed before installing `rosaceAA`. Below is a concise installation command; for complete details, please refer to the official [website](https://mc-stan.org/cmdstanr/).
```{r eval=FALSE}
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 4)  # The number of CPU cores on your device
```

`rosaceAA` also uses the library `impute`. Refer to the Bioconductor [site](https://bioconductor.org/packages/release/bioc/html/impute.html) for installation instructions.

To install `rosaceAA` start [R](https://www.r-project.org) and first install `devtools` by typing
```{r eval=FALSE}
install.packages("devtools")
```

and install `rosaceAA` by typing
```{r eval=FALSE}
devtools::install_github("pimentellab/rosace-aa")
```

For a more complete guide and potential problems during installation, refer to `vignettes/intro_rosaceAA.Rmd`.

## Getting started

```{r eval=FALSE}
library("rosaceAA")
```

We recommend starting with `vignettes/intro_rosaceAA.Rmd`.

## Further help

You may submit an issue [on the GitHub repo](https://github.com/pimentellab/rosace-aa/issues) or email [roserao@ucla.edu](mailto:roserao@ucla.edu?subject=rosaceAA%20Inquiry).

## Citing rosace

Please cite the following publication if you use `rosaceAA`: [under review]. [bioRxiv link](https://www.biorxiv.org/content/10.1101/2025.01.09.632281v1)
