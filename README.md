# rosaceAA

<p align="left">
  <img src="man/figures/rosace_logo.png" width="150">
</p>


## Overview

__rosaceAA__ is an R package for analyzing growth-based deep mutational scanning screen data. It is built on   __rosace__.

## Installation

__rosaceAA__ uses [cmdstanr](https://mc-stan.org/cmdstanr/) to run inference. Please ensure that __cmdstanr__ is properly installed before installing __rosaceAA__. Below is a concise installation command; for complete details, please refer to the official website. 
```{r eval=FALSE}
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 4)
```

__rosaceAA__ also uses the library __inpute__. Refer to the Bioconductor [site](https://bioconductor.org/packages/release/bioc/html/impute.html) for installation instructions.

To install __rosaceAA__ start [R](https://www.r-project.org) and first install devtools by typing
```{r eval=FALSE}
install.packages("devtools")
```

and install __rosaceAA__ by typing
```{r eval=FALSE}
devtools::install_github("pimentellab/rosace-aa")
```

## Getting started

```{r eval=FALSE}
library("rosaceAA")
```

We recommend starting with the vignette (intro_rosaceAA.Rmd).

## Further help

You may submit a bug report here on GitHub as an issue or you could send an email to roserao@ucla.edu.

## Citing rosace

Please cite the following publication if you use __rosaceAA__: ...
