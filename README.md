<!-- README.md is generated from README.Rmd. Please edit that file -->
Black-swan events in animal populations
=======================================

This repository contains code for an analysis of heavy-tailed population dynamics for animal populations. It accompanies the paper:

Anderson, S.C., T.B. Branch, A.B. Cooper, N.K. Dulvy. Black-swan events in animal populations. In press at Proceedings of the National Academy of Sciences.

To recreate the analysis, first you will need the following R packages installed:

``` r
install.packages(c("rstan", "dplyr", "plyr", "reshape2", "ggplot2", "gridExtra", 
  "RColorBrewer", "grImport", "TeachingDemos", "metRology", "xtable", "devtools",
  "skewt", "foreach", "Rcpp", "png"))
devtools::install_github("sckott/rphylopic")
```

Theoretically, you should be able to run the entire analysis by sourcing the following R file:

``` r
source("analysis/make.R")
```

However, many of the analyses will take a very long time to run (multiple days), and you might run into minor issues with changes to R packages over time. Therefore, it is probably more useful to open the file `analysis/make.R` and work through the analysis files as they are called. Alternatively, the `.R` files are named in the order they should be run. Files starting with the same number can be run in any order. These numbered files assume that your R working directory is the `analysis` folder. Output from a number of the Stan models is cached in `.rds` files that will be available with a freshly cloned repository. Delete these files to force the models to be fit again.

Data
----

The `analysis/gpdd/` folder contains data from the Global Population Dynamics Database:

NERC Centre for Population Biology, Imperial College (2010). *The Global Population Dynamics Database Version 2*. <http://www3.imperial.ac.uk/cpb/databases/gpdd>

The file `analysis/brook-etal.csv` contains a copy of the data from the supplemental Excel spreadsheet in:

Brook, B.W., Traill, L.W. & Bradshaw, C.J.A. (2006). Minimum viable population sizes and global extinction risk are unrelated. *Ecol. Lett.* 9, 375-382. <http://doi.org/10.1111/j.1461-0248.2006.00883.x>

Package versions
----------------

It's possible that some portions of the analysis were run with earlier versions of R packages, but the versions of R packages installed at the time of publication were:

    #>          package  version
    #> 1       devtools   1.12.0
    #> 2          dplyr    0.5.0
    #> 3        foreach    1.4.3
    #> 4        ggplot2    2.2.1
    #> 5      gridExtra    2.2.1
    #> 6       grImport    0.9-0
    #> 7      metRology 0.9-23-2
    #> 8           plyr    1.8.4
    #> 9            png    0.1-7
    #> 10  RColorBrewer    1.1-2
    #> 11          Rcpp   0.12.9
    #> 12      reshape2    1.4.2
    #> 13         rstan   2.14.1
    #> 14         skewt      0.1
    #> 15 TeachingDemos     2.10
    #> 16        xtable    1.8-2
