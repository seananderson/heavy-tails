# Black-swan events in animal populations

This repository contains code for an analysis of heavy-tailed
population dynamics for animal populations. It accompanies the paper:

Anderson, S.C., T.B. Branch, A.B. Cooper, N.K. Dulvy. Black-swan
events in animal populations. In review.

To recreate the analysis, first you will need the following R
packages installed:

```R
install.packages(c("rstan", "dplyr", "plyr", "reshape2", "ggplot2", "gridExtra", 
  "RColorBrewer", "grImport", "TeachingDemos", "metRology", "xtable", "devtools",
  "skewt", "foreach", "Rcpp"))
devtools::install_github("sckott/rphylopic")
```

Then source the following file:

```R
source("analysis/make.R")`
```

Alternatively, in Bash:

```sh
cd analysis
make
```

## Data

The `analysis/gpdd/` folder contains data from the Global Population Dynamics
Database:

NERC Centre for Population Biology, Imperial College (2010). *The Global
Population Dynamics Database Version 2*.
<http://www3.imperial.ac.uk/cpb/databases/gpdd>

The file `analysis/brook-etal.csv` contains a copy of the data from the
supplemental Excel spreadsheet in:

Brook, B.W., Traill, L.W. & Bradshaw, C.J.A. (2006).
Minimum viable population sizes and global extinction risk are unrelated.
*Ecol. Lett.* 9, 375-382. <http://doi.org/10.1111/j.1461-0248.2006.00883.x>

