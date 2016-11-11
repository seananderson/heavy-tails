# Black-swan events in animal populations

This repository contains code for an analysis of heavy-tailed
population dynamics for animal populations. It accompanies the paper:

Anderson, S.C., T.B. Branch, A.B. Cooper, N.K. Dulvy. Black-swan
events in animal populations. In review.

The analysis can be re-created by running the following command in R (version 3.3.1):

```R
source("analysis/make.R")`
```

Alternatively, you can work through the `.R` files in
sequence (`0-...R`, `1-...R`, etc.) in an R console. Files that start with the
same number can be sourced in any order.

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
