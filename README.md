# Black swans in ecological time series

This repository holds code for an in-progress analysis of heavy-tailed
population dynamics. This analysis is for a PhD thesis chapter and paper by
Sean Anderson, Trevor Branch, Andy Cooper, and Nick Dulvy.

The analysis can be run by sourcing the `.R` files in sequence
(`1-...R`, `2-...R`, etc.). Files that start with the same number can be
sourced in any order.

The `gpdd` folder contains data from the Global Population Dynamics Database:

NERC Centre for Population Biology, Imperial College (2010). *The Global
Population Dynamics Database Version 2*.
<http://www3.imperial.ac.uk/cpb/databases/gpdd>

The file `brook-etal.csv` contains a copy of the data from the supplemental Excel spreadsheet in:

Brook, B.W., Traill, L.W. & Bradshaw, C.J.A. (2006).
Minimum viable population sizes and global extinction risk are unrelated.
*Ecol. Lett.* 9, 375-382. <http://doi.org/10.1111/j.1461-0248.2006.00883.x>

To run the analysis you will need the following R packages (package version numbers are included in parentheses):

### Statistics
- rstan (2.4.0) (with Stan 2.4)
- metRology (0.9-17)

### Data manipulation
- dplyr (0.3.0.9000)
- plyr (1.8.1)
- reshape2 (1.4)

### Plotting
- ggplot2 (1.0.0.99)
- gridExtra (0.9.1)
- RColorBrewer (1.0-5)
- grImport (0.9-0)
- TeachingDemos (2.9)
