This folder contains the "stable CpG GA clocks" that were generated in the article. 
They can be applied to DNA methylation datasets from cord blood to predict gestational age in newborns.

Below is a description of how to apply the clocks.

## How to apply "stable CpG GA clocks" to a DNA methylation dataset
The clocks are essentially gam objects generated using the {mgcv} package.
Therefore, you need to install and load this package to be able to use the clocks. 

```r
# Load packages
require(mgcv)

# Load data
DNAm_data <- load("your_DNAm_data.RData")

```

