This folder contains the "stable CpG GA clocks" that were generated in the article. 
They can be applied to DNA methylation datasets from cord blood to predict gestational age in newborns.

Below is a description of how to apply the clocks.

## How to apply "stable CpG GA clocks" to a DNA methylation dataset
The clocks are essentially gam objects generated using the mod.gam function in the {mgcv} package.
Therefore, you need to install and load this package to be able to use the clocks.

```r
# Load packages
require(mgcv)

# Load data
load("DNAm_data.RData") # Your DNAm data, which should contain sample IDs (rows) and CpG IDs (columns)
load("stable_clock_15_cpg.RData") # We use the 15 stable CpG clock for this example

# Extract list of CpGs in the clock
cpgs <- attr(terms(mod.gam), "term.labels")

# Predict gestational age
pred <- predict(mod.gam, newdata = as.data.frame(DNAm_data[,cpgs]))

```

