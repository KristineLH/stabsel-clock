########################################################
## Title: Stability selection of CpGs predictive of gestational age
## Description: 
## Data: Newborn cord blood EPIC DNAm data
## 
## Author: Kristine L Haftorn
## Date created: 2022.08.29
## Date modified: 2023.02.23

## Load packages
require(methods)
require(glmnet)

## Load and prepare data:
# DNAm <- Data frame with DNAm data: samples (rows) by CpGs (columns)
# GA <- Vector of GA matching the samples in DNAm

## Fit LASSO model
x <- as.matrix(DNAm)
y <- as.numeric(GA)
set.seed(1105)
mod.cv <- cv.glmnet(x = x, y = y, alpha = 1, nfold = 10)

save(mod.cv, file = "mod_cv.RData")

## Stability selection
# Define subsample size
n_sample <- length(y)
n_stab <- floor(n_sample/2)
# Define number of times the subsampling and LASSO should be run
B <- 1000

# Define function for stability selection procedure
stabsel <- function(i){
  set.seed(i)
  
  # select subsample 
  ind_sub <- sample (1:n_sample, size = n_stab, replace = FALSE)
  x_sub <- x[ind_sub,]
  y_sub <- y[ind_sub]
  # run LASSO
  res_sub <- glmnet(x = x_sub, y = y_sub, standardize = T, alpha = 1, lambda = mod.cv$lambda.1se)
  # extract CpGs selected for prediction
  res <- coef(res_sub)
  res <- res[which(res[,1]!=0),]
  return(res)
}

# Run stability selection
require(parallel)
result <- mclapply(1:B, FUN = stabsel, mc.cores = 10)

## Save results
saveRDS(result, file = "stabsel_result.RData")


