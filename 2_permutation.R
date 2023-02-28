########################################################
## Title: Permutation of CpGs predictive of gestational age
## Description: 
## Data: Newborn cord blood EPIC DNAm data
## 
## Author: Kristine L Haftorn
## Date created: 2022.08.29
## Date modified: 2022.02.23 

## Load libraries
require(methods)
require(glmnet)

## Load and prepare data:
# DNAm <- Data frame with DNAm data: samples (rows) by CpGs (columns)
# GA <- Vector of GA matching the samples in DNAm

## Define x and y, and load LASSO model saved in 1_stability_selection (mod.cv)
x <- as.matrix(DNAm)
y <- as.numeric(GA)

load(file = "mod_cv.Rdata")

## Permutation
# Define subsample size
n_sample <- length(y)
n_stab <- floor(n_sample/2)
# Define number of times the subsampling and LASSO should be run
B <- 1000

# Define function for permutation
permut <- function(i){
  set.seed(i)
  
  # select subsample 
  ind_sub <- sample (1:n_sample, size = n_stab, replace = FALSE)
  x_sub <- x[ind_sub,]
  y_sub <- y[ind_sub]
  # rearrange GA values
  y_sub <- sample(y_sub) 
  # run LASSO
  res_sub <- glmnet(x = x_sub, y = y_sub, standardize = T, alpha = 1, lambda = mod.cv$lambda.1se)
  # Extract CpGs selected for prediction
  res <- coef(res_sub)
  res <- res[which(res[,1]!=0),]
  return(res)
}

# Run permutation
require(parallel)
result <- mclapply(1:B, FUN = permut, mc.cores = 10)

## Save results
saveRDS(result, file="permut_result.Rdata")
