########################################################
## Title: Indentify CpGs that are stably predictive of gestational age
## Description: 
## Data: Newborn cord blood EPIC DNAm data
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.01
## Date modified:2023.02.23

## Load DNAm data
# DNAm <- Data frame with DNAm data: samples (rows) by CpGs (columns)
cpg_names <- colnames(DNAm)

## Load results from stability selection and permutation
stabsel_result <- readRDS("stabsel_result.Rdata")
permut_result <- readRDS("permut_result.Rdata")

## Compute average number of selected variables (q_stab)
p_list <- list()
for ( i in 1 :length(result))
{
  res <- as.numeric(result[[i]])
  p_list[[i]] <- length(res)
}
p_out <- do.call(cbind,p_list)
q_stab <- mean(p_out)
q_stab # Average number of selected variables

## Compute the selection probabilities (select_prob)
n_select <- rep(0, length(cpg_names))
t_select <- 0*n_select

for ( i in 1:length(stabsel_result))
{
  tt <- names(stabsel_result[[i]])[-1]
  if(length(tt)==0){
    
  }else{
    t_select[which(cpg_names %in% names(stabsel_result[[i]])[-1])] <- 1
  n_select <- n_select+t_select
  t_select <- 0*n_select
  }
}

out <- n_select/length(stabsel_result)
df_out <- data.frame(cpg_names = cpg_names, pi_select = out)
select_prob <- df_out[order(df_out$pi_select, decreasing = TRUE),]

# Save results
save(select_prob, file = "select_prob.RData")

## Determination of selection probability threshold
p_var <- length(cpg_names) 

finding_thresh <- function(q, p, E_v)
{
  pi_thresh = ((q^2)/(p*E_v))+1/2
  return(min(pi_thresh,1))
}

# Create table of thresholds
thresh <- numeric(10)
for (i in 1:10){
  thresh[i] <- finding_thresh(q=q_stab, p=p_var, E_v=i)
}

n_cpgs <- numeric(10)
for (i in 1:10){
  n_cpgs[i] <- nrow(select_prob[which(select_prob$pi_select > thresh[i]),])
}

table <- data.frame(
  EV = c(1:10), 
  thresh = thresh,
  n_cpgs = n_cpgs)

# Threshold when allowing for a maximum of 2 false discoveries (E_v)
# Comment: The optimal threshold should be chosen based on the data, balancing E(V) and number of CpGs
thresh_EV2 <- finding_thresh(q=q_stab, p=p_var, E_v=2)

# Subset results above threshold
stable_set <- select_prob[which(select_prob$pi_select > thresh_EV2),]
# Save threshold and stable set
save(thresh_EV2, file="thresh_EV2.RData")
write.csv(stable_set, file="stable_cpgs_EV2.csv")
