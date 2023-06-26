########################################################
## Description: GOmeth analysis of CpGs that are stably predictive of GA
## Data: Stably selected CpGs and all CpGs analysed
## 
## Author: Kristine L Haftorn
## Date created: 2023.06.22
## Date modified: 2023.06.26
## 

## Load packages
require(missMethyl)

## Load data:
res <- read.csv("MAIN_ANALYSES/DATA/Supp_Data_1_stabsel.csv")
all_cpgs <- res$CpG_ID
sign_cpgs <- res$CpG_ID[1:24]

# Run GOmeth
# GO
go_res <- gometh(sig.cpg = sign_cpgs, all.cpg = all_cpgs, collection = "GO", plot.bias = TRUE, prior.prob = TRUE)
# Total number of GO categories significant at 5\% FDR
table(go_res$FDR<0.1)
# Table of top GO results
topGO(kegg_res)

# KEGG
kegg_res <- gometh(sig.cpg = sign_cpgs, all.cpg = all_cpgs, collection = "KEGG", plot.bias = TRUE, prior.prob = TRUE)
# Total number of KEGG categories significant at 5\% FDR
table(kegg_res$FDR<0.1)
# Table of top KEGG results
topGO(kegg_res)
