## SCRIPTS

### `1_stability selection.R`

Stability selection procedure (lasso combined with subsampling).
Output: 
mod.cv.RData (lasso model)
stabsel_result.RData (stability selection result)

### `2_permutation.R`

Run stability selection agan with the same lasso model but with permutated data.
Output:
permut_result.RData (permutation result)

### `3_identify_stable_cpgs.R`

Compute selection probabilities and determine selection probability threshold.
Output:
select_prob.RData (selection probabilities for all CpGs included in the analysis)
thresh_EV2.RData (threshold when allowing for a maximum of 2 false discoveries)
stable_cpgs_EV2.csv (set of stably predictive CpGs as defined by the above threshold)

### `4_stabsel_figure1.R`

Create stability selection plot.
Output:
stabsel_figure1.jpeg (Stability selection plot)

### `5_clock_comparison_figure2.R`

Comparing results with three established GA clocks: Haftorn, Bohlin and Knight clocks.
Create plots showing the selection probabilities for the CpGs in each clock.
Output:
Figure2_comb.jpeg (Stability selection plots for the three clocks)

### `6_DNAm-GA_figure3.R`

Run GAM on gestational age and create GAM plot for each of the stably predictive CpGs.
Output: 
Figure3_DNAm-GA.jpeg (GAM plots for each of the stably predictive CpGs)
cpg_order_r2.RData (list of stably selected CpGs ordered by decreasing R2)

### `7_create_clocks_compare_figure4_5.R`

Create gestational age clocks from stably predictive CpGs and compare their performance to a standard lasso clock.
Output:
lasso_clock_mod_cv.RData (lasso clock model)
lasso_clock_result (lasso clock)
stable clock_1-15_cpg.RData (stable clocks)
Figure4_comb.jpeg (R2 and MAD plots comparing the stable clocks and the lasso clock)
Figure5_comb.jpeg (Plot showing the predictive performance of the "5 stable CpG clock" and the "15 stable CpG clock"

### `8_compare_gam_lasso.R`

Comparison of the performance in predicting gestational age using a GAM model versus a lasso model.
Output:
comparison_gam_lasso.jpeg (plot showing the predictive performance of the "15 stable CpG clock" and the standard lasso clock)

