########################################################
## Title: DNAm-GA relationships for stably predictive CpGs
## Description: Run robust regression on GA and create regression plot for each of the stably predictive CpGs
## Data: Newborn cord blood EPIC DNAm data (training set)
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.29
## Date modified: 2023.02.23
## 

## Load packages
require(robustbase)
require(latex2exp)
require(ggplot2)
require(cowplot)

## Load data
# DNAm <- Data frame with DNAm data: samples (rows) by CpGs (columns) - training data if used for creating a clock/prediction model 
# GA <- Vector of GA matching the samples in DNAm
stable_set <- (file="stable_cpgs_EV2.csv") # Generated in 3_identify_stable_cpgs.R

stable_DNAm <- DNAm[,stable_set$cpg_names]

## Run robust regression and create plot for each cpg

stable_DNAm <- as.data.frame(stable_DNAm)
dnam_list <- as.list(stable_DNAm)

cpg_list <- list()
for (i in 1:length(dnam_list)){
  cpg_list[[i]] <- data.frame(dnam = dnam_list[[i]],
                              GA = GA_tr)
}
names(cpg_list) <- names(dnam_list)

cpg_plots <- list()
r2_list   <- numeric(length(cpg_list))
for (i in 1:length(cpg_list)) {
  lmrob <- lmrob(cpg_list[[i]]$dnam~cpg_list[[i]]$GA)
  intercept <- summary(lmrob)$coefficients[1,1]
  slope <- summary(lmrob)$coefficients[2,1]
  pval <- signif(summary(lmrob)$coefficients[2,4], digits = 3)
  r2_list[i] <- round(summary(lmrob)$adj.r.squared, digits = 3)
  
  cpg_plots[[i]] <- ggplot(data=cpg_list[[i]], aes(x=GA, y=dnam)) +
    geom_point() +
    xlab("Gestational age (ultrasound)")+
    ylab("DNAm") +
    ylim(0, 1) +
    xlim(215,300) +
    geom_abline(intercept = intercept, slope = slope, size = 1.5, color = "#D55E00") +
    ggtitle(names(cpg_list)[i]) +
    annotate("text", x=270, y=1, label=TeX(paste("$R^2 = ",r2_list[i],", p < ",pval,"$", sep="")), size = 5, fontface="bold") +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      title = element_text(face="bold", size=16),
      axis.text=element_text(face="bold", colour = "black"),
      axis.line = element_line(size=1.2),
      axis.ticks= element_line(size=1.2)
    )
}

## Create combined plot
r2_order <- order(r2_list, decreasing = TRUE)
p <- cpg_plots[r2_order]
plot_comb <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],
                       p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]], 
                       labels = letters[1:15],
                       label_size = 20,
                       ncol = 4, nrow = 4) 

## Save plot and R2 order  
ggsave(plot_comb, file="Figure3_DNAm-GA.jpeg", width=16, height=16, dpi=300) 
cpg_order_r2 <- names(cpg_list[r2_order])
save(cpg_order_r2, file="cpg_order_r2.RData")


