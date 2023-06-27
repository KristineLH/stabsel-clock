########################################################
## Title: Comparison between GAM and lasso
## Description: Comparison of the performance in predicting gestational age using a GAM model versus a lasso model
## Data: Newborn cord blood EPIC DNAm data (test set)
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.29
## Date modified: 2023.05.09
##

# Load packages
require(glmnet)
require(robustbase)
require(mgcv)
require(ggplot2)

## Load data
# DNAm_test <- Data frame with DNAm data for the test set: samples (rows) by CpGs (columns)
# GA_test <- Vector of GA matching the samples in DNAm_test
load(file="cpg_order_r2.RData") # Generated in 6_DNAm-GA_figure3.R

## Load clocks
# Lasso clock
load(file="lasso_clock_mod.cv.RData") # Generated in 7_create_clocks_compare_figure4_5.R
# GAM clock (15 stable CpG GA clock)
load(file="stable_clock_15_cpg.RData") # Generated in 7_create_clocks_compare_figure4_5.R

## Prediction with lasso clock
newx=as.matrix(DNAm_test)
pred_lasso=predict(mod.cv, newx=newx, s=mod.cv$lambda.1se)
summary(pred_lasso)
lasso <- cbind(GA_test, pred_lasso)
lasso <- as.data.frame(lasso)

lasso_mod=lmrob(lasso$GA_test~lasso$pred_lasso)
summary(lasso_mod)

## Prediction with GAM clock
stable_test <- as.data.frame(DNAm_test[,cpg_order_r2])
pred_gam <- predict(mod.gam, newdata=stable_test)
summary(pred_gam)
test <- cbind(lasso, pred_gam)
test <- as.data.frame(test)

gam_mod=lmrob(test$GA_test~test$pred_gam)
summary(gam_mod)

## Combine data (here you insert the number of samples in your test set instead of 429)
df <- data.frame(y = rep(test$GA_test,2),
                 x = c(test$pred_gam, test$pred_lasso),
                 method = as.factor(c(rep("GAM model (15 CpGs)",429), rep("Lasso model (233 CpGs)",429))))

### Create plot
mycolors <- c("GAM model (15 CpGs)"="#D55E00", "Lasso model (233 CpGs)"="#0072B2", "Ideal fit"="black")

comp <- ggplot(data=df, aes(x=x, y=y, col = method)) +
  geom_point(alpha = 0.4) +
  xlab("Gestational age (DNAm)")+
  ylab("Gestational age (ultrasound)") +
  ylim(220, 301) +
  xlim(220,301) +
  geom_abline(aes(intercept=0, slope=1, colour="Ideal fit"), size=2, show.legend=TRUE) +
  geom_abline(aes(intercept= 3.45777, slope=0.98691, colour="GAM model (15 CpGs)"), size=2, show.legend=TRUE) + 
  geom_abline(aes(intercept=-46.05759, slope=1.16295, colour="Lasso model (233 CpGs)"), size=2, show.legend=TRUE) +
  scale_colour_manual(values = mycolors) +
  theme_classic(base_size = 30) +
  theme(
    legend.position = c(0.75,0.15),
    legend.title = element_blank(),
    axis.text=element_text(face="bold", colour = "black"),
    axis.line = element_line(size=1.5),
    axis.ticks= element_line(size=1.5)
  )
ggsave(comp, file="comparison_gam_lasso.jpeg", width=10, height=10, dpi=300)
