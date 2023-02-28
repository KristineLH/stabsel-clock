########################################################
## Title: Stably predictive clocks
## Description: Create clocks from stably predictive CpGs and compare their performance to a standard LASSO clock
## Data: Newborn cord blood EPIC DNAm data (training and test set)
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.29
## Date modified: 2023.02.27

# Load packages
require(glmnet)
require(robustbase)
require(mgcv)
require(ggplot2)
require(cowplot)

## Load data
# DNAm_tr <- Data frame with DNAm data for the training set: samples (rows) by CpGs (columns)
# GA_tr <- Vector of GA matching the samples in DNAm_tr
# DNAm_test <- Data frame with DNAm data for the test set: samples (rows) by CpGs (columns)
# GA_test <- Vector of GA matching the samples in DNAm_test
stable_set <- (file="stable_cpgs_EV2.csv") # Generated in 3_identify_stable_cpgs.R
load(file="cpg_order_r2.RData") # Generated in 6_DNAm-GA_figure3.R

stable_tr <- DNAm_tr[,cpg_order_r2]
stable_test <- DNAm_test[,cpg_order_r2]

## Create standard LASSO clock
x <- as.matrix(DNAm_tr)
y <- as.numeric(GA_tr)
set.seed(1105)
mod.cv <- cv.glmnet(y = y,x = x, alpha = 1, nfold = 10)
result <- glmnet(x = x, y = y, standardize = T, alpha = 1, lambda = mod.cv$lambda.1se)
save(mod.cv, file="lasso_clock_mod.cv.Rdata")
save(result, file="lasso_clock_result.Rdata")

lasso_clock <- data.frame(as.matrix(coef(result)))
lasso_clock$cpgs <- rownames(lasso_clock)
lasso_clock <- subset(lasso_clock,lasso_clock$s0 != 0)
n_cpgs_lasso <- nrow(lasso_clock)-1

# Predict GA with standard LASSO clock
newx=as.matrix(DNAm_test)

pred <- predict(mod.cv, newx = newx, s = mod.cv$lambda.1se)
summary(pred)
test <- cbind(GA_test, pred)
test <- as.data.frame(test)
colnames(test)[2]<-"pred"

# Check predctive power
lmrob <- lmrob(test$GA_test~test$pred)

r2_list <- numeric(length(cpg_order_r2)+1)
mad_list <- numeric(length(cpg_order_r2)+1)

r2_list[16] <- round(summary(lmrob)$adj.r.squared, digits = 3)
mad_list[16] <- round(median(abs(test$GA_test-test$pred)), digits = 3)

## Create stable CpG clocks
tr <- cbind(stable_tr, GA_tr)
tr = as.data.frame(tr)

cpg.s <- paste("-s(", cpg_order_r2,")", sep="")
mod <- as.formula(paste("GA_tr", paste("s(", cpg_order_r2,")", sep="", collapse = "+"), sep="~"))

plots <- list()
for (i in length(cpg_order_r2):1){
  # Fit model to training set
  mod.gam <- gam(mod, data = tr)

  # Predict in test set
  pred <- predict(mod.gam, newdata = stable_test[,1:i,drop=FALSE])
  test <- cbind(GA_test, pred)
  test <- as.data.frame(test)
  
  # Test predictive power
  lmrob=lmrob(test$GA_test~test$pred, k.max = 500)
  intercept <- summary(lmrob)$coefficients[1,1]
  slope <- summary(lmrob)$coefficients[2,1]
  r2_list[i] <- round(summary(lmrob)$adj.r.squared, digits = 3)
  mad_list[i] <- format(round(median(abs(test$GA_test-test$pred)), digits = 2), nsmall=2)
  
  # Make plot
  plots[[i]] <- ggplot(data=test, aes(x=pred, y=GA_test)) +
    geom_point(shape=1, size=3, stroke=2) +
    xlab("Gestational age (DNAm)")+
    ylab("Gestational age (ultrasound)") +
    ylim(220, 301) +
    xlim(220,301) +
    geom_abline(intercept = intercept, slope = slope, size = 1.5, color = "#D55E00") +
    ggtitle(paste(i,"stable CpG GA clock (n = 429)", sep=" ")) +
    annotate("text", x = 240, y = 300, label=TeX(paste("$R^2 = ",r2_list[i],", MAD = ",mad_list[i],"$", sep="")), size = 10, fontface="bold") +
    theme_classic(base_size = 30) +
    theme(
      legend.position = "none",
      title = element_text(face="bold"),
      axis.text=element_text(face="bold", colour = "black"),
      axis.line = element_line(size=1.5),
      axis.ticks= element_line(size=1.5)
    )

  # Update formula
  mod <- update.formula(mod, paste("~ .",cpg.s[i], sep=" ", collapse = ""))    
}

## Combine results in data.frame
cpg_nr <- c(1:16)
df <- data.frame(cpg_nr = c(1:16), r2 = r2_list, mad = mad_list)
df$col <- "black"
df$col[16] <- "red"

# Make r2 and mad plots
r2 <- ggplot(data=df, aes(x=cpg_nr, y=r2, col=col)) +
  geom_point(size=3) +
  xlab("Number of CpGs in clock")+
  ylab(bquote('R'^2)) +
  ylim(0, 1) +
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "233")) +
  scale_colour_manual(values=c("black"="black", "red"="red")) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text=element_text(face="bold", colour = df$col),
    axis.line = element_line(size=1.5),
    axis.ticks= element_line(size=1.5)
  )
ggsave(r2, file="Figure4a_R2.jpeg", width=7, height=5, dpi=300)

mad <- ggplot(data=df, aes(x=cpg_nr, y=mad, col=col)) +
  geom_point(size=3) +
  xlab("Number of CpGs in clock")+
  ylab("MAD (days)") +
  ylim(3, 6) +
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "233")) +
  scale_colour_manual(values=c("black"="black", "red"="red")) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text=element_text(face="bold", colour = df$col),
    axis.line = element_line(size=1.5),
    axis.ticks= element_line(size=1.5)
  )
ggsave(mad, file="Figure4b_MAD.jpeg", width=7, height=5, dpi=300)

## Create and save combined plot of r2 and mad
plot_comb <- plot_grid(r2, mad,
                       labels = c("a", "b"),
                       label_size = 25,
                       ncol = 2, nrow = 1)
ggsave(plot_comb, file="Figure4_comb.jpeg", width=14, height=5, dpi=300)     

# Create and save combined plot of "4 stable CpG clock" and "15 stable CpG clock"
plot_clock_comb <- plot_grid(plots[[4]],plots[[15]], 
                       labels = c("a", "b"),
                       label_size = 40,
                       ncol = 2, nrow = 1) 
ggsave(plot_clock_comb, file="Figure5_comb.jpeg", width=20, height=10, dpi=300) 
