########################################################
## Title: Comparison between stably predictive CpGs and GA clocks
## Description: Comparing our results with three established GA clocks: Haftorn, Bohlin and Knight clocks
## Data: Newborn cord blood EPIC DNAm data
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.29
## Date modified: 2023.02.23

## Load libraries
require(data.table)
require(ggplot2)
require(dplyr)
require(robustbase)
require(cowplot)

## Load data generated in 3_identify_stable_cpgs
load(file="thresh_EV2.RData") # Chosen selection probability threshold
load(file = "select_prob.RData") # CpGs with selection probabilities
cpgs <- select_prob
colnames(cpgs)[1] <- "cpgs"

# Load clock data
knight <- read.csv(file="Knight_clock.csv")
knight <- knight[2:149,]
bohlin <- read.csv(file="Bohlin_clock.csv")
bohlin <- bohlin[2:97,2:3]
haftorn <- read.csv(file="EPIC_GA_clock.csv")
haftorn <- haftorn[2:177,2:3]

# Load DNAm data
# DNAm <- Data frame with DNAm data: samples (rows) by CpGs (columns)

# Compute DNAm variance for each clock
knight_DNAm <- DNAm[,intersect(colnames(DNAm), knight$cpgs)]
knight_var <- data.frame(cpgs = colnames(knight_DNAm), var = sapply(knight_DNAm, var))
rownames(knight) <- knight$cpgs
knight <- knight[knight_var$cpgs,]
knight <- merge(knight, knight_var, by = "cpgs")
knight$comb <- knight$s0*knight$var

bohlin_DNAm <- DNAm[,intersect(colnames(DNAm), bohlin$cpgs)]
bohlin_var <- data.frame(cpgs = colnames(bohlin_DNAm), var = sapply(bohlin_DNAm, var))
rownames(bohlin) <- bohlin$cpgs
bohlin <- bohlin[bohlin_var$cpgs,]
bohlin <- merge(bohlin, bohlin_var, by = "cpgs")
bohlin$comb <- bohlin$s0*bohlin$var

haftorn_DNAm <- DNAm[,intersect(colnames(DNAm), haftorn$cpgs)]
haftorn_var <- data.frame(cpgs = colnames(haftorn_DNAm), var = sapply(haftorn_DNAm, var))
rownames(haftorn) <- haftorn$cpgs
haftorn <- haftorn[haftorn_var$cpgs,]
haftorn <- merge(haftorn, haftorn_var, by = "cpgs")
haftorn$comb <- haftorn$s0*haftorn$var

# Merge stabsel and clock data
knight <- merge(knight, cpgs, by = "cpgs")
bohlin <- merge(bohlin, cpgs, by = "cpgs")
haftorn <- merge(haftorn, cpgs, by = "cpgs")

## Create plots

# Haftorn clock
summary(haftorn$comb)
haftorn$stable <- "no"
haftorn$stable[haftorn$pi_select > thresh_EV2]<- "yes"
Haftorn_p <- ggplot(data = haftorn, aes(x = comb, y = pi_select, col = stable)) +
  geom_point(shape = 8) +
  scale_colour_manual(values=c("no" = "#0072B2", "yes" = "#D55E00")) +
  xlab("Beta coefficient x DNAm variance")+
  ylab("Selection probability") +
  ylim(0, 1) +
  xlim(-0.15,0.15) +
  geom_hline(yintercept = thresh_EV2, size=1.2, linetype="dashed") +
  geom_hline(yintercept = 0.5, size=1) +
  ggtitle("The Haftorn clock (176 CpGs)") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    title = element_text(face = "bold", size = 14),
    axis.text=element_text(face="bold", colour = "black"),
    axis.line = element_line(size = 1.2),
    axis.ticks= element_line(size = 1.2)
  )
ggsave(Haftorn_p, file="Figure2a_Haftorn.jpeg", width=4, height=4, dpi=300)

# Bohlin clock
bohlin$stable <- "no"
bohlin$stable[bohlin$pi_select > thresh_EV2]<- "yes"
Bohlin_p <- ggplot(data = bohlin, aes(x = comb, y = pi_select, col = stable)) +
  geom_point(shape = 8) +
  scale_colour_manual(values=c("no" = "black", "yes" = "#D55E00")) +
  xlab("Beta coefficient x DNAm variance")+
  ylab("Selection probability") +
  ylim(0, 1) +
  xlim(-0.2,0.2) +
  geom_hline(yintercept = thresh_EV2, size=1.2, linetype="dashed") +
  geom_hline(yintercept = 0.5, size=1) +
  ggtitle("The Bohlin clock (86 CpGs)") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    title = element_text(face = "bold", size = 14),
    axis.text=element_text(face="bold", colour = "black"),
    axis.line = element_line(size = 1.2),
    axis.ticks= element_line(size = 1.2)
  )
ggsave(Bohlin_p, file="Figure2b_Bohlin.jpeg", width = 4, height = 4, dpi = 300)

# Knight clock
knight$stable <- "no"
knight$stable[knight$pi_select > thresh_EV2]<- "yes"
Knight_p <- ggplot(data = knight, aes(x = comb, y = pi_select, col = stable)) +
  geom_point(shape = 8) +
  scale_colour_manual(values=c("no"="black", "yes" = "#D55E00")) +
  xlab("Beta coefficient x DNAm variance")+
  ylab("Selection probability") +
  ylim(0, 1) +
  xlim(-0.05,0.05) +
  geom_hline(yintercept = thresh_EV2, size=1.2, linetype="dashed") +
  geom_hline(yintercept = 0.5, size=1) +
  ggtitle("The Knight clock (140 CpGs)") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    title = element_text(face = "bold", size = 14),
    axis.text=element_text(face="bold", colour = "black"),
    axis.line = element_line(size = 1.2),
    axis.ticks= element_line(size = 1.2)
  )
ggsave(Knight_p, file="Figure2c_Knight.jpeg", width = 4, height = 4, dpi = 300)

# Create and save combined plot
plot_comb <- plot_grid(Haftorn_p, Bohlin_p, Knight_p,
                       labels = c("a", "b", "c"),
                       label_size = 20,
                       ncol = 3, nrow = 1)
ggsave(plot_comb, file="Figure2_comb.jpeg", width=14, height=4, dpi=300) 
