########################################################
## Title: Stability selection plot
## Description: 
## Data: Newborn cord blood EPIC DNAm data
## 
## Author: Kristine L Haftorn
## Date created: 2022.09.27
## Date modified: 2023.02.23
##

## Load libraries
require(data.table)
require(ggplot2)
require(dplyr)

## Load data generated in 3_identify_stable_cpgs
load(file="thresh_EV2.RData") # Chosen selection probability threshold
load(file = "select_prob.RData") # CpGs with selection probabilities
cpgs <- select_prob
colnames(cpgs)[1] <- "Name"

# Load manifest file and merge with stabsel results
anno <- fread("anno_MethylationEPIC_v-1-0_B4.csv", data.table = F)
anno$Name <- as.character(anno$Name)
result <- subset(anno,anno$Name %in% cpgs$Name)
result <- merge(x = result, y = cpgs[,c("Name","pi_select")], by = "Name")
result <- anno
result$CHR <- as.numeric(result$CHR)

# Compute cumulative position of CpGs
result$MAPINFO_cum <- as.numeric(result$MAPINFO)/100000
max_table <- aggregate(result$MAPINFO_cum, by = list(result$CHR),FUN = "max")
for(i in 2:22){result$MAPINFO_cum[result$CHR==i] <- result$MAPINFO_cum[result$CHR==i] + sum(max_table$x[1:(i-1)])}

# Prepare x axis
xlabel_pos <- aggregate(result$MAPINFO_cum, by = list(result$CHR),"mean")
xlabel_pos$Group.1[xlabel_pos$Group.1==23] <- "X"
xlabel_pos$Group.1[xlabel_pos$Group.1==24] <- "Y"
new_x_labels <- as.character(1:22)
new_x_labels[seq.int(11, 22, 2)] <- ""

# Create variables for coloring
# Load clock CpGs
knight <- read.csv(file="Knight_clock.csv")
knight <- knight[2:149,1]
bohlin <- read.csv(file="Bohlin_clock.csv")
bohlin <- bohlin[2:97,2]
haftorn <- read.csv(file="EPIC_GA_clock.csv")
haftorn <- haftorn[2:177,3]
clock <- c(haftorn, bohlin, knight)
clock <- unique(clock)

result$clock <- result$Name %in% clock
result$col <- "grey"
result$col[result$CHR %in% c(2,4,6,8,10,12,14,16,18,20,22)] <- "darkgrey"
result$col[result$clock == TRUE] <- "clock"
result$col[result$pi_select > thresh_EV2] <- "sign"
table(result$col)

result$trans <- 1
result$trans[result$col %in% c("grey","darkgrey")] <- 0.2
table(result$trans)

result$shape <- "circle"
result$shape[result$clock == TRUE] <- "star"
result$shape <- as.factor(result$shape)
table(result$shape)

## Make plot
plot <- ggplot(result %>% arrange(trans), aes(x=MAPINFO_cum, y=pi_select, col=col)) +
  geom_point(aes(shape=shape)) +
  scale_colour_manual(values=c("grey"="grey", "darkgrey"="grey50", "sign"="#D55E00", "clock"="#0072B2")) + #"clock"="#0072B2"
  scale_shape_manual(values=c("circle"=16, "star"=8)) +
  geom_hline(yintercept=0.729, size=1.2, linetype="dashed") +
  geom_hline(yintercept=0.5, size=1) +
  scale_x_continuous(label=new_x_labels, breaks=xlabel_pos$x) +
  scale_y_continuous(limits=c(0,1) , expand = c(0,0.1)) +
  ggtitle("Stability selection of CpGs predictive of GA") +
  xlab("Chromosome") +
  ylab("Selection probability") +
  theme_classic(base_size=18) +
  theme(
    legend.position = "none",
    title = element_text(face="bold", size=16),
    axis.text=element_text(face="bold", colour = "black"),
    axis.line = element_line(size=1.2),
    axis.ticks= element_line(size=1.2)
  )
# Save plot
ggsave(plot, file="stabsel_figure1.jpeg", width=8, height=5, dpi=900)
