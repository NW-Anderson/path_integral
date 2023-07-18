library(data.table)
library(ggraptR)
library(gridExtra)
library(viridis)
library(matrixStats)
library(tidyr)
library(patchwork)
library(ggbreak)
dev.off()
data <- fread("alpha_VG.csv")
data$popalpha = 1000 * data$alpha
ggraptR(data)
ggplot(data = data, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_break(c(2,3), scales = 1.5) 
