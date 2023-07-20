library(ggraptR)
library(poolSeq)
library(patchwork)
library(dplyr)
library(ggridges)
library(tidyr)
library(viridis)
library(colorspace)
library(gganimate)
library(gifski)
library(av)
library(png)
library(magick)
library(scico)
library(plotrix)
library(ggdensity)
library(ggplot2) # needs to be version ≥ 2.1.0
library(scales)
library(devtools)
library(network)
library(sna)
library(GGally)
library(geomnet)
library(ggnetwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(vcfR)
library(ggplotify)
library(pheatmap)
library(ggbreak)
dev.off()

#######################################

data <- fread("alpha_VG.csv")
data$popalpha = 1000 * data$alpha
tmp <- data[1:4,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
data <- rbind(data,tmp)
rm(tmp)
# ggraptR(data)
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
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.25,0.5,0.75,1,2,5), 
                     limits = c(0,5),
                     labels = c(0,0.25,0.5,0.75,1,2,5)) +
  scale_y_continuous("P(detected)",
                     breaks = c(0.01,0.02,0.04,0.06,0.08)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",4),"red",rep("black",2))))

#####################################################
data <- fread("alpha_mult.csv")
data$popalpha = 1000 * data$alpha
tmp <- data[1:3,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
data <- rbind(data,tmp)
rm(tmp)
# ggraptR(data)
ggplot(data = data, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(mult)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(mult)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Starting freq multiplier",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.25,0.5,0.75,1,2,5), 
                     limits = c(0,5),
                     labels = c(0,0.25,0.5,0.75,1,2,5)) +
  scale_y_continuous("P(detected)",
                     breaks = c(0.01,0.02,0.04,0.06)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",4),"red",rep("black",2))))

############################################3

data <- fread("alpha_start.csv")
data$popalpha = 1000 * data$alpha
tmp <- data[1:5,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
data <- rbind(data,tmp)
rm(tmp)
# ggraptR(data)
ggplot(data = data, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(start)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(start)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Starting freq",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.25,0.5,0.75,1,2,5), 
                     limits = c(0,5),
                     labels = c(0,0.25,0.5,0.75,1,2,5)) +
  scale_y_continuous("P(detected)",
                     breaks = c(0.01,0.02,0.04,0.06,0.08,0.10)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",4),"red",rep("black",2))))

###############################################

data <- fread("alpha_VG.csv")
data$popalpha = 1000 * data$alpha
tmp <- data[1:4,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
data <- rbind(data,tmp)
rm(tmp)

neutdata <- data %>% filter(popalpha == 0, VG == 0.000101)
neutdata$reps <- 1
neutdata$Pdetectedatleastonce <- neutdata$Pdetected 
og <- neutdata
for(reps in 2:10){
  tmp <- og
  Pdetectedatleastonce <- 1-(1-tmp$Pdetected)^reps
  tmp$reps <- reps
  tmp$Pdetectedatleastonce <- Pdetectedatleastonce
  neutdata <- rbind(neutdata, tmp)
}

data <- data %>% filter(VG == 0.000101)
data$reps <- 1
data$Pdetectedatleastonce <- data$Pdetected 
og <- data
for(reps in 2:10){
  tmp <- og
  Pdetectedatleastonce <- 1-(1-tmp$Pdetected)^reps
  tmp$reps <- reps
  tmp$Pdetectedatleastonce <- Pdetectedatleastonce
  data <- rbind(data, tmp)
}

ggplot(data = data, aes(x = reps, y = Pdetectedatleastonce)) + 
  geom_line(aes(color = as.factor(popalpha)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(popalpha)),
             alpha=1,
             size = 2.5) +
  geom_line(data = neutdata, 
            aes(x = reps, y = Pdetectedatleastonce), 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Population Scaled\n Selection Coefficient",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Number of Replicates", 
                     breaks = 1:10, 
                     labels = 1:10) +
  scale_y_continuous("P(detected in at least one replicate)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
