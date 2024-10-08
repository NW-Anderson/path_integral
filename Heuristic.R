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
library(ggimage)
library(rsvg)
library(ggrepel)
library(shadowtext)
library(Cairo)
library(extrafont)
dev.off()

####################################
####### pDetection Alpha VG #######
####################################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

pintDf <- master %>% filter(statDist < 0.05)

numDf <- master %>% filter(VG == 1e-4,
                           popalpha >= 15)
numDf$yval <- c(numDf[1,]$pintDetected,
                numDf[2,]$numDetected)

pintDf$VG <- 1000 * pintDf$VG
numDf$VG <- 1000 * numDf$VG

setwd("~/Documents/GitHub/path_integral/results/heuristic/Selection")
heurDf <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  heurDf <- bind_rows(heurDf, tmp)
}

heurDf <- heurDf %>% rename(s = V1,
                            pdet = V2,
                            pdetlog = V3) %>%
  mutate(s = 1000 * s) %>%
  reshape2::melt(., measure.vars = c("pdet",
                           "pdetlog"))

p1 <- ggplot(data = pintDf, aes(x = popalpha, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             size  = 0.75,
             alpha = 0.6) + 
  theme_bw() +
  scale_x_continuous(bquote(italic(2 * N[e] * alpha * Lambda * "/" *  V[S])), 
                     breaks = c(0, 1, 5, 10, 15, 20), 
                     labels = c(0, 1, 5, 10, 15, 20),
                     expand = c(0.05,0)) +
  scale_y_continuous("Probability detected\n(Q)",
                     # breaks = seq(0,0.3,by = 0.05),
                     expand = c(0,0), 
                     limits = c(0,max(heurDf$value)*1.05)) +
  geom_point(data = numDf, aes(x = popalpha, 
                               y = yval,
                               color = as.factor(VG)),
             alpha=c(0,1),
             size = 2.5,
             shape = 17,
             show.legend = F) +
  geom_line(data = numDf, aes(x = popalpha, 
                              y = yval,
                              color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5,
            linetype = "dashed")  + 
  geom_point(data = heurDf, aes(x = s,
                                y = value)) +
  geom_line(data = heurDf, aes(x = s,
                               y = value,
                               linetype = variable)) + 
  scale_linetype_discrete("Heuristic",
                           labels = c("Exp.",
                                      "Log.")) + 
  guides(color=guide_legend(title = bquote(italic(2 * N[e] * V[G] * "/" *  V[S])),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.position = c(.01,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 

p1

#########################
#### pDetection Time ####
#########################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedTimeVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

pintDf <- master %>% filter(statDist < 0.05)

numDf <- master %>% filter(VG == 1e-4,
                           time >= 0.15)
numDf$yval <- c(numDf[1,]$pintDetected,
                numDf[2,]$numDetected)

# pintDf <- pintDf %>% filter(VG == 0.1)

setwd("~/Documents/GitHub/path_integral/results/heuristic/time")
heurDf <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  heurDf <- bind_rows(heurDf, tmp)
}

heurDf <- heurDf %>% rename(t = V1,
                            pdet = V2,
                            pdetlog = V3) %>%
  mutate(t = t / 1000) %>%
  reshape2::melt(., measure.vars = c("pdet",
                                     "pdetlog"))

p2 <- ggplot(data = pintDf, aes(x = time, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             size  = 0.75,
             alpha = 0.6) + 
  theme_bw() +
  scale_x_continuous("Time\n(genomic units)",
                     expand = c(0.05,0)) +
  scale_y_continuous("P(detected)",
                     expand = c(0,0),
                     limits = c(0,max(heurDf$value) * 1.05)) +
  geom_point(data = numDf, aes(x = time,
                               y = yval,
                               color = as.factor(VG)),
             alpha=c(0,1),
             size = 2.5,
             shape = 17,
             show.legend = F) +
  geom_line(data = numDf, aes(x = time,
                              y = yval,
                              color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5,
            linetype = "dashed") +
  geom_point(data = heurDf, aes(x = t,
                                y = value)) +
  geom_line(data = heurDf, aes(x = t,
                               y = value,
                               linetype = variable)) + 
  scale_linetype_discrete("Heuristic",
                          labels = c("Exp.",
                                     "Log.")) + 
  guides(color = F,
         shape = F,
         linetype = F) +
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 

p2

##########################
#### pDetection Start ####
##########################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedStartVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

# pintDf <- master %>% filter(statDist < 0.05)
# 
# numDf <- master %>% filter(VG == 1e-4,
#                            time >= 0.15)
# numDf$yval <- c(numDf[1,]$pintDetected,
#                 numDf[2,]$numDetected)

setwd("~/Documents/GitHub/path_integral/results/heuristic/start")
heurDf <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  heurDf <- bind_rows(heurDf, tmp)
}

heurDf <- heurDf %>% rename(start = V1,
                            pdet = V2,
                            pdetlog = V3) %>%
  reshape2::melt(., measure.vars = c("pdet",
                                     "pdetlog"))

p3 <- ggplot(data = master, aes(x = start, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed",
             size  = 0.75) + 
  theme_bw() +
  scale_x_continuous("Starting frequency", 
                     breaks = c(0.025,0.05,0.1,0.15,0.20), 
                     labels = c(0.025,0.05,0.1,0.15,0.20),
                     expand = c(0.05,0)) +
  scale_y_continuous("P(detected)",
                     expand = c(0,0),
                     limits = c(0,max(heurDf$value) * 1.05)) +
  geom_point(data = heurDf, aes(x = start,
                                y = value)) +
  geom_line(data = heurDf, aes(x = start,
                               y = value,
                               linetype = variable)) + 
  scale_linetype_discrete("Heuristic",
                          labels = c("Exp.",
                                     "Log.")) + 
  guides(color = F,
         linetype = F) + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) +
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 


p3

########################
#### pDetected Main ####
########################

p1 + p2 + p3 + 
  plot_annotation(tag_levels = 'A')  & 
  theme(plot.tag = element_text(size = 12))

#########################