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
dev.off()

##################
#### Error VG ####
##################

setwd("~/Documents/GitHub/path_integral/results/numErrorAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "pintAUC", "numAUC")
master$popalpha = 1000 * master$alpha
master$error = abs(master$numAUC - master$pintAUC)
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = popalpha, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Population Scaled Selection Coefficient") +
  scale_y_continuous("Absolute Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

# mster <- master %>% filter(VG<=0.001)
# ggplot(data = mster, aes(x = popalpha, y = error)) +
#   geom_line(aes(color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5) +
#   geom_point(aes(color = as.factor(VG)),
#              alpha=1,
#              size = 2.5) +
#   theme_bw() +
#   guides(color=guide_legend(title = "Genetic Variance",
#                             override.aes = list(alpha=1))) +
#   scale_x_continuous("Population Scaled Selection Coefficient") +
#   scale_y_continuous("Absolute Error") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_color_viridis(discrete = T)

####################
#### Error Time ####
####################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorTimeVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("time", "VG", "pintAUC", "numAUC")
master$error = abs(master$numAUC - master$pintAUC)
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = time, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Time (Genomic Units)") +
  scale_y_continuous("Absolute Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

# mster <- master %>% filter(VG<=0.001)
# ggplot(data = mster, aes(x = time, y = error)) +
#   geom_line(aes(color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5) +
#   geom_point(aes(color = as.factor(VG)),
#              alpha=1,
#              size = 2.5) +
#   theme_bw() +
#   guides(color=guide_legend(title = "Genetic Variance",
#                             override.aes = list(alpha=1))) +
#   scale_x_continuous("Population Scaled Selection Coefficient") +
#   scale_y_continuous("Absolute Error") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_color_viridis(discrete = T)

#################
#### Error k ####
#################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorKVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("k", "VG", "pintAUC", "numAUC")
master$error = abs(master$numAUC - master$pintAUC)
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = k, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("kmax") +
  scale_y_continuous("Absolute Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

#####################
#### Error Start ####
#####################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorStartVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("start", "VG", "pintAUC", "numAUC")
master$error = abs(master$numAUC - master$pintAUC)
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = start, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Starting Frequency") +
  scale_y_continuous("Absolute Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

############################
#### Convergence 20 Gen ####
############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergence20gen")
df <- data.frame(sel=c(),
                 geg=c(),
                 fig=c())
for(sel in c(1, 5,10)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(sel,
                                      geg,
                                      path.expand(paste("2na",sel,"_", geg,"geg","_0.02time.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na...sel..._...geg...geg...._0.02time.svg...,
                    Gegenbauers = geg,
                    PopScaledSelection = sel)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(Gegenbauers),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

##############################
#### Convergence Alpha VG ####
##############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergenceAlphaVG")
df <- data.frame(sel=c(),
                 genVar=c(),
                 fig=c())
for(sel in c(1, 5,10)){
  for(genVar in c("0.0001","0.001","0.01")){
    df <- dplyr::bind_rows(df,
                           data.frame(sel,
                                      genVar,
                                      path.expand(paste("2na",sel,"._50geg","_",genVar,"VG.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na...sel...._50geg...._...genVar...VG.svg...,
                    VG = genVar,
                    PopScaledSelection = sel)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(VG),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

############################
#### Convergence 2Na 5 ####
############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergence2Na5")
list.files()
df <- data.frame(tme=c(),
                 geg=c(),
                 fig=c())
for(tme in c(0.5, 0.25, 0.75)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(tme,
                                      geg,
                                      path.expand(paste("2na5_", geg,"geg","_",tme,"time.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na5_...geg...geg...._...tme...time.svg...,
                    Gegenbauers = geg,
                    Time = tme)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(Gegenbauers),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

############################
#### Convergence 2Na 10 ####
############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergence2Na10")
list.files()
df <- data.frame(tme=c(),
                 geg=c(),
                 fig=c())
for(tme in c(0.05, 0.1, 0.15)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(tme,
                                      geg,
                                      path.expand(paste("2na10_", geg,"geg","_",tme,"time.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na10_...geg...geg...._...tme...time.svg...,
                    Gegenbauers = geg,
                    Time = tme)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(Gegenbauers),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

#######################
#### pDetection VG ####
#######################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
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
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)",
                     breaks = c(0.01,0.02,0.04,0.06,0.08)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

#########################
#### pDetection Time ####
#########################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaTime")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "time", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(time)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(time)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Time (Genomic Units)",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T)

##########################
#### pDetection Start ####
##########################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaStart")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
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
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))
