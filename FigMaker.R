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
master$error = (master$pintAUC - master$numAUC) / master$numAUC
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
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

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
master$error = (master$pintAUC - master$numAUC) / master$numAUC
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
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

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
master$error = (master$pintAUC - master$numAUC) / master$numAUC
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
  scale_x_continuous("k_max") +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

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
master$error = (master$pintAUC - master$numAUC) / master$numAUC
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
  scale_x_continuous("Starting Frequency",
                     breaks = (1:9 / 10)) +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

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
                    m_max = geg,
                    PopScaledSelection = sel)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 
ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") + 
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))



# ggplot(pointdf, aes(x = x, y = y, color = cols)) +
#   geom_line(aes(linetype = cols),
#              alpha = 0) + 
#   guides(color = guide_legend(override.aes = list(alpha = 1),
#                               nrow = 1)) +
#   theme_bw() + 
#   labs(color = "k_max",
#        linetype = "k_max") + 
#   scale_color_manual(values = the_colors, name = "k_max") +
#   scale_linetype_manual(values=c(rep("solid",6), "dashed"), name="k_max") + 
#   theme(legend.position="bottom") 
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

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 
ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(VG),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))

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
                    m_max = geg,
                    Time = tme)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))

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
                    m_max = geg,
                    Time = tme)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))

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

#######################
#### pDetection Ne ####
#######################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaNe")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

ggplot(data = master, aes(x = selCoef, y = Pdetected)) + 
  geom_line(aes(color = as.factor(Ne)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(Ne)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Ne",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Selection Coefficient", 
                     breaks = c(0, 5e-04, 1e-03, 5e-03, 1e-02), 
                     limits = c(0, 1e-02),
                     labels = c(0, 5e-04, 1e-03, 5e-03, 1e-02)) +
  scale_y_continuous("P(detected)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

#############
#### RFS ####
#############

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaNe")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

master <- master %>% filter(selCoef == 0.01)
master$reps = 1000 / master$Ne

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(j in 0:tmp$reps){
    foo <- data.table(Ne = tmp$Ne,
                      reps = paste(tmp$reps, " Replicates"),
                      reps_number = tmp$reps,
                      bin = j,
                      prob = choose(tmp$reps, j) * 
                        tmp$Pdetected^j * 
                        (1 - tmp$Pdetected)^(tmp$reps - j))
    df <- dplyr::bind_rows(df, foo)
  }
}
rm(foo, master, tmp, file, i, j)

dfdetected <- df %>% filter(bin > 0) %>% 
  group_by(reps) %>% 
  mutate(Pdetected = sum(prob)) %>% 
  mutate(prob = prob / Pdetected) %>% ungroup() 
dfnotdetected <- df %>% filter(bin == 0)
dfnotdetected$reps_number <- c(1,3,2.5,1.5)

ggplot(dfdetected, aes(y=prob, x=as.factor(bin))) + 
  geom_bar(stat="identity",
           fill = viridis(4)[2],
           alpha = 0.5) + 
  facet_wrap(~ reps, scales="free_x") + 
  theme_bw() + 
  xlab("Number of Replicates Detected") + 
  ylab("Probability Given Detected at least Once") + 
  geom_text(data = dfnotdetected, aes(x = reps_number, 
                                      y = 1.1, 
                                      label = paste("P(detected): ",round(1- prob,3))),
            hjust = 0, vjust = 1) + 
  theme(panel.grid.minor.x = element_blank(),
        axis.title = element_text(size=18),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))


######################
#### model figure ####
######################

rm(list=ls())

x <- seq(from = -2, to = 3, length.out = 500)

phenodist <- dnorm(x, sd = sqrt(10^-2))
phenodist <- phenodist / max(phenodist)

olddist <- dnorm(x, sd = 1)
olddist <- olddist / max(olddist)

newdist <- dnorm(x, mean = 1, sd = 1)
newdist <- newdist / max(newdist)

df <- dplyr::bind_rows(data.table(x = x,
                                  y = phenodist,
                                  class = "pheno",
                                  cols = 1),
                       data.table(x = x,
                                  y = olddist,
                                  class = "old",
                                  cols = 2),
                       data.table(x = x,
                                  y = newdist,
                                  class = "new",
                                  cols = 2))

areadf <- data.table(x = x,
                     y = phenodist,
                     class = "pheno",
                     cols = 1)

the_colors <- c(rgb(0.831964, 0.810543, 0.372854),
                rgb(0.35082, 0.595178, 0.853742))

ggplot(df, aes(x = x, y = y)) + 
  geom_line(aes(color = class,
                linetype = class,
                group = class),
            linewidth = 2,
            alpha = 0.75) +
  scale_color_manual(values = c(rep(the_colors[1],2), the_colors[2]),
                     labels=c("Shifted Fitness Function",
                              "Initial Fitness Function",
                              "Trait Distribution"),
                     name = "") +
  scale_linetype_manual(values = c("dotted", rep("solid",2)),
                        labels=c("Shifted Fitness Function",
                                 "Initial Fitness Function",
                                 "Trait Distribution"),
                        name = "") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size=18),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12)) + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("0" = "Old\nOptimum", 
                                "1" = "New\nOptimum"),
                     name = "Trait Value") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 2),
                              nrow = 1)) +
  geom_area(data = areadf, 
            fill = the_colors[2],
            alpha = 0.5)
