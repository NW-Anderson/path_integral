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

#######################################
setwd("~/Documents/GitHub/path_integral/alpha_VG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
tmp <- master[1:4,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
master <- rbind(master,tmp)
rm(tmp)

master <- master %>% filter(!(popalpha == 10 & VG == 0.0101 & Pdetected <= 0.01))

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
  scale_color_viridis(discrete = T)

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

############################################

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

############################################

setwd("~/Documents/GitHub/path_integral/alpha_start")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
tmp <- master[1:4,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
master <- rbind(master,tmp)
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
  scale_color_viridis(discrete = T)

############################################

setwd("~/Documents/GitHub/path_integral/alpha_time")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "time", "start", "thresh", "totalP", "Pdetected")
master <- filter(master, time < 0.25,
                 time > 0.025)
master$popalpha = 1000 * master$alpha
tmp <- master[1:4,]
tmp$popalpha <- 0
tmp$Pdetected <- 0.01
master <- rbind(master,tmp)
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
                                   colour = c(rep("black",4),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T)

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

#################################################

#set parameters
dt = .001
n = 1000
t = cumsum(rep(dt,n))

#Brownian motion
B = cumsum(rnorm(n,0,sqrt(dt)))

plot(t,B,type="l")

# master <- data.frame("class","rep", "X", "t")

for(rep in 1:10){
  #neutral Wright-Fisher
  Edx = function(x,t){0}
  Vdx = function(x,t){x*(1-x)}
  X = numeric(n)
  X[1] = .3
  for (i in 2:n) {
    X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt + sqrt(Vdx(X[i-1],t[i-1]))*rnorm(1,0,sqrt(dt))
    X[i] = max(0,X[i])
    X[i] = min(1,X[i])
  }
  
  # plot(t,X,type="l",ylim = c(0,1))
  if(rep == 1){
    master <- data.frame(class = 0,rep = rep, X = X,t = t)
  }else{
    master <- dplyr::bind_rows(master,
                               data.frame(class = 0,
                                          rep = rep, 
                                          X = X,
                                          t = t))
  }
  for(Numinteractions in c(10,100,200,500)){
    interactions = sample(2:n,Numinteractions) %>% sort()
    Edx = function(x,t){5*x*(1-x)}
    Vdx = function(x,t){x*(1-x)}
    X = numeric(n)
    X[1] = .3
    for (i in 2:n) {
      if(i %in% interactions){
        X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt 
        X[i] = max(0,X[i])
        X[i] = min(1,X[i])
      } else{
        X[i] = X[i-1] + sqrt(Vdx(X[i-1],t[i-1]))*rnorm(1,0,sqrt(dt))
        X[i] = max(0,X[i])
        X[i] = min(1,X[i])
      }
    }
    master <- dplyr::bind_rows(master,
                               data.frame(class = Numinteractions,
                                          rep = rep, 
                                          X = X,
                                          t = t))
  }
  
}

Edx = function(x,t){5*x*(1-x)}
Vdx = function(x,t){x*(1-x)}
X = numeric(n)
X[1] = .3
for (i in 2:n) {
  X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt
  X[i] = max(0,X[i])
  X[i] = min(1,X[i])
}
master <- dplyr::bind_rows(master,
                           data.frame(class = 1000,
                                      rep = 1, 
                                      X = X,
                                      t = t))
# ggplot(master, aes(x = t, y = X)) + 
#   geom_line(aes(color = as.factor(class),
#                 alpha = as.factor(rep)))

tmp <- master %>% # filter(class %in% c(0,500,1000)) %>%
  filter(rep <= 3) %>%
  group_by(class,rep) %>%
  mutate(lineID = cur_group_id()) %>%
  ungroup() %>% 
  mutate(deterministic = (class == 1000))

detline <- tmp %>% filter(deterministic)
stochline <- tmp %>% filter(!deterministic)

# #genic selection Wright-Fisher
for(rep in 1:10){
  #neutral Wright-Fisher
  Edx = function(x,t){0}
  Vdx = function(x,t){x*(1-x)}
  X = numeric(n)
  X[1] = .3
  for (i in 2:n) {
    X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt + sqrt(Vdx(X[i-1],t[i-1]))*rnorm(1,0,sqrt(dt))
    X[i] = max(0,X[i])
    X[i] = min(1,X[i])
  }
  
  # plot(t,X,type="l",ylim = c(0,1))
  if(rep == 1){
    masterdiff <- data.frame(class = 0,rep = rep, X = X,t = t)
  }else{
    masterdiff <- dplyr::bind_rows(masterdiff,
                                   data.frame(class = 0,
                                              rep = rep, 
                                              X = X,
                                              t = t))
  }
}

diffline <- masterdiff %>% filter(rep <= 3)
# Edx = function(x,t){5*x*(1-x)}
# Vdx = function(x,t){x*(1-x)}
# X = numeric(n)
# X[1] = .3
# for (i in 2:n) {
#   X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt + sqrt(Vdx(X[i-1],t[i-1]))*rnorm(1,0,sqrt(dt))
#   X[i] = max(0,X[i])
#   X[i] = min(1,X[i])
# }
# 
# ggplot(stochline, aes(x = t, y = X)) +
#   geom_line(aes(color = class,
#                 group = as.factor(lineID)),
#                 alpha = 0.6,
#                 linewidth = 0.5) +
#   theme_bw() +
#   guides(color=guide_legend(title = "Number of\n interactions",
#                             override.aes = list(alpha=1))) +
#   # scale_x_break(c(1,2), scales = 0.9) 
#   scale_x_continuous("time (genomic units)") +
#   scale_y_continuous("Allele Frequency") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_color_viridis(end = 0.8) + 
#   geom_line(data = detline,
#             aes(x = t, y = X),
#             color = viridis(11)[11],
#             size = 1.5) +
#   geom_line(data = diffline, 
#             aes(x = t, y = X,
#                 group = rep), 
#             linetype="dashed", 
#             color = viridis(11)[11], size  = 0.75) 
# 
# firstplot <- ggplot(detline,aes(x = t, y = X)) + 
#   theme_bw() +
#   guides(color=element_blank()) +
#   # scale_x_break(c(1,2), scales = 0.9) 
#   scale_x_continuous("time (genomic units)") +
#   scale_y_continuous("Allele Frequency",
#                      limits = c(0,1)) +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   geom_line(color = viridis(11)[11],
#             size = 1.5) 

firstplot

secondplot <- firstplot +
  geom_line(data = diffline, 
            aes(x = t, y = X,
                group = rep), 
            linetype="dashed", 
            color = viridis(11)[11], size  = 0.75) 

secondplot

neutline <- stochline %>% filter(class == 0)
thirdplot <- secondplot + 
  geom_line(data = neutline,
            aes(group = as.factor(lineID)),
            alpha = 0.6,
            linewidth = 0.5,
            color = viridis(11)[1]) 

thirdplot

tmp <- master %>% # filter(class %in% c(0,500,1000)) %>%
  filter(rep %in% sample(1:10,3)) %>%
  group_by(class,rep) %>%
  mutate(lineID = cur_group_id()) %>%
  ungroup() %>% 
  mutate(deterministic = (class == 1000))
pertline <- tmp %>% filter(!deterministic,
                           class == 200) 


thirdplot + geom_line(data = pertline,
                      aes(group = as.factor(lineID)),
                      alpha = 0.6,
                      linewidth = 0.5,
                      color = viridis(11)[6]) 

##################################
# #genic selection Wright-Fisher
# Edx = function(x,t){5*x*(1-x)}
# Vdx = function(x,t){x*(1-x)}
# X = numeric(n)
# X[1] = .3
# for (i in 2:n) {
#   X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt + sqrt(Vdx(X[i-1],t[i-1]))*rnorm(1,0,sqrt(dt))
#   X[i] = max(0,X[i])
#   X[i] = min(1,X[i])
# }
# 
# lines(t,X,type="l", col = "red")

# genic selection
# Edx = function(x,t){5*x*(1-x)}
# Vdx = function(x,t){x*(1-x)}
# X = numeric(n)
# X[1] = .3
# for (i in 2:n) {
#   X[i] = X[i-1] + Edx(X[i-1],t[i-1])*dt
#   X[i] = max(0,X[i])
#   X[i] = min(1,X[i])
# }


# lines(t,X,type="l",col = "blue")


#########################################################

traj<-function(p,w11,w12,w22,tgen=200,plot.it=TRUE,add.it=FALSE,col="red"){
  
  #w11<-.1
  #w12<-1
  #w22<-.1
  p.array<-p
  for(i in 1:tgen){
    wbar<- w11*p^2 +w12*2*p*(1-p) + w22 * (1-p)^2
    margin<-(w11*p + w12*(1-p)) - (w12*p + w22*(1-p))
    d_p<-p*(1-p)*(margin) /(wbar)
    
    p<- p+d_p
    p.array<-c(p.array,p)
  }
  if(!add.it) plot(p.array,xlab="generations",ylab="Frequency of allele 1",type="l",lwd=3,col=col,ylim=c(0,1),cex.lab=1.5,cex.axis=1.5)
  if(add.it) lines(p.array,lwd=3,col=col)
  
  #if(plot.it==FALSE) return(p.array)
}

traj(p=0.03,w11=1, w12=.95, w22=.9,tgen=2000,col="red")

#########################################################

df <- data.frame(gen=c(),
                 geg=c(),
                 fig=c())
for(gen in c(50, 100,150)){
  for(geg in c(10,20,30)){
    df <- dplyr::bind_rows(df,
                           data.frame(gen,
                                      geg,
                                      path.expand(paste("2Na10_geg",geg,"_gen",gen,".png",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2Na10_geg...geg..._gen...gen....png...sep.......,
                    Generation = gen,
                    Gegenbauers = geg)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(Gegenbauers),
             cols=vars(Generation),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

#########################################################

df <- data.frame(gen=c(),
                 geg=c(),
                 fig=c())
for(gen in c(250, 500,750)){
  for(geg in c(10,20,30)){
    df <- dplyr::bind_rows(df,
                           data.frame(gen,
                                      geg,
                                      path.expand(paste("2na5_geg",geg,"_gen",gen,".png",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na5_geg...geg..._gen...gen....png...sep.......,
                    Generation = gen,
                    Gegenbauers = geg)

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(Gegenbauers),
             cols=vars(Generation),
             scales = "free",
             labeller = label_both) +
  theme_blank() +
  xlab("Ending Frequency") +
  ylab("Density")

#########################################################

df <- data.frame(sel=c(),
                 geg=c(),
                 fig=c())
for(sel in c(1, 5,10)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(sel,
                                      geg,
                                      path.expand(paste("2na",sel,"_geg",geg,"_gen20.png",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na...sel..._geg...geg..._gen20.png...sep.......,
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
