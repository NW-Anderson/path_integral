library(dplyr)
library(ggplot2)
###############################################
#set parameters
dt = .0001
n = 1000
t = cumsum(rep(dt,n))

#Brownian motion
master <- data.frame()
for(rep in 1:20){
  print(rep)
  B = cumsum(rnorm(n,0,sqrt(dt)))

  C = (1-t/last(t))*0.3 + t*0.7 / last(t) + (B - t*last(B) / last(t))
  
  
  if(min(C) > 0 & max(C) < 1){
    p <- sum(log(dnorm(C[-1],
                       mean = C[-length(C)], 
                       sd = sqrt(C[-length(C)]*(1-C[-length(C)]) * dt))*dt))
    master <- dplyr::bind_rows(master,
                               data.frame(t = t,
                                          C = C,
                                          log_amp = p,
                                          rep = rep))
  }
}


master %>% filter(rep %in% sample(unique(rep),5)) %>%
  ggplot(., aes(x = t, y = C, 
                group = rep,
                color = log_amp #,
                # alpha = log_amp
                )) + 
  geom_point(size = 0.4) +
  geom_line(linewidth = 1) +  
  # geom_hline(yintercept = 1) + 
  guides(color = F,
         alpha = F) +
  scale_x_continuous(name = "Time",
                     labels = c("0","t"),
                     breaks = c(0,last(t)),
                     expand = c(0,0),
                     limits = c(0,last(t))) + 
  scale_y_continuous(name = "Frequency",
                     breaks = c(0,1),
                     expand = c(0,0),
                     limits = c(0,1)) +
  scale_colour_gradientn(colours = c("black", "grey80"),
                         trans = "reverse") +
  # scale_color_continuous(trans = "reverse")
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm")) + 
  theme_bw() + 
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))  


################################################################


#set parameters
dt = .0001
n = 1000
t = cumsum(rep(dt,n))

#Brownian motion
master <- data.frame()
for(rep in 1:3){
  # print(rep)
  
  B1 <- cumsum(rnorm(n/2,0,sqrt(dt)))
  t1 <- t[1:500]
  
  C1 = (1-t1/last(t1))*0.3 + t1*0.2 / last(t1) + (B1 - t1*last(B1) / last(t1))
  
  B2 <- cumsum(rnorm(n/2,0,sqrt(dt)))
  t2 <- t1
  
  C2 = (1-t2/last(t2))*0.2 + t2*0.7 / last(t2) + (B2 - t2*last(B2) / last(t2))
  
  C <- c(C1,C2)
  
  plot(t,C,type = "l")
  
  if(min(C) > 0 & max(C) < 1){
    p1 <- sum(log(dnorm(C1[-1],
                       mean = C1[-length(C1)], 
                       sd = sqrt(C1[-length(C1)]*(1-C1[-length(C1)])* dt))*dt))
    
    p2 <- sum(log(dnorm(C2[-1],
                       mean = C2[-length(C2)], 
                       sd = sqrt(C2[-length(C2)]*(1-C2[-length(C2)])* dt))*dt))
    
    t2 <- t[501:1000]
    # tmp <- dplyr::bind_rows(data.frame(t = t1,
    #                                    C = C1,
    #                                    log_amp = p1,
    #                                    rep = rep,
    #                                    scats = 1),
    #                         data.frame(t = t2,
    #                                    C = C2,
    #                                    log_amp = p2,
    #                                    rep = rep,
    #                                    scats = 1))
    # master <- dplyr::bind_rows(master,
    #                            tmp)
    master <- dplyr::bind_rows(master, data.frame(t = t,
                                                  C = C,
                                                  log_amp = c(rep(p1,500),
                                                              rep(p2,500)),
                                                  rep = rep))
  }
}

print(unique(master$log_amp))
linedf <- data.frame(x = rep(t[500],2),
                     y = c(1,0.2))


master %>% # filter(rep %in% sample(unique(rep),1)) %>%
  ggplot(., aes(x = t, y = C, 
                group = rep,
                color = log_amp #,
                # alpha = log_amp
                )) + 
  # geom_point(size = 0.4) +
  geom_line(linewidth = 0.75) +
  guides(color = F,
         alpha = F) +
  scale_x_continuous(name = "Time",
                     labels = c("0","t"),
                     breaks = c(0,last(t)),
                     expand = c(0,0),
                     limits = c(0,last(t))) + 
  scale_y_continuous(name = "Frequency",
                     breaks = c(0,1),
                     expand = c(0,0),
                     limits = c(0,1)) + 
  scale_colour_gradientn(colours = c("black", "grey80"),
                         trans = "reverse") +
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm")) + 
  theme_bw() + 
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))  +
  geom_line(data = linedf, aes(x=x,y=y), 
            inherit.aes = F,
            linetype = "dashed",
            linewidth = 0.75,
            alpha = 0.5)
















