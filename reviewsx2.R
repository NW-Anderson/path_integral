library(dplyr)
library(ggplot2)

start <- 0.01
s <- 0.01
N <- 10^3

# 2Ns = 20

getTraj <- function(start,s,N){
  x <- start
  while(x != 1){
    traj <- c()
    x <- start
    while(!(x %in% c(0,1))){
      xPrime <- x + s * x* (1 - x)
      x <- rbinom(1, 2*N, xPrime) / (2 * N)
      traj <- c(traj, x)
    }
    
  }
  return(traj)
}

master <- data.frame()
for(i in 1:5){
  traj <- getTraj(start, s, N)
  master <- dplyr::bind_rows(master,
                             data.frame(rep = i,
                                        freq = traj,
                                        gen = 1:length(traj)))
}

ggplot(master) + 
  geom_line(aes(x = gen,
                y = freq,
                color = as.factor(rep))) + 
  theme_bw()

P=(1-exp(-4 * N * s * start))/(1-exp(-4* N * s))

df <- data.frame(x = seq(0,1, by = 0.001)) %>%
  mutate(dens = dexp(x,  P / start))

ggplot(df) + 
  geom_line(aes(x = x,
                y = dens)) + 
  theme_bw()

######################################################

start = 0.1
s = 0.01
ne = 500

pfix = (1 - exp(-4 * ne * s * start)) / (1 - exp(-4 * ne * s))

stochStart = rexp(10^2, pfix / start)

t <- 0:100

master <- data.frame()
count <- 0
for(x in stochStart){
  f <- x * exp(s * t)
  count <- count + 1
  master <- bind_rows(master, 
                      data.frame(traj = f,
                                 rep = count,
                                 gen = 0:100))
}

ggplot(master) + 
  geom_line(aes(x = gen,
                y = traj,
                color = as.factor(rep),
                group = as.factor(rep))) +
  theme_bw()


data.frame(x = seq(0,2,by = 0.01),
           y = dexp(seq(0,2,by = 0.01), start / pfix)) %>%
  ggplot() + 
  geom_line(aes(x=x,y=y)) + 
  theme_bw() + 
  xlab("Starting Freq") + 
  ylab("Density")


