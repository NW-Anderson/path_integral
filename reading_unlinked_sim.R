library(data.table)
library(dplyr)
library(ggplot2)

df <- fread("e_and_r.mu_0.0025.a_0.005.txt")
df <- df %>% filter(x0 < 0.11,
             x0 > 0.09,
             a > 0)
end_freqs <- c()
for(i in 5:104){
  end_freqs <- c(end_freqs, df[,..i][[1]])
}
end_freqs <- data.frame(end_freqs)
ggplot(end_freqs, aes(x = end_freqs)) + geom_density()
