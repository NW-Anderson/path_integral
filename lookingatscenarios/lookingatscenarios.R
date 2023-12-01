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
dev.off()

setwd("~/Documents/GitHub/path_integral/lookingatscenarios/data")
list.files()
count = 0
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  seed <- as.numeric(strsplit(file,split = "p")[[1]][1])
  tmp$seed <- seed
  master <- dplyr::bind_rows(master, tmp)
}
# rm(tmp,file, seed)
df <- master
#%>% filter(gen >= 19900)

p1 <- ggplot(df, aes(x = gen, y = mean_pheno)) + 
  geom_line(aes(group = seed,
                color = scenario),
            size = 1, alpha = 1) 
  
p2 <- ggplot(df, aes(x = gen, y = mean_fit)) + 
  geom_line(aes(group = seed,
                color = scenario),
            size = 1, alpha = 1) 

p3 <- ggplot(df, aes(x = gen, y = var_pheno)) + 
  geom_line(aes(group = seed,
                color = scenario),
            size = 1, alpha = 1) 

p1 + p2 / p3 + plot_layout(guides = "collect")

