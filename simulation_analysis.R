library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_positive_eff_09_11.csv.gz")  %>% mutate(clr = "1")

labels <- unique(allele_freqs$group_id)

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("linked_sim_pint_densities.csv") %>% 
  mutate(group_id = case_when(Scenario == 1 ~ "U=0.0025_a=0.005",
                              Scenario == 2 ~ "U=0.0025_a=0.01",
                              Scenario == 3 ~ "U=0.025_a=0.005",
                              Scenario == 4 ~ "U=0.025_a=0.01",),
         
         clr = "2")

pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end <= 0.5, end > 0)

pint_means <- pints %>% group_by(group_id) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq))

group_means <- allele_freqs %>% group_by(group_id) %>%
  summarize(mn = mean(end))

ggplot(data = allele_freqs, aes(x = end, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Ending Frequency for Alleles Starting Between 0.09 and 0.10") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_wrap(vars(group_id)) + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[2]) + 
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11], turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Path Integral")) + 
  theme(legend.position = "bottom") + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[11])



# ggplot(data = pintDf, aes(x = popalpha, y = pintDetected)) + 
#   geom_line(aes(color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5) + 
#   geom_point(aes(color = as.factor(VG)),
#              alpha=1,
#              size = 2.5) +
#   geom_hline(yintercept=0.01, 
#              linetype="dashed", 
#              color = turbo(11)[11], size  = 0.75) + 
#   geom_vline(xintercept=1, 
#              linetype="dashed", 
#              color = turbo(11)[11], size  = 0.75) + 
#   theme_bw() +
#   guides(color=guide_legend(title = "Genetic Variance",
#                             override.aes = list(alpha=1))) +
#   scale_x_continuous("2 Ne \u03b1 \u039b / W", 
#                      breaks = c(0, 1, 5, 10, 15, 20), 
#                      labels = c(0, 1, 5, 10, 15, 20)) +
#   scale_y_continuous("P(detected)") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
#                                    colour = c(rep("black",1),"red",rep("black",4)))) + 
#   scale_color_viridis(discrete = T) + 
#   theme(axis.title = element_text(size=18)) + 
#   geom_point(data = numDf, aes(x = popalpha, 
#                                y = yval,
#                                color = as.factor(VG)),
#              alpha=c(0,1),
#              size = 2.5,
#              shape = 17) +
#   geom_line(data = numDf, aes(x = popalpha, 
#                               y = yval,
#                               color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5,
#             linetype = "dashed") 










# allele_freqs <- allele_freqs[start > .09 & start < .11 ]
# gc() 

# params <- fread(
#   "~/Documents/GitHub/path_integral/simulations/params.txt",
#   col.names = c("seed", "mut_rate", "effect_size")
# )
# setorder(params, mut_rate, effect_size)
# params[, group_id := .GRP, by = .(mut_rate, effect_size)]
# params[, group_id := as.factor(.SD$group_id)]


# allele_freqs <-allele_freqs[params, on="seed"]
# 
# allele_freqs <- allele_freqs %>% select(-"effect_size")
# effect_sizes <- fread("effect_sizes.csv")
# allele_freqs <- merge(allele_freqs, effect_sizes, by = c("seed", "site"))
# rm(effect_sizes)

# pos_allele_freqs <- allele_freqs[effect_size > 0]
# neg_allele_freqs <- allele_freqs[effect_size < 0]
# rm(allele_freqs)

p1 <- ggplot(pos_allele_freqs, aes(x = end, fill = group_id, color = group_id)) + 
  geom_density(alpha = 0) + ggtitle("positive")

p2 <- ggplot(neg_allele_freqs, aes(x = end, fill = group_id, color = group_id)) + 
  geom_density(alpha = 0) + ggtitle("negative")
p1+p2



