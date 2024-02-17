library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

############
## Linked ##
############
setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_positive_eff_09_11.csv.gz")  %>% mutate(clr = "1")

allele_freqs$group_id = case_when(allele_freqs$group_id == "U=0.0025_a=0.005" ~ "U = 0.0025 \u03b1 = 0.005",
                                  allele_freqs$group_id == "U=0.0025_a=0.01" ~ "U = 0.0025 \u03b1 = 0.01",
                                  allele_freqs$group_id == "U=0.025_a=0.005" ~ "U = 0.025 \u03b1 = 0.005",
                                  allele_freqs$group_id == "U=0.025_a=0.01" ~ "U = 0.025 \u03b1 = 0.01")

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("linked_sim_pint_densities.csv") %>% 
  mutate(group_id = case_when(Scenario == 1 ~ "U = 0.0025 \u03b1 = 0.005",
                              Scenario == 2 ~ "U = 0.0025 \u03b1 = 0.01",
                              Scenario == 3 ~ "U = 0.025 \u03b1 = 0.005",
                              Scenario == 4 ~ "U = 0.025 \u03b1 = 0.01"),
         
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
    title = "Linked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_wrap(vars(group_id)) + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
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
             color = turbo(11)[2]) + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) 

##############
## Unlinked ##
##############

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("unlinked_positive_eff_09_11.csv.gz")  %>% mutate(clr = "1")
allele_freqs$group_id = case_when(allele_freqs$group_id == "mu_0.0025.a_0.005" ~ "U = 0.0025 \u03b1 = 0.005",
                                  allele_freqs$group_id == "mu_0.0025.a_0.01" ~ "U = 0.0025 \u03b1 = 0.01",
                                  allele_freqs$group_id == "mu_0.025.a_0.005" ~ "U = 0.025 \u03b1 = 0.005",
                                  allele_freqs$group_id == "mu_0.025.a_0.01" ~ "U = 0.025 \u03b1 = 0.01")

labels <- unique(allele_freqs$group_id)

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("unlinked_sim_densities.csv.gz") %>% 
  mutate(group_id = case_when(scenario == 1 ~ labels[1],
                              scenario == 2 ~ labels[2],
                              scenario == 3 ~ labels[3],
                              scenario == 4 ~ labels[4]),
         
         clr = "2")

pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freqs <= 0.5, end_freqs > 0)

pint_means <- pints %>% group_by(group_id) %>% 
  mutate(total_dens = sum(density)) %>% 
  mutate(norm_dens = density / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq))

group_means <- allele_freqs %>% group_by(group_id) %>%
  summarize(mn = mean(end_freqs))

ggplot(data = allele_freqs, aes(x = end_freqs, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Unlinked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_wrap(vars(group_id)) + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
  geom_line(data = pints, aes(x = end_freq, y = density, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11], turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Path Integral")) + 
  theme(legend.position = "bottom") + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) +
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) 

#############
## Neutral ##
#############

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_neutral_09_11.csv.gz")  %>% mutate(clr = "1")

allele_freqs$par = case_when(allele_freqs$par == "U=0.0025_a=0.005" ~ "U = 0.0025 \u03b1 = 0.005",
                                  allele_freqs$par == "U=0.0025_a=0.01" ~ "U = 0.0025 \u03b1 = 0.01",
                                  allele_freqs$par == "U=0.025_a=0.005" ~ "U = 0.025 \u03b1 = 0.005",
                                  allele_freqs$par == "U=0.025_a=0.01" ~ "U = 0.025 \u03b1 = 0.01")

labels <- unique(allele_freqs$par)

setwd("/media/nathan/T7/path_integral/simulations/out/KimuraComparison")

pints <- fread("linked_sim_kimura_densities.csv") %>% 
  mutate(clr = "2")

setwd("/media/nathan/T7/path_integral/trueNeutralSims")

true_neut <- fread("trueNeutralMaster.csv.gz") %>%
  mutate(clr = "3")


pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freq <= 0.5, end_freq > 0)
true_neut <- true_neut %>% filter(freq <= 0.5, freq > 0)

pint_means <- pints %>% group_by(Scenario) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq))

group_means <- allele_freqs %>% group_by(par) %>%
  summarize(mn = mean(end_freq))

true_neut_means <- true_neut %>% summarize(mn = mean(freq))

ggplot(data = allele_freqs, aes(x = end_freq, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 0.75) +
  geom_vline(data = group_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(11)[11]) +
  xlim(0,1) + 
  labs(
    title = "Neutral: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_wrap(vars(par)) + 
  theme_bw() +
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) +
  geom_density(data = true_neut, aes(x = freq, color =clr),
            size = 1.25,
            alpha = 0.75,
            show.legend = F) +
  geom_vline(data = true_neut_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(11)[7]) +
  scale_color_manual(values = c(turbo(11)[11],  
                                turbo(11)[7], 
                                turbo(11)[2]),
                     name = "",
                     labels = c("Linked Selection\n Simulation",
                                "True Neutral\n Simulation",
                                "Kimura's Solution")) + 
  theme(legend.position = "bottom") + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) 

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

# p1 <- ggplot(pos_allele_freqs, aes(x = end, fill = group_id, color = group_id)) + 
#   geom_density(alpha = 0) + ggtitle("positive")
# 
# p2 <- ggplot(neg_allele_freqs, aes(x = end, fill = group_id, color = group_id)) + 
#   geom_density(alpha = 0) + ggtitle("negative")
# p1+p2



