library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

################
#### Linked ####
################
setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_positive_eff_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(group_id %in% c("U=0.0025_a=0.005", "U=0.0025_a=0.01") ~ "U = 0.0025",
                          group_id %in% c("U=0.025_a=0.005", "U=0.025_a=0.01") ~ "U = 0.025"),
         selCoef = case_when(group_id %in% c("U=0.025_a=0.01", "U=0.0025_a=0.01") ~ "\u03b1 = 0.01",
                             group_id %in% c("U=0.025_a=0.005", "U=0.0025_a=0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("linked_sim_pint_densities.csv") %>% 
  mutate(clr = "2",
         bigU = case_when(Scenario %in% c(1,2) ~ "U = 0.0025",
                          Scenario %in% c(3,4) ~ "U = 0.025"),
         selCoef = case_when(Scenario %in% c(1,3) ~ "\u03b1 = 0.005",
                             Scenario %in% c(2,4) ~ "\u03b1 = 0.01"))

setwd("/media/nathan/T7/path_integral/genicSimComparison")

genic <- data.frame()
for(file in list.files()){
  params <- strsplit(file, split = "_")[[1]] 
  end <- params[2] %>% 
    gsub("end", "", .) %>% 
    as.numeric()
  selCoef <- params[3] %>% 
    gsub("selCoef.csv", "", .) %>% 
    as.numeric()
  
  tmp <- fread(file) %>%
    unlist()
  genic <- dplyr::bind_rows(genic, data.frame(dens = tmp,
                                              end = end, 
                                              selCoef = selCoef))
}

genic <- genic %>% 
  mutate(selCoef = case_when(selCoef == 5 ~ "\u03b1 = 0.005",
                             selCoef == 10 ~ "\u03b1 = 0.01"),
         clr = "3")

# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end <= 0.5, end > 0)
# genic <- genic %>% filter(end <= 0.5, end > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end > 0)
genic <- genic %>% filter(end > 0)

pint_means <- pints %>% group_by(bigU, selCoef) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq)) %>%
  ungroup()

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
  summarize(mn = mean(end)) %>%
  ungroup()

genic_means <- genic %>% group_by(selCoef) %>%
  mutate(total_dens = sum(dens)) %>%
  mutate(norm_dens = dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end)) %>%
  ungroup()

ggplot(data = allele_freqs, aes(x = end, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Linked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_grid(rows = vars(bigU),
             cols = vars(selCoef)) + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11],  
                                turbo(11)[7], 
                                turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Hayward",
                                "Genic")) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[7]) + 
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") 

########
# main #
########

allele_freqs <- allele_freqs %>% 
  filter(group_id == "U=0.025_a=0.01")

pints <- pints %>% filter(Scenario == 4)

genic <- genic %>% filter(selCoef == "α = 0.01")
genic_means <-  genic_means %>% filter(selCoef == "α = 0.01")
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")
pint_means <- pint_means %>% filter(bigU == "U = 0.025",
                                    selCoef == "α = 0.01")

p1 <- ggplot(data = allele_freqs, aes(x = end, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Linked") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11],  
                                turbo(11)[7], 
                                turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Hayward",
                                "Genic")) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[7]) + 
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") 

p1

##################
#### Unlinked ####
##################

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("unlinked_positive_eff_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(group_id %in% c("mu_0.0025.a_0.005", "mu_0.0025.a_0.01") ~ "U = 0.0025",
                          group_id %in% c("mu_0.025.a_0.005", "mu_0.025.a_0.01") ~ "U = 0.025"),
         selCoef = case_when(group_id %in% c("mu_0.025.a_0.01", "mu_0.0025.a_0.01") ~ "\u03b1 = 0.01",
                             group_id %in% c("mu_0.025.a_0.005", "mu_0.0025.a_0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("unlinked_sim_densities.csv.gz") %>% 
  mutate(clr = "2",
         bigU = case_when(scenario %in% c(1,2) ~ "U = 0.0025",
                          scenario %in% c(3,4) ~ "U = 0.025"),
         selCoef = case_when(scenario %in% c(1,3) ~ "\u03b1 = 0.005",
                             scenario %in% c(2,4) ~ "\u03b1 = 0.01"))

setwd("/media/nathan/T7/path_integral/genicSimComparison")

genic <- data.frame()
for(file in list.files()){
  params <- strsplit(file, split = "_")[[1]] 
  end <- params[2] %>% 
    gsub("end", "", .) %>% 
    as.numeric()
  selCoef <- params[3] %>% 
    gsub("selCoef.csv", "", .) %>% 
    as.numeric()
  
  tmp <- fread(file) %>%
    unlist()
  genic <- dplyr::bind_rows(genic, data.frame(dens = tmp,
                                              end = end, 
                                              selCoef = selCoef))
}

genic <- genic %>% 
  mutate(selCoef = case_when(selCoef == 5 ~ "\u03b1 = 0.005",
                             selCoef == 10 ~ "\u03b1 = 0.01"),
         clr = "3")

# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end_freqs <= 0.5, end_freqs > 0)
# genic <- genic %>% filter(end <= 0.5, end > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freqs > 0)
genic <- genic %>% filter(end > 0)

pint_means <- pints %>% group_by(bigU, selCoef) %>% 
  mutate(total_dens = sum(density)) %>% 
  mutate(norm_dens = density / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq)) %>% 
  ungroup()

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
  summarize(mn = mean(end_freqs))

genic_means <- genic %>% group_by(selCoef) %>%
  mutate(total_dens = sum(dens)) %>%
  mutate(norm_dens = dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end)) %>%
  ungroup()

ggplot(data = allele_freqs, aes(x = end_freqs, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Unlinked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU)) + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
  geom_line(data = pints, aes(x = end_freq, y = density, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11],  
                                turbo(11)[7], 
                                turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Hayward",
                                "Genic")) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[7]) + 
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) +
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") 

########
# main #
########


allele_freqs <- allele_freqs %>% 
  filter(group_id == "mu_0.025.a_0.01")

pints <- pints %>% filter(scenario == 4)

genic <- genic %>% filter(selCoef == "α = 0.01")
genic_means <-  genic_means %>% filter(selCoef == "α = 0.01")
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")
pint_means <- pint_means %>% filter(bigU == "U = 0.025",
                                    selCoef == "α = 0.01")

p2 <- ggplot(data = allele_freqs, aes(x = end_freqs, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 1.75) + 
  xlim(0,1) + 
  labs(
    title = "Unlinked") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  theme_bw() +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(11)[11]) + 
  geom_line(data = pints, aes(x = end_freq, y = density, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_color_manual(values = c(turbo(11)[11],  
                                turbo(11)[7], 
                                turbo(11)[2]),
                     name = "",
                     labels = c("Simulation",
                                "Hayward",
                                "Genic")) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[7]) + 
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[2]) +
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") 

p2

# p1 + p2 + plot_layout(guides = "collect") + 
#   plot_annotation(tag_levels = 'A')  & 
#   theme(plot.tag = element_text(size = 24),
#         legend.position = "bottom")

#############
## Neutral ##
#############

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_neutral_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(par %in% c("U=0.0025_a=0.005", "U=0.0025_a=0.01") ~ "U = 0.0025",
                          par %in% c("U=0.025_a=0.005", "U=0.025_a=0.01") ~ "U = 0.025"),
         selCoef = case_when(par %in% c("U=0.025_a=0.01", "U=0.0025_a=0.01") ~ "\u03b1 = 0.01",
                             par %in% c("U=0.025_a=0.005", "U=0.0025_a=0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/KimuraComparison")

pints <- fread("linked_sim_kimura_densities.csv") %>% 
  mutate(clr = "2")

setwd("/media/nathan/T7/path_integral/trueNeutralSims")

true_neut <- fread("trueNeutralMaster.csv.gz") %>%
  mutate(clr = "3")


# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end_freq <= 0.5, end_freq > 0)
# true_neut <- true_neut %>% filter(freq <= 0.5, freq > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freq > 0)
true_neut <- true_neut %>% filter(freq > 0)

pint_means <- pints %>% group_by(Scenario) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq))

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
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
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU)) + 
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

########
# main #
########

allele_freqs <- allele_freqs %>% filter(par == "U=0.025_a=0.01")
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")

p3 <- ggplot(data = allele_freqs, aes(x = end_freq, color = clr)) + 
  geom_density(show.legend = F, size = 1.25,
               alpha = 0.75) +
  geom_vline(data = group_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(11)[10]) +
  xlim(0,1) + 
  labs(
    title = "Neutral") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  theme_bw() +
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(11)[1]) +
  geom_density(data = true_neut, aes(x = freq, color =clr),
               size = 1.25,
               alpha = 0.75,
               show.legend = F) +
  geom_vline(data = true_neut_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(11)[6]) +
  scale_color_manual(values = c(turbo(11)[10],  
                                turbo(11)[6], 
                                turbo(11)[1]),
                     name = "",
                     labels = c("Neutral + Linked\n Selection Simulation",
                                "Neutral\n Simulation",
                                "Kimura's Solution")) + 
  theme(legend.position = "bottom") + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) 

p3

p4 <- data.frame(x = 1,
                       y = 1,
                       txt = "LD FIG") %>%
  ggplot(., aes(x = x, y = y)) + 
  geom_text(aes(label = txt),
            size = 24) + 
  theme_bw() + 
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

(p1 + p2) / (p4 + p3) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 24),
        legend.position = "bottom")







