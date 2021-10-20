library(tidyverse)
library(here)
library(rcartocolor)
source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_community_composition.R"))

# Format data -------------------------------------------------------------

# Convert to relative abundances
counts_f <- counts %>%
  group_by(sample) %>%
  mutate(count1=sum(count.correct)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq=count.correct/count1) %>% select(-count1) %>%
  ungroup() %>%
  left_join(., tax)


# Collapse species with < 1% mean abundance into an "other" category

counts_f.99 <- counts_f %>%
  filter(day > 0) %>%
  group_by(strainID) %>%
  mutate(sum=sum(count.correct)) %>%
  distinct(strainID, sum) %>%
  ungroup() %>%
  mutate(sum1=sum(sum)) %>%
  mutate(rank = sum/sum1*100) %>%
  filter(rank >= 1) %>%
  pull(strainID) 


# make data frame for plotting
full <- counts_f %>%
  mutate(genus=ifelse(strainID %in% counts_f.99, genus, "Other")) %>%
  mutate(species=ifelse(strainID %in% counts_f.99, species, "Other")) %>%
  mutate(strainID=ifelse(strainID %in% counts_f.99, strainID, "Other")) %>%
  group_by(strainID, genus, species, day, microcosmID) %>%
  mutate(freq=sum(freq)) %>%
  ungroup() %>%
  distinct(freq, strainID, genus, species, day, microcosmID, .keep_all=TRUE) %>%
  mutate(strainID=factor(strainID, levels = c("Other", "HAMBI-0006", "HAMBI-0105", "HAMBI-1287",
                                              "HAMBI-1923", "HAMBI-1972", "HAMBI-2659",
                                              "HAMBI-0403", "HAMBI-1292", "HAMBI-1977")))

# Plot data ---------------------------------------------------------------


p <- ggplot(data = full, aes(x=day)) +
  geom_area(aes(y=freq, fill=strainID, group=interaction(replicate, strainID)),
            color="black", size=0.25) +
  facet_grid(replicate ~ treatment) + 
  scale_fill_manual(values=c("#A5AA99", "#E58606", "#5D69B1", "#52BCA3", 
                             "#99C945", "#CC61B0", "#24796C", "#DAA51B", 
                             "#2F8AC4", "#764E9F")) +
  labs(x="Day", y="Scaled relative abundance", fill="Strain") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

ggsave(here::here("figs", "figS1a.svg"), p, width=17.8, height=11.8, units="cm",
       device="svg")
