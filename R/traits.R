source(here::here("r", "utils_generic.R"))

library(here)
library(tidyverse)
library(ggfortify)
library(patchwork)
library(rcartocolor)
library(scales)
library(withr)

# Read data ---------------------------------------------------------------

# predator <- read_rds(here::here("data", "formatted_predator_prey_density.rds")) %>%
#   mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
#   mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
#          microcosmID=factor(microcosmID))
# 
# counts <- read_rds(here::here("data", "normalized_corrected_species_counts.rds"))

traits <- read_tsv(here::here("data", "HAMBI-traits.tsv"))


# Plot scaled traits ------------------------------------------------------

p <- scale(column_to_rownames(traits, var = "strain")) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "strainID") %>%
  pivot_longer(cols=-strainID) %>%
  filter(!(name %in% c("auc_carrying_capacity", "auc_pred", "ncarbon_total"))) %>%
  mutate(name=factor(name, levels=c("D", "biofilm_formation",  "n_carbon_ecoplate", 
                                    "r_specific_growth_rate", "k_carrying_capacity", "g_doubling_time",
                                    "motile", "facultative", 
                                    "growth_carboxylic_acid", "growth_sugar", "growth_peptide",
                                    "glucose_fermenter"))) %>%
  ggplot(aes(x=name, y=strainID)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(name="Scaled\nTrait Value", low = "#009392", high="#cf597e") +
  labs(x="", y="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks=element_blank(),
        plot.background=element_blank(),
        rect = element_blank())

# Trait redundancy --------------------------------------------------------

# Meed to reduce redundancy of traits so that models are most effective. 
# Many of the traits are colinear/correlated. So including them all in the
# model makes estimating the their coefficients for their effect on species 
# niches difficult. Strategy here is to use PCA to reduce traits down to 
# a handful of dimensions.

# PCA

withr::with_seed(12378, pca_input 
                 <- scale(column_to_rownames(traits, var = "strain")))

withr::with_seed(12378, pca_res <- prcomp(pca_input))

# Using the first 6 components accounts for 
# about 90% of the variance in the trait data
summary(pca_res)


# PCA loadings

# Defense traits and high doubling time mostly on PC1. 
# Growth traits mostly on PC2

# plot loadings
load <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column(var="trait") %>%
  as_tibble() %>%
  pivot_longer(PC1:PC15, names_to="PC") %>%
  mutate(PC=factor(PC, levels=c("PC1", "PC2", "PC3", "PC4",
                                "PC5", "PC6", "PC7", "PC8",
                                "PC9", "PC10", "PC11", "PC12",
                                "PC13", "PC14", "PC15")))

pload <- ggplot(load, aes(y=trait, x=PC, fill=value)) +
  geom_tile() +
  labs(y="", x="Principal Component") +
  scale_fill_carto_c(palette = "Tropic") + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=here::here("figs", "traits_PCA_loadings.svg"), plot=pload, device="svg",
       units="cm", height=14, width=12)

pload

# PCA species scores
TrData <- as.data.frame(pca_res$x[,1:6])

# save loadings for later (just in case)
readr::write_rds(TrData, here::here("data", "TrData.rds"))

# plot PCA
ppca <- autoplot(pca_res, 
                 x=1,
                 y=2,
                 loadings = TRUE, loadings.label = TRUE, label.repel=TRUE, 
                 loadings.label.repel=TRUE, color="strain",
                 data = traits) + 
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())


#final plot

pfinal <- (ppca + p) + plot_layout(guides = 'collect')
pfinal

ggsave(filename=here::here("figs", "figS4.svg"), plot=pfinal, device="svg",
       units="cm", height=12, width=22.8)
