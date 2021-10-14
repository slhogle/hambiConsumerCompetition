library(here)
library(tidyverse)
library(metacal)
library(Polychrome)

# Read data ---------------------------------------------------------------
counts <- read_tsv(here("data", "speciesCounts.tsv"), col_types="ccccd")
expdesign <- read_tsv(here("data", "metadata.tsv"), col_types="ccddc")

# Prepare -----------------------------------------------------------------
# We will use the metacal package for estimating bias and performing 
# calibration in the special case where the bias of all the taxa of interest 
# can be directly measured from the control sample. Since we know the exact 
# concentration of cells added in the begining of the experiment we can 
# estimate a corrected/unbiased relative abundance.

# These are the control samples
controls <- c("HAMBI1","HAMBI2", "HAMBI3", "HAMBI4", "HAMBI5")

# Make a matrix of observed counts
observed_mat <- counts %>%
  filter(sample %in% controls) %>%
  select(sample, strainID, count) %>% 
  pivot_wider(names_from="strainID", values_from="count") %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

# Make a matrix of true proportions
actual_mat <- counts %>%
  filter(sample %in% controls) %>%
  mutate(prop=round(1/24, 3)) %>%
  select(sample, strainID, prop) %>% 
  pivot_wider(names_from="strainID", values_from="prop") %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

# Bias estimation ---------------------------------------------------------
mc_fit <- estimate_bias(observed_mat, actual_mat, 1, boot = TRUE) %>% print()
bias <- coef(mc_fit) %>% print

mc_fit.summary <- summary(mc_fit)
coef_tb <- mc_fit.summary$coefficients
print(mc_fit.summary)

# Plot bias estimation
coef_tb %>%
  mutate(taxon = fct_reorder(taxon, estimate)) %>%
  ggplot(aes(taxon, estimate, 
             ymin = estimate / gm_se^2, ymax = estimate * gm_se^2)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_pointrange() +
  scale_y_log10() +
  coord_flip() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot model fit
mypal <- unname(createPalette(24, c("#F3874AFF", "#FCD125FF"), M=5000))

observed.fitted <- fitted(mc_fit)

a <- as_tibble(observed.fitted) %>% 
  pivot_longer(everything()) %>% 
  mutate(type="Fitted")
b <- as_tibble(actual_mat) %>% 
  pivot_longer(everything()) %>% 
  mutate(type="Actual")
c <- as_tibble(observed_mat/rowSums(observed_mat)) %>% 
  pivot_longer(everything(), values_to="observed")

bind_rows(a,b) %>%
  left_join(., c) %>%
  ggplot(aes(x=observed, y=value, color = name)) +
  geom_abline(color = "darkgrey") +
  geom_point(show.legend = T) +
  scale_color_manual(values=mypal) +
  facet_wrap(~type) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# HAMBI-1842 is off by so much because it actually has 4 copies of the 16S
# rRNA gene!! Metacal procedure accounts for this! 

# Calibrate ---------------------------------------------------------------

# Make a matrix of observed counts
ps.mat <- counts %>%
  select(sample, strainID, count) %>% 
  pivot_wider(names_from="strainID", values_from="count") %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

ps.cal <- calibrate(ps.mat, bias, margin=1)

# convert to long
counts_final <- rownames_to_column(as.data.frame(ps.cal), var="sample") %>% 
  as_tibble() %>%
  pivot_longer(`HAMBI-0006`:last_col(), names_to="strainID", values_to="f.correct") %>%
  left_join(., counts) %>%
  group_by(sample) %>%
  mutate(total=sum(count)) %>%
  ungroup() %>%
  mutate(count.correct=total*f.correct,
         f.obs=count/total)

samples <- c("HAMBI1", "HN2-d5", "HN2-d9", "HN2-d13")

counts_final %>%
  filter(sample %in% samples) %>%
  dplyr::select(sample, strainID, f.correct, f.obs) %>%
  pivot_longer(cols=c(f.correct, f.obs), names_to="type", values_to="prop") %>%
  ggplot(aes(x=type, prop, fill = strainID)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values=mypal) +
  facet_wrap(~sample, ncol = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

# Output ------------------------------------------------------------------
# What I am doing here is multiplying each species relative abundance by the 
# total number of reads per sample so you get corrected counts. Counts are 
# necessary for DEseq normalization etc...

# Also make sure to remove the negative contols
counts_final %>%
  mutate(count.correct = round(count.correct)) %>%
  filter(!str_detect(sample, "negCTRL")) %>%
  dplyr::select(sample, strainID, genus, species, count, count.correct) %>%
  write_tsv(., here("processed_data", "correctedSpeciesCounts.tsv"))

#sessioninfo::session_info()