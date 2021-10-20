source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_community_composition.R"))

library(withr)
library(phyloseq)
library(breakaway)
library(DivNet)

# Make phyloseq object ----------------------------------------------------

otumat <- counts %>%
  dplyr::select(sample, strainID, count.correct) %>%
  pivot_wider(names_from=sample, values_from=count.correct) %>%
  column_to_rownames(var="strainID") %>%
  as.matrix()

taxmat <- tax %>%
  column_to_rownames(var="strainID") %>%
  as.matrix()

metadf <- counts %>% 
  dplyr::select(sample=sample, microcosmID, treatment, day, cil, wrm) %>%
  distinct() %>%
  mutate(tpool=factor(paste(treatment, day, sep="_"))) %>%
  column_to_rownames(var="sample")

myphyseq <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE), 
                     tax_table(taxmat),
                     sample_data(metadf))

# Estimate diversity ------------------------------------------------------

# standard plugin estimate

plugin <- estimate_richness(myphyseq, measures = "Shannon") %>% 
  rownames_to_column(var="sample") %>%
  dplyr::rename(estimate=Shannon) %>%
  mutate(div="Shannon", inference="plugin")


# Divnet estimate
with_seed(24526, dv.cov1 <- divnet(myphyseq, ncores=1, tuning="fast",
                  X=model.matrix(~ tpool, data = metadf)))

# combine plugin and divnet
summary <- summary(dv.cov1$shannon) %>%
  dplyr::select(sample=sample_names, estimate, lower, upper) %>%
  mutate(div="Shannon", inference="divnet") %>%
  bind_rows(., plugin)%>%
  left_join(., dplyr::select(counts, sample, microcosmID, treatment, replicate, day, cil, wrm) %>% distinct())

# Save results ------------------------------------------------------------
write_rds(dv.cov1, here::here("data", "divnet_out.rds"))

write_rds(summary, here::here("data", "shannon_summary.rds"))