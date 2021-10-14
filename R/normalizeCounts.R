library(here)
library(tidyverse)
library(ggforce)
library(vegan)
library(DESeq2)
library(withr)

# HMSC contains only a couple of different probability distributions for the 
# response variable - "normal", "probit", "poisson", "log normal poisson." For 
# sequencing data the best option is usually something like negative binomial 
# or beta binomial which can account for extra dispersion and implicitly 
# account for sequencing depth (Num trials parameter). HMSC is not really 
# designed for sequencing data so we can't use these count models and instead
# need to do our best to transform to something that resembles a normal distribution.
# The easiest way is to put counts to the log scale.

# Load data ---------------------------------------------------------------
counts.corrected <-
  read_tsv(here("data", "correctedSpeciesCounts.tsv"), col_types = "ccccdd") %>%
  filter(sample != "HAMBI1")

expdesign <-
  read_tsv(here("data", "metadata.tsv"), col_types = "ccddc")

predator <-
  read_tsv(here("data", "predatorDensity.tsv"), col_types = "dccdddddd") %>%
  filter(replicate != 1)

counts.corrected %>%
  group_by(sample) %>%
  mutate(rab = count.correct / sum(count.correct))

# Formatting --------------------------------------------------------------
# Using rlog transform from DESeq2 package to correct all the samples 
# for different sequencing library size, to shrink variances in low-count 
# species, and to transform read count abundances to log scale for use in HMSC and LGPR.

# We need to include predator counts and filter out replicate 1 and day 0.
# We drop the first replicate because we have no bacteria 16S sequencing data 
# from that replicate - only predator densities. For the Day 0 samples we will 
# have to wrangle these separately by using the size factors estimated from 
# the full dataset.

meta1 <- left_join(expdesign, predator) %>%
  filter(day > 0) %>%
  dplyr::select(-date, -ciliate_count_raw, -worm_count_raw)

counts1 <- left_join(expdesign, counts.corrected) %>%
  filter(day > 0)

# Rlog transform ----------------------------------------------------------

# Make a DESeq2 formatted object: First format data frames for input
count.df <- counts1 %>%
  arrange(sample) %>%
  dplyr::select(sample, strainID, count.correct) %>%
  pivot_wider(names_from = sample, values_from = count.correct) %>%
  column_to_rownames(var = "strainID") %>%
  as.data.frame()

meta.df <- meta1 %>%
  arrange(sample) %>%
  mutate(day = factor(day)) %>%
  mutate(treatment = factor(treatment)) %>%
  column_to_rownames(var = "sample") %>%
  as.data.frame()

# Generate the DESeq2 dataset
# We specify a model of formula `~ trt`. This just makes sure variance across 
# conditions/times is not pooled. Doesn't actually matter for rlog transform 
# only neg binomial regression

dds <- DESeqDataSetFromMatrix(
  countData = count.df,
  colData = meta.df,
  design = ~ treatment + day + treatment:day,
  tidy = FALSE
)

## Perform the rlog transformation
with_preserve_seed(rld <- rlog(dds, blind = FALSE))


### NOTE: Blind dispersion estimation
# The two functions, vst and rlog have an argument blind, for whether the 
# transformation should be blind to the sample information specified by the 
# design formula. When blind equals TRUE (the default), the functions will 
# re-estimate the dispersions using only an intercept. This setting should be 
# used in order to compare samples in a manner wholly unbiased by the 
# information about experimental groups, for example to perform sample QA 
# (quality assurance) as demonstrated below.

# However, blind dispersion estimation is not the appropriate choice if one 
# expects that many or the majority of genes (rows) will have large differences 
# in counts which are explainable by the experimental design, and one wishes to 
# transform the data for downstream analysis. In this case, using blind 
# dispersion estimation will lead to large estimates of dispersion, as it 
# attributes differences due to experimental design as unwanted noise, and will 
# result in overly shrinking the transformed values towards each other. 
# By setting blind to FALSE, the dispersions already estimated will be used to 
# perform transformations, or if not present, they will be estimated using the 
# current design formula. Note that only the fitted dispersion estimates from 
# mean-dispersion trend line are used in the transformation (the global 
# dependence of dispersion on mean for the entire experiment). So setting blind 
# to FALSE is still for the most part not using the information about which 
# samples were in which experimental group in applying the transformation.

# Convert to transformed assay
rldass <- assay(rld)

rldtb <- rldass %>%
  as_tibble(rownames = "strainID") %>%
  pivot_longer(-strainID, names_to = "sample", values_to = "log2fc") %>%
  mutate(logcount = log(2 ^ log2fc))

## Add control samples back
count.df.ctrl <- left_join(expdesign, counts.corrected) %>%
  arrange(sample) %>%
  dplyr::select(sample, strainID, count.correct) %>%
  distinct() %>%
  pivot_wider(names_from = sample, values_from = count.correct) %>%
  column_to_rownames(var = "strainID") %>%
  as.data.frame()

meta.df.ctrl <-  left_join(expdesign, predator) %>%
  mutate(
    treatment = ifelse(str_detect(sample, "HAMBI"), "ctrl", treatment),
    description = ifelse(str_detect(sample, "HAMBI"), "controlT0", description),
    day = ifelse(str_detect(sample, "HAMBI"), 65, day)
  ) %>%
  distinct() %>%
  arrange(sample) %>%
  mutate(day = factor(day)) %>%
  mutate(treatment = factor(treatment)) %>%
  mutate(trt = factor(paste(treatment, day, sep = "."))) %>%
  select(sample, trt) %>%
  column_to_rownames(var = "sample") %>%
  as.data.frame()

# Make new deseq2 object
dds.ctrl <- DESeqDataSetFromMatrix(
  countData = count.df.ctrl,
  colData = meta.df.ctrl,
  design = ~ trt,
  tidy = FALSE
)

# calculate size factors for T0 samples
sizefactors <-
  as_tibble(sizeFactors(estimateSizeFactors(dds.ctrl)), rownames = "sample") %>%
  dplyr::rename(sf = value)

ctrllogcounts <- left_join(sizefactors, counts.corrected) %>%
  mutate(sflogcount = log((count.correct / sf) + 1)) %>%
  dplyr::select(sample, strainID, sflogcount)

# Join control and experimental samples
# Rlog transform has a difficult time with some low abundance strains and 
# gives very high values (like exp(80)). Here we set really high/low rlog 
# counts to size factor scaled values from DESeq2 if the difference is larger 
# than 2 log units. To remove any zero inflation we set any remaining zeros to 
# a random value centered on the mean of the lowest log abundance (-2).

norm1 <- left_join(ctrllogcounts, rldtb) %>%
  left_join(expdesign) %>%
  mutate(sample1 = ifelse(
    str_detect(sample, "HAMBI"),
    paste(treatment,
          paste(replicate,
                paste("d", day, sep = ""), sep = "_"), sep =
            ""),
    sample
  )) %>%
  left_join(counts.corrected) %>%
  dplyr::select(
    sample,
    sample1,
    strainID,
    treatment,
    day,
    replicate,
    sflogcount,
    rlogcount = logcount,
    count.correct,
    count
  ) %>%
  mutate(
    rlogcount = ifelse(is.na(rlogcount), sflogcount, rlogcount),
    diff = abs(rlogcount - sflogcount)
  ) %>%
  mutate(rlogcount = ifelse(diff > 3, sflogcount, rlogcount)) %>%
  mutate(rlogcount = ifelse(rlogcount == 0,
                            rnorm(100, mean = -2, sd = 0.5), rlogcount))

# Other normalization -----------------------------------------------------

# We compare the DESeq2 normalization with two other approaches.
# 1. rarefying to minimum library depth. Minimum library size = 6098 reads
# 2. Using log relative abundance
# 3. Presence/absence for use with probit models

with_preserve_seed(rarefied <- counts.corrected %>%
  dplyr::select(-count) %>%
  dplyr::select(-genus,-species) %>%
  pivot_wider(names_from = "strainID", values_from = "count.correct") %>%
  column_to_rownames(var = "sample") %>%
  as.data.frame() %>%
  rrarefy(., sample = 6098) %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(`HAMBI-0006`:`HAMBI-3237`,
               names_to = "strainID",
               values_to = "rarecount") %>%
  mutate(lograrecount = log(rarecount + 1)))

RA <- counts.corrected %>%
  mutate(PA = ifelse(count.correct > 0, 1, 0)) %>%
  group_by(sample) %>%
  mutate(log10rab = log10((count.correct + 1) / sum(count.correct + 1))) %>%
  ungroup() %>%
  dplyr::select(sample, strainID, PA, log10rab)

norm2 <- left_join(norm1, RA) %>%
  left_join(., rarefied) %>%
  dplyr::select(-sample) %>%
  dplyr::rename(sample = sample1) %>%
  mutate(
    cil = factor(case_when(
      str_detect(treatment, "Pevo") ~ 2,
      str_detect(treatment, "Panc") ~ 1,
      TRUE ~ 0
    )),
    wrm = factor(case_when(str_detect(treatment, "HN") ~ 1,
                           TRUE ~ 0)),
    replicate = factor(
      case_when(
        replicate == 2 ~ "A",
        replicate == 3 ~ "B",
        replicate == 4 ~ "C",
        replicate == 5 ~ "D"
      )
    )
  ) %>%
  mutate(
    microcosmID = factor(paste(treatment, replicate, sep = "_")),
    treatment = factor(
      treatment,
      levels = c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")
    ),
    sample = factor(sample),
    strainID = factor(strainID)
  ) %>%
  dplyr::select(
    sample,
    microcosmID,
    treatment,
    replicate,
    day,
    cil,
    wrm,
    strainID,
    count,
    count.correct,
    lograrecount,
    log10rab,
    sflogcount,
    rlogcount,
    PA
  ) %>%
  arrange(microcosmID, day)

write_tsv(norm2, here::here("data", "normalizedCorrectedSpeciesCounts.tsv"))
saveRDS(norm2, here::here("data", "normalizedCorrectedSpeciesCounts.rds"))

# Plotting ----------------------------------------------------------------
norm2 %>%
  pivot_longer(
    c(count.correct, lograrecount, log10rab, rlogcount),
    names_to = "metric",
    values_to = "val"
  ) %>% 
  ggplot() + 
  geom_point(aes(x=day, y=val, color=treatment)) + 
  geom_line(aes(x=day, y=val, color=treatment, group=microcosmID)) + 
  facet_wrap_paginate(strainID ~ metric, ncol=4, nrow=4, page=4, scales="free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#sessioninfo::session_info()