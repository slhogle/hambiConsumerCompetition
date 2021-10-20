source(here::here("R", "utils_generic.R"))

library(ape)
library(Hmsc)
library(withr)

# Format quasi equilibrium, phase and prepare models for HMSC on puhti cluster

# Read data ---------------------------------------------------------------

counts <- readr::read_rds(here::here("data", "normalized_corrected_species_counts.rds")) %>%
  mutate(cil=factor(cil), wrm=factor(wrm)) %>%
  group_by(sample) %>%
  mutate(effort=log(sum(count.correct))) %>%
  ungroup() %>%
  mutate(Y=log(count.correct+1)) %>%
  mutate(Y=ifelse(Y==0, NA_real_, Y))


# Prepare abundance and covariate data ------------------------------------

# prepare Y, X, SD matrices
# Includes days 5, 9, 13. 
# This is the sorting phase identified by the change point analysis


# These species are too rare to model accurately
rm_sort <- c("HAMBI-0262")

counts_sort <- counts %>%
  filter(day < 17 & day != 0) %>%
  filter(!(strainID %in% rm_sort)) %>%
  mutate(strainID=factor(strainID)) 

# Y matrix (these are species abundances)

# abundance conditional upon presence (COP)
# Ya is an abundance matrix but we set 0 values == NA
Ya <- counts_sort %>%
  select(sample, strainID, Y=Y) %>%
  pivot_wider(names_from="strainID", values_from="Y") %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

# presence/absence
# Yp in the richness matrix of just 0/1
Yp <- counts_sort %>%
  select(sample, strainID, Y=PA) %>%
  pivot_wider(names_from="strainID", values_from="Y") %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

# X matrix
XData <- counts_sort %>%
  select(sample, day, cil, wrm, effort) %>% #replicate, treatment
  distinct() %>%
  column_to_rownames(var="sample") %>%
  as.data.frame()

# study design
studyDesign <- counts_sort %>% 
  select(sample, day.coord=day, microcosmID.coord=microcosmID) %>%
  distinct() %>%
  mutate(microcosmID.unit=factor(paste0("microcosmID_", microcosmID.coord)),
         day.unit=factor(paste0("day_",day.coord))) %>%
  column_to_rownames(var='sample') %>%
  as.data.frame()

# trait data
TrData <- read_tsv(here::here("data", "HAMBI-traits.tsv")) %>%
  select(strain, D, n_carbon_ecoplate, r_specific_growth_rate, biofilm_formation) %>%
  arrange(strain) %>%
  column_to_rownames(var="strain") %>% 
  mutate(across(everything(), scale2)) %>% 
  as.data.frame()

TrData <- TrData[!(row.names(TrData) %in% rm_sort), ]

# phylogenetic tree
phyloTree <- read.tree(here::here("data", "phylotree.nwk"))
phyloTree <- drop.tip(phyloTree, rm_sort)

# random effects
time <- as.matrix(data.frame(day=unique(studyDesign$day.coord)))
rownames(time) <- unique(studyDesign$day.unit)
rLday <- HmscRandomLevel(sData = time)

microcosmID <- as.matrix(data.frame(microcosmID=unique(studyDesign$microcosmID.coord)))
rownames(microcosmID) <- unique(studyDesign$microcosmID.unit)
rLmicrocosmID <- HmscRandomLevel(units = microcosmID)

studyDesign$day <- studyDesign$day.unit
studyDesign$microcosmID <- studyDesign$microcosmID.unit
studyDesign <- studyDesign[c(5,6)]

# Make model objects ------------------------------------------------------

# Here we will fit three different types of models.

# 1. A "FULL" model that includes both predation regime (environmental) and 
# sequencing effort as fixed covariates as well as temporal random effects 
# and hierarchical/repeated measure random effects.

# 2. A "ENVI" model that includes predation regime (environmental) and 
# sequencing effort as fixed covariates but does not account for random effects.

# 3. A "TIME" model that excludes environmental covariates 
# (sequencing effort only) and only models baseline abundance, but 
# includes temporal and hierarchical/repeated measure random effects.

# Hurdle models -----------------------------------------------------------

# Ecological count data are often dominated by zeros, i.e. they are 
# zero-inflated. One standard model that can be fitted to such data is 
# the Zero-Inflated Poisson model. While HMSC-R 3.0 does not include this or 
# other zero-inflated models, we note that it is always possible to apply the 
# closely-related Hurdle approach. In the Hurdle model, one analyses separately 
# occurrence data 1*(Y>0) and abundance conditional on presence Y[Y==0] = NA. 
# We note that these two aspects of the data (occurrence and abundance 
# conditional on presence) are statistically independent of each other, and 
# thus it is interesting to compare results obtained for them, e.g. to see if 
# variation in occurrence and variation in abundance are explained by the same 
# environmental covariates. In HMSC-R 3.0, occurrence data are analysed by the
# probit model, whereas the best model for abundance (conditional on presence) 
# depends on the type of the data.

# Hurdle model: abundance -------------------------------------------------
ma_full <- Hmsc(Y           = Ya, 
                XData       = XData, 
                XFormula    = ~day*cil*wrm + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                ranLevels   = list(day=rLday, microcosmID=rLmicrocosmID),
                distr       = "normal",
                YScale      = TRUE)

ma_envi <- Hmsc(Y           = Ya, 
                XData       = XData, 
                XFormula    = ~day*cil*wrm + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                distr       = "normal",
                YScale      = TRUE)

ma_time <- Hmsc(Y           = Ya, 
                XData       = XData, 
                XFormula    = ~day + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                distr       = "normal",
                YScale      = TRUE)

saveRDS(ma_full, here::here("data", "JSDM_unfit", "sort_ma_full.rds"))
saveRDS(ma_envi, here::here("data", "JSDM_unfit", "sort_ma_envi.rds"))
saveRDS(ma_time, here::here("data", "JSDM_unfit", "sort_ma_time.rds"))

# Hurdle model: probit ----------------------------------------------------
mp_full <- Hmsc(Y           = Yp, 
                XData       = XData, 
                XFormula    = ~day*cil*wrm + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                ranLevels   = list(day=rLday, microcosmID=rLmicrocosmID),
                distr       = "probit",
                YScale      = TRUE)

mp_envi <- Hmsc(Y           = Yp, 
                XData       = XData, 
                XFormula    = ~day*cil*wrm + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                distr       = "probit",
                YScale      = TRUE)

mp_time <- Hmsc(Y           = Yp, 
                XData       = XData, 
                XFormula    = ~day + effort, 
                TrData      = TrData,
                TrFormula   = formula(paste0("~",paste(names(TrData), collapse = "+"))),
                phyloTree   = phyloTree,
                studyDesign = studyDesign,
                distr       = "probit",
                YScale      = TRUE)

saveRDS(mp_full, here::here("data", "JSDM_unfit","sort_mp_full.rds"))
saveRDS(mp_envi, here::here("data", "JSDM_unfit","sort_mp_envi.rds"))
saveRDS(mp_time, here::here("data", "JSDM_unfit","sort_mp_time.rds"))

# Test run HMSC -----------------------------------------------------------

# run only a couple samples to make sure models are correct format

m1 <- sampleMcmc(mp_time, thin = 1, samples = 100, 
                 transient = 50, nChains = 2)
