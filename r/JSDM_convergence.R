source(here::here("r", "utils_generic.R"))
source(here::here("r", "utils_JSDM.R"))

library(ape)
library(Hmsc)
library(withr)

# Sample names ------------------------------------------------------------

mynames <- crossing(
  n = c(100, 1000),
  p = c("sort", "qeq"),
  d = c("ma", "mp"),
  m = c("full", "envi", "time")) %>%
  mutate(model = paste(p, d, m, "thin", n, "samples_250_chains_4.rds", sep = "_")) %>%
  pull(model)
  
# JSDM mcmc convergence ---------------------------------------------------

# The Gelman/Rubin's diagnostic should be 1, ideally less than 1.05.

mybetas  <- map_dfr(mynames, psrf_ess, "beta")
mygammas <- map_dfr(mynames, psrf_ess, "gamma")
myomegas <- map_dfr(mynames, psrf_ess, "omega")
myrhos   <- map_dfr(mynames, psrf_ess, "rho")

# save convergence statistics
myrhos %>%
  group_by(model, thin) %>%
  summarize(p.est.m=mean(p.est, na.rm=T),
            p.est.sd=sd(p.est, na.rm=T)) %>%
    write_tsv(here::here("tables", "JSDM_convergence", "GR_rho.tsv"))

myrhos %>%
  group_by(model, thin) %>%
  summarize(ess.m=mean(ess, na.rm=T),
            ess.sd=sd(ess, na.rm=T)) %>%
  ungroup() %>%
  write_tsv(here::here("tables", "JSDM_convergence", "ESS_rho.tsv"))

mybetas %>%
  group_by(model, thin) %>%
  summarize(p.est.m=mean(p.est, na.rm=T),
            p.est.sd=sd(p.est, na.rm=T)) %>%
  write_tsv(here::here("tables", "JSDM_convergence", "GR_betas.tsv"))

# JSDM R2 (variance explained) --------------------------------------------

myr2final <- map_dfr(mynames, R2fit)

myr2final %>%
  ungroup() %>%
  group_by(model, thin, phase, distribution) %>%
  summarize(mean_R2=mean(R2, na.rm=T),
         sd_R2=sd(R2, na.rm=T),
         mean_RMSE=mean(RMSE, na.rm=T),
         sd_RMSE=sd(RMSE, na.rm=T),
         mean_AUC=mean(AUC, na.rm=T),
         sd_AUC=sd(AUC, na.rm=T)) %>%
  write_tsv(here::here("tables", "JSDM_fits", "R2_summary_results.tsv"))

# JSDM WAIC ---------------------------------------------------------------
# WAIC stand for widely applicable information criterion

mywaicfinal <- map_dfr(mynames, WAICfit)

mywaicfinal %>%
  group_by(phase, distribution, model, thin) %>%
  summarize(min=min(WAIC_full)) %>%
  write_tsv(here::here("tables", "JSDM_fits", "min_waic.tsv"))


# 5 fold cross validation -------------------------------------------------

# 5f cv takes a long time so we only fit for a subset of the models
finalnames <- c("qeq_ma_full_thin_1000_samples_250_chains_4.rds", 
                "qeq_mp_full_thin_1000_samples_250_chains_4.rds",
                "sort_ma_full_thin_1000_samples_250_chains_4.rds",
                "sort_mp_full_thin_1000_samples_250_chains_4.rds")

mycv5final <- map_dfr(finalnames, cv5fit)

mycv5final %>% 
  filter(R2 > 0) %>%
  group_by(phase, distribution) %>%
  mutate(mean_R2=mean(R2, na.rm=T),
         sd_R2=sd(R2, na.rm=T),
         mean_RMSE=mean(RMSE, na.rm=T),
         sd_RMSE=sd(RMSE, na.rm=T),
         mean_AUC=mean(AUC, na.rm=T),
         sd_AUC=sd(AUC, na.rm=T)) %>%
  write_tsv(here::here("tables", "JSDM_fits", "cv5f_results.tsv"))

# Table S9 ----------------------------------------------------------------

tcv5f <- mycv5final %>%
  ungroup() %>%
  filter(R2 > 0) %>%
  group_by(phase, distribution) %>%
  summarize(mean_predictive_R2=mean(R2, na.rm=T),
            sd_predictive_R2=sd(R2, na.rm=T)) %>%
  mutate(model="full")

tr2 <-myr2final %>%
  filter(thin==1000) %>%
  group_by(phase, distribution, model) %>%
  summarize(mean_explanatory_R2=mean(R2, na.rm=T),
            sd_explanatory_R2=sd(R2, na.rm=T))

# make and format final table
left_join(tr2, tcv5f) %>%
  dplyr::rename(Included_effects=model,
                Model=distribution,
                Phase=phase) %>%
  dplyr::mutate(Phase=if_else(Phase=="sorting", "Sorting", "Equil"),
                Model=if_else(Model=="normal", "COP", "PA"),
                Included_effects=case_when(Included_effects == "full" ~ "Full",
                                           Included_effects == "no_fixed_effects" ~ "Random only",
                                           Included_effects == "no_random_effects" ~ "Fixed only")) %>%
  mutate(Phase=factor(Phase, levels=c("Sorting", "Equil")),
         Model=factor(Model, levels=c("COP", "PA")),
         Included_effects=factor(Included_effects, levels=c("Full", "Fixed only", "Random only"))) %>%
  arrange(Model, Phase, Included_effects) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S9.tex"))
