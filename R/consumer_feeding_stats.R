source(here::here("r", "utils_consumer_feeding.R"))
source(here::here("r", "utils_generic.R"))

library(patchwork)

library(see)
library(insight)
library(performance)
library(modelbased)
library(emmeans)

library(rstanarm)
library(bayesplot)
library(tidybayes)
library(bayestestR)
library(withr)

# Grazing specificity from 16S amplicon -----------------------------------

# model
withr::with_seed(123784, 
    mfull <- stan_glm(freq ~ H*treatment,
                  data = df3,
                  iter = 4000,
                  cores = 4))


# describe posterior
# Table S11A
describe_posterior(mfull, ci = 0.95, rope_ci = 0.95,
                   test = c("p_direction", "rope")) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S11a.tex"))

# Estimate ROPE
slopeROPE <- pairs(emtrends(mfull, "treatment", var = "H", 
               max.degree = 1, transform="response")) %>%
  as.data.frame() %>%
  summarize(rope=0.1*sd(estimate)) %>% pull()

# Estimate slopes
# Table S11B
estimate_slopes(mfull, 
                trend="H", 
                levels="treatment", 
                test=c("p_direction", "rope"), 
                ci = 0.95,
                rope_range = c(-slopeROPE, slopeROPE),
                rope_ci = 1,
                standardize_robust=T) %>%
  as_tibble() %>%
  mutate(Median=round(Coefficient, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=slopeROPE,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(treatment, Median, CI, pd, ROPE_low, ROPE_Percentage)# %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S11b.tex"))

# Contrasts of slopes
# Table S11C
describe_posterior(pairs(emtrends(mfull, 
                                  "treatment", 
                                  var = "H",
                                  max.degree = 1,
                                  data=df3,
                                  transform="response")),
                   ci = 0.95,
                   rope_range = c(-slopeROPE, slopeROPE),
                   rope_ci = 1,
                   test = c("p_direction", "rope")) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=slopeROPE,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S11c.tex"))

# Consumer clearance ------------------------------------------------------

# visualize distribution
traits2 %>%
  ggplot(aes(prey_clearance, color=consumertype)) + 
  geom_density() + 
  scale_x_log10()

# coefficient of variation
traits2 %>%
  group_by(consumertype) %>%
  summarize(mn=mean(prey_clearance),
            sd=sd(prey_clearance),
            cv=cv(prey_clearance))

# log transform response variable
traits4 <- traits2 %>% mutate(lprey_clearance=log(prey_clearance))

# model
withr::with_seed(13547, 
                 mtraits <- stan_glm(lprey_clearance ~ consumertype,
                    data = traits4,
                    iter = 2000,
                    cores = 4,
                    seed = 12345))


# describe posterior
# Table S12A
describe_posterior(mtraits, ci = 0.95, rope_ci = 0.95,
                   test = c("p_direction", "rope")) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S12a.tex"))

# Estimate ROPE
traitROPE <- pairs(emmeans(mtraits, "consumertype", transform="response")) %>%
  as.data.frame() %>%
  summarize(rope=0.1*sd(estimate)) %>% pull()

# Model contrasts
# Table S12B
describe_posterior(pairs(emmeans(mtraits, 
                                  "consumertype", 
                                  data=traits4,
                                  transform="response")),
                   ci = 0.95,
                   rope_range = c(-traitROPE, traitROPE),
                   rope_ci = 1,
                   test = c("p_direction", "rope")) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=traitROPE,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S12b.tex"))
