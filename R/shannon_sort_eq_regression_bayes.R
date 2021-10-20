source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_community_composition.R"))

library(patchwork)
library(performance)
library(modelbased)
library(insight)
library(see)
library(ggeffects)
library(emmeans)
library(withr)
library(rstanarm)
library(bayestestR)

# Load data ---------------------------------------------------------------

# load divnet estimate
div.sum <- read_rds(here::here("data", "shannon_summary.rds"))

# Format data -------------------------------------------------------------

sh1 <- div.sum %>% 
  filter(day < 17) %>% 
  filter(inference=="plugin")

sh2 <- div.sum %>% 
  filter(day > 16) %>% 
  filter(inference=="plugin")


# Sorting models ----------------------------------------------------------

withr::with_seed(12376,
                 m2stan <- stan_glm(estimate~poly(day, df=2)*treatment,
                                    data = sh1,
                                    family=inverse.gaussian(link = "1/mu^2"),
                                    iter = 4000,
                                    cores = 4))


# Sorting model results ---------------------------------------------------


# slope rope
slopeROPE <- pairs(emtrends(m2stan, "treatment", var = "day", 
         max.degree = 1, transform="response")) %>%
  as.data.frame() %>%
  summarize(rope=0.1*sd(estimate)) %>% pull()


# model posterior
# Table S6a
describe_posterior(m2stan, test = c("pd", "rope"), ci = 0.95, rope_range = c(-0.1,0.1), rope_ci=1) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S6a.tex"))


# estimate slope
# Table S6b
estimate_slopes(m2stan, 
                trend="day", 
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
  select(treatment, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S6b.tex"))


# contrasts of slopes
# Table 6c

# Can use estimate_slopes from `modelbased` package
#estimate_slopes(m2stan, trend="day", levels="treatment", test="pd", standardize_robust=T) 

# This is uses emmeans under the hood
describe_posterior(pairs(emtrends(m2stan, 
                                  "treatment", 
                                   var = "day",
                                  max.degree = 1,
                                  data=sh1,
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
  arrange(desc(abs(Median))) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S6c.tex"))


# Sorting plots -----------------------------------------------------------

# model predictions
grid <- visualisation_matrix(sh1, target = c("day", "treatment"), length=100)
estim <- estimate_relation(m2stan, data=grid, modulate = "day", ci = 0.95, length=1000) #modulate = "day", ci = 0.95, length=100

p1 <- ggplot(estim, aes(x = day, y = Predicted)) +
  geom_ribbon(aes(group=treatment, fill = treatment, ymin = CI_low, ymax = CI_high), alpha = 0.1) +
  geom_line(aes(group=treatment, color = treatment)) +
  geom_jitter(data = sh1, aes(y = estimate, color=treatment), width = 0.2, height = 0.1) +
  ylab("Shannon Diversity") +
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  theme_modern()

p1

# model contrasts
constrasts <- estimate_contrasts(m2stan, 
                                 modulate = "day", 
                                 length = 100, 
                                 transform="response")

constrasts$Contrast <- paste0(constrasts$Level1, "-", constrasts$Level2)

# visualize the changes in the differences
p2 <- constrasts %>%
  filter(pd >= 0.97) %>%
  ggplot(aes(x = day, y = Difference)) +
  geom_ribbon(aes(fill = Contrast, ymin = CI_low, ymax = CI_high), alpha = 0.1) +
  geom_line(aes(colour = Contrast), size=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_carto_d(palette="Vivid") + 
  scale_fill_carto_d(palette="Vivid") + 
  ylab("Difference") +
  theme_modern()

p2


# Equilibrium models ------------------------------------------------------


withr::with_seed(12378, 
            m2Bstan <- stan_glm(estimate~treatment,
                  data = sh2,
                  family=inverse.gaussian(link = "1/mu^2"),
                  iter = 4000,
                  cores = 4,
                  seed = 12345))



# Equilibrium model results -----------------------------------------------

describe_posterior(m2Bstan, test = c("pd", "rope"), ci = 0.95, rope_range = c(-0.1,0.1), rope_ci=1) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S7a.tex"))

# estimate means on response scale
estimate_means(m2Bstan)


constrasts$Contrast <- paste0(constrasts$Level1, "-", constrasts$Level2)

estimate_contrasts(m2Bstan, transform="response") %>%
  as_tibble() %>%
  mutate(Difference=round(Difference, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=slopeROPE,
         ROPE_Percentage=ROPE_Percentage*100,
         Contrast=paste0(Level1, "-", Level2)) %>%
  select(Contrast, Difference, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S7b.tex"))



# Equilibrium plots -------------------------------------------------------


estim <- estimate_relation(m2Bstan, modulate = "treatment", ci = 0.95)

p3 <- ggplot(estim, aes(x = treatment, y = Predicted)) +
  geom_pointrange(aes(x=treatment, y=Predicted, ymin=CI_low, ymax=CI_high,
                      group=treatment)) +
  geom_pointrange(data=sh2, aes(x=treatment, y=estimate, ymin=lower, ymax=upper, 
                      color=treatment, group=treatment), 
                  fatten=1, position=position_jitter(width=0.1, height = 0)) +
  ylab("Shannon Diversity") +
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  theme_modern()

p3


# Fig. S3 -----------------------------------------------------------------

pfinal <- p1 + p2 + p3 +
  plot_layout(guides="collect") + 
  plot_annotation(tag_levels = 'A')

ggsave(here("figs", "figS3.svg"), pfinal, width=25.8, height=8, units="cm",
       device="svg")
