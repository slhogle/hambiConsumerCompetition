source(here::here("r", "utils_generic.R"))
source(here::here("r", "utils_community_composition.R"))

library(patchwork)
library(performance)
library(modelbased)
library(insight)
library(see)
library(ggeffects)
library(emmeans)
library(withr)

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


# fit glm for the first non-equilibrium phase. Response variable is 
# bounded `0 -> Inf` so `Gamma` and `inverse.gaussian` seem like
# good error distribution families to try

withr::with_seed(123784,
          m1 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=gaussian(link='identity'),
                    data=sh1))
withr::with_seed(123784,          
          m1.1 <- glm(formula=estimate~poly(day, df=2)*treatment,
                      family=gaussian(link='log'),
                      data=sh1))
withr::with_seed(123784,          
          m2 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=inverse.gaussian(link = "1/mu^2"),
                    data=sh1))
withr::with_seed(123784,         
          m3 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=inverse.gaussian(link='identity'),
                    data=sh1))
withr::with_seed(123784,          
          m4 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=inverse.gaussian(link='inverse'),
                    data=sh1))
withr::with_seed(123784,          
          m5 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=inverse.gaussian(link='log'),
                    data=sh1))
withr::with_seed(123784,         
          m6 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=Gamma(link="identity"),
                    data=sh1))
withr::with_seed(123784,          
          m7 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=Gamma(link="inverse"),
                    data=sh1))
withr::with_seed(123784,         
          m8 <- glm(formula=estimate~poly(day, df=2)*treatment,
                    family=Gamma(link="log"),
                    data=sh1))


# Sorting performance -----------------------------------------------------

compare_performance(m1, m1.1)
compare_performance(m1, m1.1, m2, m3, m4, m5, m6, m7, m8, rank=T)

# So Gaussian is bad (as expected).

# Let's check inverse.gaussian only

compare_performance(m2, m3, m4, m5, rank=T)

# seems like inverse.gaussian w/ inverse-quadratic link is best.
# Has lowest RMSE and highest R2

# diagnostics are mostly ok. homogeneity of variance is weird because at 
# day 0 everything starts at the same diversity more or less

check_model(m2)

# Do we need the quadratic term?
m2_noq <- glm(formula=estimate~day*treatment,
    family=inverse.gaussian(link = "1/mu^2"),
    data=sh1)

compare_performance(m2, m2_noq, rank=T)

check_model(m2_noq)



# Sorting emmeans ---------------------------------------------------------


# Effect of the quadratic term

sh1.emt.deg <- emtrends(m2_noq, ~ degree | treatment, var = "day", max.degree = 2)
summary(sh1.emt.deg, infer = c(TRUE, TRUE))

# Quadratic term is important in some cases, 
# but role seems small compared to linear term

# Testing for different linear trends (slopes) between conditions

sh1.emt <- emtrends(m2, specs = pairwise ~ treatment, var = "day", max.degree = 1)
pairs(sh1.emt, adjust = "bonferroni", type = "response")


# Sorting Plots -----------------------------------------------------------

dat <- ggemmeans(m2_noq, terms = c("day [all]", "treatment")) %>%
  rename(treatment=group, day=x)

dat2 <- sh1 %>% select(estimate, lower, upper, treatment, day) %>% distinct()

dat3 <- left_join(dat, dat2) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo",
                                              "HN", "HNPanc", "HNPevo")))

p1 <- ggplot(dat3) + 
  geom_pointrange(aes(x=day, y=estimate, ymin=lower, ymax=upper, 
                      color=treatment, group=treatment), fatten=2, 
                  position=position_jitter(width=0.2, height = 0)) +
  geom_ribbon(aes(x=day, ymin=conf.low, ymax=conf.high,
                  fill=treatment, group=treatment), alpha=0.15) +
  geom_line(aes(x=day, y=predicted, color=treatment, group=treatment)) +
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  theme_modern()

p1

dat <- ggemmeans(m2, terms = c("day [all]", "treatment")) %>%
  rename(treatment=group, day=x)

dat2 <- sh1 %>% select(estimate, lower, upper, treatment, day) %>% distinct()

dat3 <- left_join(dat, dat2)

p1 <- ggplot(dat3) + 
  geom_pointrange(aes(x=day, y=estimate, ymin=lower, ymax=upper, 
                      color=treatment, group=treatment), fatten=2, 
                  position=position_jitter(width=0.2, height = 0)) +
  geom_ribbon(aes(x=day, ymin=conf.low, ymax=conf.high,
                  fill=treatment, group=treatment), alpha=0.15) +
  geom_line(aes(x=day, y=predicted, color=treatment, group=treatment)) +
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  theme_modern()

p1


# Equilibrium models ------------------------------------------------------

withr::with_seed(24671,
                 m2B <- glm(
                   formula = estimate ~ treatment,
                   family = inverse.gaussian(link = "1/mu^2"),
                   data = sh2
                 ))


### EMMEANS

# regrid to put back on orginal scale
m2B.rg <- ref_grid(m2B)
m2B.rg.remm.s <- emmeans(regrid(m2B.rg), "treatment")

sh2.emm <- emmeans(regrid(m2B.rg), specs = pairwise ~ treatment, type = "response", adjust = "mvt", data=sh2)

sh2.emm.pairs <- pairs(sh2.emm, adjust = "fdr", type = "response")

sh2.emm.pairs$emmeans

dat <- ggemmeans(m2B, terms = c("treatment")) %>%
  rename(treatment=x)

dat2 <- sh2 %>% select(estimate, lower, upper, treatment, day) %>% distinct()

dat3 <- left_join(dat, dat2) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo",
                                              "HN", "HNPanc", "HNPevo")))
p2 <- ggplot(dat3) + 
  geom_pointrange(aes(x=treatment, y=predicted, ymin=conf.low, ymax=conf.high,
                      group=treatment)) +
  geom_pointrange(aes(x=treatment, y=estimate, ymin=lower, ymax=upper, 
                      color=treatment, group=treatment), 
                  fatten=1, position=position_jitter(width=0.1, height = 0)) +
  scale_color_manual(values=my_colors) + 
  scale_fill_manual(values=my_colors) + 
  theme_modern()

p2