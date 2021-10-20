library(tidyverse)
library(here)
library(mgcv)
library(itsadug)
library(performance)
library(gratia)
library(rcartocolor)
library(withr)
library(emmeans)

# Load data ---------------------------------------------------------------
predator <- readr::read_rds(here("data", "formatted_predator_prey_density.rds")) %>%
  dplyr::select(-worm_per_ml) %>%
  dplyr::rename(worm_per_ml=worm_per_ml_imp)

ggplot(predator, aes(x=day, y=worm_per_ml, color=treatment, group=microcosmID)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=10) + # starting concentration
  scale_fill_carto_d(palette = "Vivid") + 
  scale_color_carto_d(palette = "Vivid") + 
  labs(y="Ciliates per mL", x="Day") + 
  scale_y_log10() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# GAM formulas ------------------------------------------------------------

# f1 - most complex model

# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations
# c) addition of a microcosm-specific smooth with different wiggliness for each microcosms
# d) random effect smooth for each of the 24 microcosms

f1 <- formula(worm_per_ml ~ treatment + 
                s(day, k=9, m=2, bs="tp") +
                s(day, k=9, by=treatment, bs="tp", m=1) + 
                s(microcosmID, bs="re", k=15) +
                1)

# f2 - no random effect smooth

# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations
# c) addition of a microcosm-specific smooth with different wiggliness for each sample

f2 <- formula(worm_per_ml ~ treatment + 
                s(day, k=9, m=2, bs="tp") +
                s(day, k=9, by=treatment, bs="tp", m=1) +
                1)

# f3 - no sample specific time smooth

# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations

f3 <- formula(worm_per_ml ~ treatment + 
                s(day, k=9, m=2, bs="tp") +
                1)


# GAM fits ----------------------------------------------------------------

# Fit the models. Comparing multiple probability distributions

input <- dplyr::select(predator, -OD, -ciliate_per_ml, -ciliate_per_ml_imp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HN", "HNPanc", "HNPevo")))

# Gaussian w/ identity link

#control = gam.control(trace = TRUE)
with_seed(13513,
          m01 <- bam(f1,
                     data = input,
                     family = gaussian(link = "identity"),
                     method = "fREML",
                     discrete = TRUE
          ))

with_seed(13513,
          m02 <- bam(f2,
                     data = input,
                     family = gaussian(link = "identity"),
                     method = "fREML",
                     discrete = TRUE
          ))

with_seed(13513,
          m03 <- bam(f3,
                     data = input,
                     family = gaussian(link = "identity"),
                     method = "fREML",
                     discrete = TRUE))

# Gaussian w/ log-link

#control = gam.control(trace = TRUE)
with_seed(13513,
          m04 <- bam(f1,
                     data = input,
                     family = gaussian(link = "log"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513,
          m05 <- bam(f2, 
                     data = input, 
                     family = gaussian(link = "log"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m06 <- bam(f3,
                     data = input,
                     family = gaussian(link = "log"),
                     method = "fREML",
                     discrete = TRUE))

# Gamma w/ identity

#control = gam.control(trace = TRUE)
with_seed(13513,
          m07 <- bam(f1, 
                     data = input,
                     family = Gamma(link = "identity"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m08 <- bam(f2,
                     data = input, 
                     family = Gamma(link = "identity"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m09 <- bam(f3, 
                     data = input, 
                     family = Gamma(link = "identity"), 
                     method = "fREML",
                     discrete = TRUE))

# Gamma w/ inverse link

#control = gam.control(trace = TRUE)
with_seed(13513,
          m10 <- bam(f1, 
                     data = input,
                     family = Gamma(link = "inverse"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m11 <- bam(f2, 
                     data = input,
                     family = Gamma(link = "inverse"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513,
          m12 <- bam(f3, 
                     data = input,
                     family = Gamma(link = "inverse"),
                     method = "fREML",
                     discrete = TRUE))

# Gamma w/ inverse link

#control = gam.control(trace = TRUE)
with_seed(13513,
          m10 <- bam(f1, 
                     data = input,
                     family = Gamma(link = "inverse"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m11 <- bam(f2, 
                     data = input,
                     family = Gamma(link = "inverse"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513,
          m12 <- bam(f3, 
                     data = input,
                     family = Gamma(link = "inverse"),
                     method = "fREML",
                     discrete = TRUE))

# Gamma w/ log link

#control = gam.control(trace = TRUE)
with_seed(13513,
          m13 <- bam(f1, 
                     data = input, 
                     family = Gamma(link = "log"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m14 <- bam(f2, 
                     data = input, 
                     family = Gamma(link = "log"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m15 <- bam(f3, 
                     data = input,
                     family = Gamma(link = "log"),
                     method = "fREML",
                     discrete = TRUE))

# Negative binomial w/ log-link

with_seed(13513, 
          m16 <- bam(f1, 
                     data = input, 
                     family = nb(link = "log"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m17 <- bam(f2, 
                     data = input, 
                     family = nb(link = "log"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m18 <- bam(f3, 
                     data = input, 
                     family = nb(link = "log"), 
                     method = "fREML",
                     discrete = TRUE))


# Negative binomial w/ identity-link

with_seed(13513, 
          m19 <- bam(f1, 
                     data = input, 
                     family = nb(link = "identity"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m20 <- bam(f2, 
                     data = input, 
                     family = nb(link = "identity"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m21 <- bam(f3, 
                     data = input, 
                     family = nb(link = "identity"), 
                     method = "fREML",
                     discrete = TRUE))

# Negative binomial w/ identity-link
with_seed(13513, 
          m22 <- bam(f1, 
                     data = input, 
                     family = nb(link = "sqrt"), 
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m23 <- bam(f2, 
                     data = input, 
                     family = nb(link = "sqrt"),
                     method = "fREML",
                     discrete = TRUE))

with_seed(13513, 
          m24 <- bam(f3, 
                     data = input, 
                     family = nb(link = "sqrt"), 
                     method = "fREML",
                     discrete = TRUE))

# Compare models ----------------------------------------------------------
compare_performance(m01, m04, m07, 
                    m10, m13, m16,
                    m19, m22,
                    rank = TRUE)

compare_performance(m02, m05, m08, 
                    m11, m14, m17,
                    m20, m23,
                    rank = TRUE)

compare_performance(m03, m06, m09, 
                    m12, m15, m18,
                    m21, m24,
                    rank = TRUE)

# Looks like negative binomial is best
# Which formula is best?

compare_performance(m16, m17, rank = TRUE)

test_performance(m16, m17)

# So including random effects term is better.
# Proceed with formula f1 using the negative binomial distribution with log-link.

# Diagnostics -------------------------------------------------------------

# After doing a bunch of gam.checks with different knot values (k) it 
# appears that k=9 is the best. Or at least it indicates there is no longer 
# any additional nonlinearity or structure in the residuals that can be 
# explained by a further smooth of time.

# classic gam.check from mgcv
gam.check(m16)

# visualize smooths using gratia
draw(m16)

# visualize diagnostics using gratia
appraise(m16)

# Supplementary table S4 --------------------------------------------------
gamtabs(m16, caption='Summary of m1lnb', type = "latex") %>%
  write_lines(here::here("tables", "table_S4a.tex"))

# emmeans model contrasts
emmeans(m16, ~ treatment, data=input) %>%
  pairs(type = "response", adjust = "bonf") %>%
  xtable::xtable(type = "response") %>%
  print() %>%
  write_lines(here::here("tables", "table_S4b.tex"))

# Save model --------------------------------------------------------------
write_rds(m16, here::here("data", "GAM_nematode.rds"))

#sessioninfo::package_info(pkgs = NULL)