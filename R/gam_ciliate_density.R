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
  dplyr::select(-ciliate_per_ml) %>%
  dplyr::rename(ciliate_per_ml=ciliate_per_ml_imp)

ggplot(predator, aes(x=day, y=ciliate_per_ml, color=treatment, group=microcosmID)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=10^4) + # starting concentration
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

f1 <- formula(ciliate_per_ml ~ treatment + 
                s(day, k=13, m=2, bs="tp") +
                s(day, k=13, by=treatment, bs="tp", m=1) + 
                s(microcosmID, bs="re", k=16) +
                1)

# f2 - no random effect smooth

# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations
# c) addition of a microcosm-specific smooth with different wiggliness for each sample

f2 <- formula(ciliate_per_ml ~ treatment + 
                s(day, k=13, m=2, bs="tp") +
                s(day, k=13, by=treatment, bs="tp", m=1) +
                1)

# f3 - no sample specific time smooth

# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations

f3 <- formula(ciliate_per_ml ~ treatment + 
                s(day, k=13, m=2, bs="tp") +
                1)

# GAM fits ----------------------------------------------------------------

# Comparing multiple probability distributions
input <- dplyr::select(predator, -OD, -worm_per_ml, -worm_per_ml_imp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HPanc", "HPevo", "HNPanc", "HNPevo")))

# GAUSSIAN ID-LINK

# control = gam.control(trace = TRUE)
with_seed(13513,
  m01 <- bam(f1, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)
  )

with_seed(13513,
  m02 <- bam(f2, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)
)

with_seed(13513,
  m03 <- bam(f3, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)
)


# GAUSSIAN LOG-LINK

#control = gam.control(trace = TRUE)
with_seed(13513,
          m04 <- bam(f1, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)
          )

with_seed(13513,
          m05 <- bam(f2, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)
          )

with_seed(13513,
          m06 <- bam(f3, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)
          )

# GAMMA ID-LINK

#control = gam.control(trace = TRUE)
with_seed(13513, 
          m07 <- bam(f1, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)
          )

with_seed(13513, 
          m08 <- bam(f2, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)
          )

with_seed(13513, 
          m09 <- bam(f3, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)
          )

# GAMMA INVERSE-LINK

#control = gam.control(trace = TRUE)
with_seed(13513, 
          m10 <- bam(f1, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m11 <- bam(f2, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m12 <- bam(f3, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE))

# GAMMA LOG-LINK

#control = gam.control(trace = TRUE)
with_seed(13513, 
          m13 <- bam(f1, data=input, family=Gamma(link="log"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m14 <- bam(f2, data=input, family=Gamma(link="log"), method="fREML",
           discrete=TRUE))

with_seed(13513, 
          m15 <- bam(f3, data=input, family=Gamma(link="log"), method="fREML",
           discrete=TRUE))

# NEGATIVE BINOMIAL LOG-LINK
with_seed(13513, 
          m16 <- bam(f1, data=input, family=nb(link="log"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m17 <- bam(f2, data=input, family=nb(link="log"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m18 <- bam(f3, data=input, family=nb(link="log"), method="fREML", 
           discrete=TRUE))

### NEGATIVE BINOMIAL ID-LINK
with_seed(13513, 
          m19 <- bam(f1, data=input, family=nb(link="identity"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m20 <- bam(f2, data=input, family=nb(link="identity"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m21 <- bam(f3, data=input, family=nb(link="identity"), method="fREML", 
           discrete=TRUE))

### NEGATIVE BINOMIAL SQRT-LINK

with_seed(13513, 
          m22 <- bam(f1, data=input, family=nb(link="sqrt"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m23 <- bam(f2, data=input, family=nb(link="sqrt"), method="fREML", 
           discrete=TRUE))

with_seed(13513, 
          m24 <- bam(f3, data=input, family=nb(link="sqrt"), method="fREML", 
           discrete=TRUE))


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

# So including random effects term is not really any better BF = 1
# The microcosm-level random effects don't appear to be necessary.
# Proceed with formula f2 using the negative binomial distribution with log-link.

# Diagnostics -------------------------------------------------------------

# After doing a bunch of gam.checks with different knot values (k) it appears
# that k=7 is the best. Or at least it indicates there is no longer any additional
# nonlinearity or structure in the residuals that can be explained by a 
# further smooth of time.

# classic gam.check from mgcv
gam.check(m17)

# visualize smooths using gratia
draw(m17)

# visualize diagnostics using gratia
appraise(m17)

# Supplementary table S2 --------------------------------------------------

gamtabs(m17, caption='Summary of m1lnb', type = "latex") %>%
  write_lines(here::here("tables", "table_S2a.tex"))

# emmeans model contrasts
emmeans(m17, ~ treatment, data=input) %>%
  pairs(type = "response", adjust = "bonf") %>%
  xtable::xtable(type = "response") %>%
  print() %>%
  write_lines(here::here("tables", "table_S2b.tex"))

# Save model --------------------------------------------------------------
write_rds(m17, here::here("data", "GAM_ciliate.rds"))

#sessioninfo::package_info(pkgs = NULL)