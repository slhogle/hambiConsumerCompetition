source(here::here("r", "utils_generic.R"))
source(here::here("r", "utils_community_composition.R"))

library(mcp)
library(rjags)
library(withr)

# Load data ---------------------------------------------------------------

# load divnet estimate
div.sum <- read_rds(here::here("data", "shannon_summary.rds"))


# Multiple change point regression ----------------------------------------

# https://lindeloev.github.io/mcp/index.html

# The idea here is that we assume each microcosm takes some amount of time to 
# reach equilibrium dynamics. We want to do further analysis (HMSC) on each of
# these phases (pre equilirbium vs equilibrium) separately. We use the 
# stabilization of species diversity as an indicator of equilibrium so we 
# fit a two-segment change point regression model. The first segment models 
# community diversity changing with time (the pre equilibrium phase). The second
# segment models diversity only as an intercept so that diversity can oscillates
# around a fixed mean value but that the mean doesn't change with time.

# population mean for changepoint of 12.2438102 +- 3.1680302	days

dpsegin <- div.sum %>% 
  filter(inference=="plugin") %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo",
                                              "HN", "HNPanc", "HNPevo")))

# Define models -----------------------------------------------------------

# model 0
# segment 1: change in intercept and slope modeled on day

# segment 2: population-level change point with treatment-varying change points
# sampled around the population mean.joined segment at cp_1. Slope = 0. 
# Intercept (Int) is estimated

model00 <- list(
 estimate ~ 1 + day, 
 estimate ~ 1 + (1|treatment) ~ 0  
)

with_seed(12378,
          fit00 <- mcp(model00, data = dpsegin, sample = "both", iter = 10000))

fit00$loo <- loo(fit00)

# Summary of the model fit
# summary(fit00)

# Summary of random effects. Treatment-specific deviation from the 
# population mean change point
# ranef(fit00)

# Plot
mcpplot00 <- plot(fit00, facet_by = "treatment", q_predict=F)
mcpplot00

# Check convergence
plot_pars(fit00, pars = "varying", type = "trace", ncol = 3)

plot_pars(fit00, regex_pars = "^cp_1_treatment", type = "dens_overlay", ncol = 3)

# model 01
# This is an alternative model where we allow for a change in slope in segment 2.

# segment 1: change in intercept and slope modeled on day (same as before)

# segment 2: population-level change point with treatment-varying change points 
# sampled around the population mean. Joined segment at cp_1. Slope modeled on 
# day. Intercept (Int) is estimated

model01 = list(
 estimate ~ 1 + day,
 estimate ~ 1 + (1|treatment) ~ 0 + day
)

with_seed(4567,
          fit01 <- mcp(model01, data = dpsegin))

fit01$loo = loo(fit01)

# summary(fit01)

# ranef(fit01)

mcpplot01 <- plot(fit01, facet_by = "treatment", q_predict=F)
mcpplot01

# model 02

# This is an alternative model where we allow for a population changepoint 
# only. No-treatment specific effects

# segment 1: change in intercept and slope modeled on day

# segment 2: population-level change point only. Joined segment at cp_1. 
# Slope modeled on day. Intercept (Int) is estimated

model02 = list(
 estimate ~ 1 + day,
 estimate ~ 1 ~ 0 + day
)

with_seed(12478, 
          fit02 <- mcp(model02, data = dpsegin))

fit02$loo = loo(fit02)

#summary(fit02)

#mcpplot02 <- plot(fit02, q_predict=F)
#mcpplot02

# Compare models ----------------------------------------------------------

# We can use the cross-validation from the loo package to compare the 
# predictive performance of mcp models. Use loo to compute 
# Widely Applicable Information Criterion (WAIC) or 
# Estimated Log Predictive Density (ELPD) for each model, and then 
# compare them using loo::loo_compare().

loo::loo_compare(fit00$loo, fit01$loo, fit02$loo)

# So the first model (“model00”) is prefered since model01 has a smaller ELPD. 
# We conclude that a non-zero slope in the second "quasi-equilibrium" 
# phase does not improve out-of-sample predictions over a zero slope 
#  (that is constant longitudinal mean diversity).

## HYPOTHESIS TEST

hypothesis(fit00, c("`cp_1_treatment[H]` = 0",
                    "`cp_1_treatment[HN]` = 0",
                    "`cp_1_treatment[HPanc]` = 0",
                    "`cp_1_treatment[HPevo]` = 0",
                    "`cp_1_treatment[HNPanc]` = 0",
                    "`cp_1_treatment[HNPevo]` = 0"
                    ))

# Final plot --------------------------------------------------------------

mcpplotfinal <- plot(fit00, facet_by="treatment", q_predict=TRUE)

thresh <- tibble::tribble(
  ~treatment, ~day,       ~event,
        "HN",  17, "stabilized",
     "HPanc",  20, "stabilized",
     "HPevo",  20, "stabilized",
    "HNPanc",  21, "stabilized",
    "HNPevo",  25, "stabilized",
    "HNPevo",  32,       "lost",
    "HNPanc",  11,       "lost"
  ) %>%
  mutate(event=factor(event))

pfinal <- mcpplotfinal +
  facet_wrap(treatment ~ ., nrow= 2) +
  geom_pointrange(data=dpsegin, aes(x=day, y=estimate, ymin = lower, ymax = upper), fatten=1) + 
  geom_vline(data=thresh, aes(xintercept=day, linetype=event)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(here::here("figs", "figS2.svg"), pfinal, width=17.8, height=7.8, units="cm",
       device="svg")
