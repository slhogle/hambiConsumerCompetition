library(here)
library(withr)
library(tidyverse)

library(see)
library(insight)
library(bayestestR)
library(performance)

library(rstanarm)
library(bayesplot)
library(tidybayes)

library(emmeans)


# Read data ---------------------------------------------------------------

dq1.f <- readr::read_rds( here::here("data", "community_dissimilarities.rds"))

# Beta regression ---------------------------------------------------------

# Beta regression from rstanarm
with_seed(2344523, 
          m <- stan_glm(value ~ index,
              data = dq1.f,
              family = mgcv::betar,
              iter = 4000,
              cores = 4))

# Plot the model results. Nothing in 95% ROPE
p1 <- m %>%
  gather_draws(`(Intercept)`, indexshannon_div) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.181))) +
  stat_halfeye() +
  geom_vline(xintercept = c(-0.181, 0.181), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue"), name="95% ROPE") +
  labs(x="Condition mean", y="Condition") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# Table results -----------------------------------------------------------

# Table 5
describe_posterior(m, ci = 0.95, rope_ci = 0.95,
                   test = c("p_direction", "rope"))

emmeans(m, ~ index) %>%
  describe_posterior(ci = 0.95, rope_ci = 0.95, rope_range = c(-0.18, 0.18),
                     test = c("p_direction", "rope")) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=0.18,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S5.tex"))

