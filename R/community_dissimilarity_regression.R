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

# posterior summary
describe_posterior(m, ci = 0.95, rope_ci = 0.95,
                   test = c("p_direction", "rope"))

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

p1

# Examine emmeans contrasts

m %>%
  emmeans( ~ index) %>%
  gather_emmeans_draws() %>%
  median_qi() %>%
  xtable::xtable() %>%
  print() %>%
  write_lines(here::here("tables", "table_S5.tex"))