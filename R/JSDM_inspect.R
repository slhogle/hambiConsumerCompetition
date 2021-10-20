source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_JSDM.R"))

library(ape)
library(Hmsc)
library(patchwork)

# SELECTED MODELS
models <-
  c(
    "qeq_ma_full_thin_1000_samples_250_chains_4.rds",
    "qeq_mp_full_thin_1000_samples_250_chains_4.rds",
    "sort_ma_full_thin_1000_samples_250_chains_4.rds",
    "sort_mp_full_thin_1000_samples_250_chains_4.rds"
  )

# Variance partitioning ---------------------------------------------------

# Check how well model explains fraction of variance of each species
vp <- map_dfr(models, vpformat)

# Express in terms of R2 
# Check how well model explains fraction variance for each species
r2 <- map_dfr(models, R2fit)

vpr2 <- left_join(vp, r2) %>%
  mutate(tvar = case_when(distribution=="normal" ~ fvar*R2,
                          distribution=="probit" ~ -fvar*R2)) %>%
  mutate(phase=factor(phase, levels=c("sorting", "quasi-equilibrium")))


# Setup a custom color palette
# sorting phase
snorcols <- c("#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#f66356", "#ee4d5a", "#826dba", "#63589f",
              "#7ccba2", "#46aea0")

snor <- vpr2 %>%
  filter(phase=="sorting" & distribution=="normal") %>%
  mutate(species=fct_reorder2(species, variable, R2)) %>%
  arrange(species) %>%
  distinct(species) %>%
  mutate(species=str_remove(species, "HAMBI-")) %>% pull()

vpr2.s <- vpr2 %>% filter(phase=="sorting" & distribution=="normal") %>% 
  mutate(species=str_remove(species, "HAMBI-")) %>%
  mutate(species=factor(species, level=snor))

# equilibrium phase
qnorcols <- c("#f0f0f0", "#d9d9d9", "#bdbdbd", "#f66356", "#826dba", "#7ccba2")

qnor <- vpr2 %>%
  filter(phase=="quasi-equilibrium" & distribution=="normal") %>%
  mutate(species=fct_reorder2(species, variable, R2)) %>%
  arrange(species) %>%
  distinct(species) %>% 
  mutate(species=str_remove(species, "HAMBI-")) %>% pull()

vpr2.q <- vpr2 %>% filter(phase=="quasi-equilibrium" & distribution=="normal") %>% 
  mutate(species=str_remove(species, "HAMBI-")) %>%
  mutate(species=factor(species, level=qnor))

# Posterior estimate of Betas ---------------------------------------------

b <- map_dfr(models, betaformat, 0.95)

beta.s <- b %>% filter(distribution=="normal") %>%
  filter(!(variable %in% c("(Intercept)", "effort"))) %>%
  filter(phase=="sorting") %>%
  mutate(species=str_remove(species, "HAMBI-")) %>%
  mutate(species=factor(species, level=snor))

beta.q <- b %>% filter(distribution=="normal") %>%
  filter(!(variable %in% c("(Intercept)", "effort"))) %>%
  filter(phase=="quasi-equilibrium") %>%
  mutate(species=str_remove(species, "HAMBI-")) %>%
  mutate(species=factor(species, level=qnor))


# Posterior estimate of Gammas --------------------------------------------
# gamma is how much traits explain species response (Betas)

g <- map_dfr(models, gammaformat, 0.95)

gamma.s <- g %>% filter(distribution=="normal") %>%
  filter(!(variable %in% c("(Intercept)", "effort"))) %>%
  filter(trait != "Intercept") %>%
  filter(phase=="sorting") 

gamma.q <- g %>% filter(distribution=="normal") %>%
  filter(!(variable %in% c("(Intercept)", "effort"))) %>%
  filter(trait != "Intercept") %>%
  filter(phase=="quasi-equilibrium")


# A negative interaction coefficient means that the effect of the combined 
# action of two predictors is less then the sum of the individual effects.
# The concrete interpretation is done best visually by inspecting an
# interaction plot.

tgr <- traitgradient("qeq_ma_full_thin_1000_samples_250_chains_4.rds") 
  
  
ptgr <- tgr %>% mutate(cil = case_when(cil == 0 ~ "None",
                               cil == 1 ~ "CLTV",
                               cil == 2 ~ "CHTV"),
               wrm = case_when(wrm == 0 ~ "None", 
                               wrm == 1 ~ "N")) %>%
  mutate(cil=factor(cil, levels=c("None", "CLTV", "CHTV")),
         wrm=factor(wrm, levels=c("None", "N"))) %>%
  ggplot() + 
  geom_pointrange(aes(x=cil, y=me, ymax=hi, ymin=lo, color=factor(wrm)),
                  position=position_jitter(width=0.2, height=0)) +
  scale_color_carto_d() +
  labs(x="Ciliate consumer", y="Scaled trait posterior", color="Nematode\nconsumer") +
  facet_wrap(~trait, scales="free_y", nrow=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=here::here("figs", "JSDM_qeq_trait_gradient.svg"), plot=ptgr, device="svg",
       units="cm", height=8, width=25.8)

# Combine plots
vplot.s <- myvarplot(vpr2.s, snorcols)
vplot.q <- myvarplot(vpr2.q, qnorcols)

bplot.s <- mybetaplot(beta.s)
bplot.q <- mybetaplot(beta.q)

gplot.s <- mygammaplot(gamma.s)
gplot.q <- mygammaplot(gamma.q)

fig4 <- vplot.s + bplot.s + gplot.s + 
  vplot.q + bplot.q + gplot.q + 
  plot_layout(nrow=2, ncol=3, widths = c(1,1,1), heights=c(1.7, 1), guides = 'collect') +
  plot_annotation(tag_levels="A")

ggsave(filename=here::here("figs", "fig4.svg"), plot=fig4, device="svg",
       units="cm", height=15, width=17.8)

# Variance explained by traits --------------------------------------------

# VP$R2T$Beta tells how much traits explain out of variation among species in 
# their response to each covariate. VP$R2T$Y tells how much traits explain out 
# of variation in species occurrences 

# How much traits explain for Beta in total
r2t <- map_dfr(models, vptraitformat)

r2t %>%
  filter(model=="full" & distribution=="normal") %>%
  filter(variable %nin% c("(Intercept)", "effort")) %>%
  mutate(Experiment_phase=ifelse(phase=="sorting", "Sorting (days 0-13)", "Equilibrium (days 17-61")) %>%
  arrange(desc(Experiment_phase)) %>%
  mutate(Covariate=str_replace_all(variable, c("day" = "time", "cil1" = "CLTV", 
                          "cil2" = "CHTV", wrm1 = "N"))) %>%
  select(Covariate, R2T_beta, R2T_y, Experiment_phase) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S10.tex"))

# Phylogenetic signal -----------------------------------------------------

print("Phylogenetic signal sorting phase")
print(summary(loadhmsccoda("sort_ma_full_thin_1000_samples_250_chains_4.rds")$Rho))

print("Phylogenetic signal equilibrium phase")
print(summary(loadhmsccoda("qeq_ma_full_thin_1000_samples_250_chains_4.rds")$Rho))
