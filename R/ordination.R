source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_community_composition.R"))

library(withr)
library(patchwork)
library(ggforce)
library(vegan)

library(rstanarm)
library(bayestestR)
library(modelbased)
library(cvequality)

# Ordination --------------------------------------------------------------

a <- counts %>%
  filter(day != 0) %>%
  select(sample, strainID, count.correct) %>%
  group_by(sample) %>%
  mutate(count.correct=count.correct/sum(count.correct)) %>%
  ungroup() %>%
  pivot_wider(names_from="strainID", values_from="count.correct") %>%
  column_to_rownames(var="sample") %>%
  data.frame()

# nMDS ordination ---------------------------------------------------------

# Bray-curtis dissimilarity 
withr::with_seed(12378, dist <- vegan::vegdist(a,  method = "bray"))      

# ordination
withr::with_seed(123784, nmds_ord <- vegan::metaMDS(dist, k = 2, trymax = 100, trace = F, autotransform = FALSE))

# Plot --------------------------------------------------------------------

# tibble of all samples
nmds2plot <- data.frame(nmds_ord$points) %>%
  rownames_to_column(var="sample") %>%
  left_join(., distinct(select(counts, sample:wrm))) %>%
  mutate(method="NMDS") %>%
  dplyr::rename(axis_1 = MDS1, axis_2 = MDS2) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))

# get centroids for each treatment/day
centroids <- nmds2plot %>% 
  group_by(treatment, day) %>%
  mutate(NMDS1=mean(axis_1), NMDS2=mean(axis_2),
         NMDS1sd=sd(axis_1), NMDS2sd=sd(axis_2)) %>%
  ungroup() %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))

# for making spider plot
spiders <- nmds2plot %>% 
  select(sample, fromx=axis_1, fromy=axis_2) %>%
  left_join(., centroids) %>%
  dplyr::rename(tox=NMDS1, toy=NMDS2) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))

centroids <- centroids %>%
  select(NMDS1, NMDS2, NMDS1sd, NMDS2sd, day, treatment) %>% distinct()

p2 <- ggplot() + 
  geom_point(data=nmds2plot, aes(x=axis_1, y=axis_2, color=factor(day)), size=0.5) +
  geom_segment(data=spiders, aes(x=fromx, xend = tox, y=fromy, yend = toy,
                                 group=treatment, color=factor(day))) +
  geom_link2(data = centroids, aes(x=NMDS1, y=NMDS2, color=factor(day), group=treatment)) +
  geom_point(data = centroids, aes(x=NMDS1, y=NMDS2, fill=factor(day)), size=1, shape=21) + #size=((NMDS1sd+NMDS2sd)/2)^-1
  facet_grid(~treatment) + 
  scale_color_brewer(type="qual", palette="Paired") +
  scale_fill_brewer(type="qual", palette="Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# save subplot
readr::write_rds(p2, here::here("data", "fig3c.rds"))

# Jump lengths ------------------------------------------------------------

# https://www.sciencedirect.com/science/article/pii/S1369527418300092#bib0070

# Are the jumps through community space larger or smaller for a given treatment?

# calculate distances by pythagorean theorem
jumps <- nmds2plot %>%
  arrange(treatment, replicate, day) %>%
  group_by(microcosmID) %>%
  mutate(axis1_lag = lag(axis_1),
         axis2_lag = lag(axis_2)) %>%
  mutate(d=sqrt((axis_1-axis1_lag)^2 + (axis_2-axis2_lag)^2)) %>%
  ungroup()


# Jump length regression --------------------------------------------------

withr::with_seed(34578, 
                 m2stan <- stan_glm(
                   d ~ treatment,
                   data = jumps,
                   family = Gamma(),
                   iter = 4000,
                   cores = 4
                 ))

# describe posterior
# Table S8a
describe_posterior(m2stan, test = c("pd", "rope"), ci = 0.95,
                   rope_range = c(-0.1, 0.1), rope_ci=1) %>%
  as_tibble() %>%
  mutate(Median=round(Median, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100,
         ROPE_Percentage=ROPE_Percentage*100) %>%
  select(Parameter, Median, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S8a.tex"))

# Estimate contrasts
# Table S8b
estimate_contrasts(m2stan, 
                   levels="treatment", 
                   test=c("p_direction", "rope"), 
                   ci = 0.95,
                   rope_range = c(-0.1, 0.1),
                   rope_ci = 1) %>%
  tibble() %>%
  filter(pd > 0.97) %>%
  mutate(Difference=round(Difference, 2),
         CI=paste0("[", round(CI_low, 2),", ",round(CI_high, 2),"]"),
         pd=pd*100) %>%
  mutate(ROPE_low=0.1,
         ROPE_Percentage=ROPE_Percentage*100,
         treatment=paste(Level1, Level2, sep=" - ")) %>%
  select(treatment, Difference, CI, pd, ROPE_low, ROPE_Percentage) %>%
  xtable::xtable(auto=TRUE) %>%
  print() %>%
  write_lines(here::here("tables", "table_S8b.tex"))

# Default community with no predator has greatest jump lengths.
# CHTV, N, and HNCHTV all have smaller jumps through community space than HPanc

# Are the treatments differentially variable? 
# Does one condition vary more than we expect by chance

jumps1 <- jumps %>% select(d, treatment) %>% drop_na() %>% data.frame()
print(with(jumps1, asymptotic_test(d, treatment)))
print(with(jumps1, mslr_test(nr = 1e4, d, treatment)))

# No evidence that they are differentially variable