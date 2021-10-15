library(tidyverse)
library(here)
library(rcartocolor)
library(patchwork)
library(mgcv)
library(emmeans)
library(ggeffects)
library(withr)
source(here::here("R", "genericUtils.R"))

# Load data ---------------------------------------------------------------
predator <- read_tsv(here("data", "formattedPredatorPreyDensity.tsv"), col_types="ccccdddd") %>%
  mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
         microcosmID=factor(microcosmID))

cilGAM <- readRDS(here::here("data", "ciliateGAM.rds"))
wrmGAM <- readRDS(here::here("data", "nematodeGAM.rds"))
OD6GAM <- readRDS(here::here("data", "OD600GAM.rds"))
compLV <- readRDS(here::here("data", "competitiveLV.rds"))

dcil <- dplyr::select(predator, -OD, -worm_per_ml, -worm_per_ml_imp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HPanc", "HPevo", "HNPanc", "HNPevo")))

dwrm <- dplyr::select(predator, -OD, -ciliate_per_ml, -ciliate_per_ml_imp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HN", "HNPanc", "HNPevo")))

dbac <- dplyr::select(predator, -worm_per_ml, -ciliate_per_ml, 
                      -ciliate_per_ml_imp, -worm_per_ml_imp) %>%
  drop_na()

# Emmeans -----------------------------------------------------------------
dcil1 <- ggemmeans(cilGAM, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low, 
                pred.high=conf.high, treatment=group)

dwrm1 <- ggemmeans(wrmGAM, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low,
                pred.high=conf.high, treatment=group) 

dbac1 <- ggemmeans(OD6GAM, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low,
                pred.high=conf.high, treatment=group)

# Plot --------------------------------------------------------------------

# jitter params
size=2
shape=16
alpha=0.5
width=0.5

# ciliate
dcil2 <- dcil %>%
  filter(ciliate_per_ml != 0)

pcil <- ggplot(dcil1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  # starting value
  annotate(geom="point", x=0, y=10^4) + 
  geom_jitter(data=dcil2, aes(x=day, y=ciliate_per_ml, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  scale_y_log10() + 
  my_theme    

# nematode
dwrm2 <- dwrm %>%
  filter(worm_per_ml_imp != 0)

pwrm <- ggplot(dwrm1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  annotate(geom="point", x=0, y=10^1) + 
  geom_jitter(data=dwrm2, aes(x=day, y=worm_per_ml, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  scale_y_log10() + 
  my_theme

# bacteria
pbac <- ggplot(dbac1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  geom_jitter(data=dbac, aes(x=day, y=OD, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  my_theme  

# competitive lotka voltera
pcompetitive <- 
  ggplot(compLV) + 
  geom_hline(yintercept=1, linetype="dashed", color="grey", size=1) + 
  geom_function(fun = function(x) (x)) +
  geom_function(fun = function(x) (x)^-1) +
  geom_linerange(aes(y=mkckn, xmin = mrho-sdrho, xmax = mrho+sdrho, color=type)) + 
  geom_linerange(aes(x=mrho, ymin = mkckn-sdkckn, ymax = mkckn+sdkckn, color=type)) + 
  geom_point(aes(x=mrho, y=mkckn, color=type)) +
  scale_y_log10() +
  scale_x_continuous(breaks=seq(0, 1, 0.25)) +
  coord_fixed(xlim=c(0,1), ylim=c(0.1,3), expand=FALSE) +
  labs(x="Niche overlap", 
       y="Competitive ratio") +
  my_theme

# final plot

pfinal <- (pbac / pcil / pwrm + plot_layout(nrow=3, ncol=1, guides = 'collect')) | 
  (pcompetitive / guide_area() + plot_layout(heights=c(2,1))) + plot_annotation(tag_levels = 'A')
pfinal

ggsave(here::here("figs", "Fig2.svg"), plot=pfinal, device="svg", units="cm", width=17.8, height=17.8)
