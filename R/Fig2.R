library(tidyverse)
library(here)
library(rcartocolor)
library(patchwork)
library(mgcv)
library(emmeans)
library(itsadug)
library(xtable)
library(ggeffects)
library(withr)


# LOAD DATA

predator <- read_tsv(here("data", "predator_density.tsv"), col_types="ccccdddd") %>%
  mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
         microcosmID=factor(microcosmID))



predator1 <- read_tsv(here("data", "predator_density1.tsv"), col_types="dccdddddd") %>%
  filter(replicate > 1) %>%
  mutate(replicate1=case_when(replicate==2 ~ "A",
                              replicate==3 ~ "B",
                              replicate==4 ~ "C",
                              replicate==5 ~ "D")) %>%
  mutate(replicate1=factor(replicate1, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo"))) %>%
  select(day, treatment, replicate=replicate1, ciliate_per_ml, worm_per_ml) %>%
  mutate(cil_interp=ifelse(ciliate_per_ml==0, 1, 0),
         wrm_interp=ifelse(worm_per_ml==0, 1, 0)) %>%
  select(day, treatment, replicate, cil_interp, wrm_interp)


predator2 <- left_join(predator, predator1)
predator <- predator2



# FIT MODELS

## CILIATE

dcil <- dplyr::select(predator, -OD, -worm_per_ml, -wrm_interp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HPanc", "HPevo", "HNPanc", "HNPevo")))

fcil <- formula(ciliate_per_ml ~ treatment + 
                  s(day, k=13, m=2, bs="tp") +
                  s(day, k=13, by=treatment, bs="tp", m=1) +
                  1)

with_preserve_seed(mcil <- bam(fcil, data=dcil, family=nb(link="log"), method="fREML", 
                               discrete=TRUE))


## WORM

dwrm <- dplyr::select(predator, -OD, -ciliate_per_ml, -cil_interp) %>%
  drop_na() %>%
  mutate(treatment=factor(treatment, levels=c("HN", "HNPanc", "HNPevo")))

fwrm <- formula(worm_per_ml ~ treatment + 
                  s(day, k=9, m=2, bs="tp") +
                  s(day, k=9, by=treatment, bs="tp", m=1) +
                  1)

with_preserve_seed(mwrm <- bam(fwrm, data=dwrm, family=nb(link="log"), method="fREML", 
                               discrete=TRUE))


## BACTERIA

dbac <- dplyr::select(predator, -worm_per_ml, -ciliate_per_ml, -cil_interp, -wrm_interp) %>%
  drop_na()

fbac <- formula(OD ~ treatment + 
                  s(day, k=12, m=2, bs="tp") +
                  s(day, k=12, by=treatment, bs="tp", m=1) +
                  1)

with_preserve_seed(mbac <- bam(fbac, data=dbac, family=gaussian(link="log"), method="fREML", 
                               discrete=TRUE))




# PLOT

#display_carto_pal(6, "Vivid")

my_colors <- carto_pal(6, "Vivid")

names(my_colors) <- c("HPanc", "HPevo", "HN", "HNPanc", "HNPevo", "H")

my_colors

my_theme <- theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# jitter params
size=2
shape=16
alpha=0.5
width=0.5


When consumers lost

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


## CILIATE

startingcil <- tibble::tribble(
  ~day,  ~pred, ~std.error, ~pred.low, ~pred.high, ~treatment,
  0L, 10^4,         NA,    10^4,     10^4,   "HPanc",
  0L, 10^4,         NA,    10^4,     10^4,   "HPevo",
  0L, 10^4,         NA,    10^4,     10^4,   "HNPanc",
  0L, 10^4,         NA,    10^4,     10^4,   "HNPevo"
)

dcil1 <- ggemmeans(mcil, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low, pred.high=conf.high, treatment=group) #%>%
#bind_rows(startingcil,. )

dcil2 <- dcil %>%
  filter(cil_interp == 0)

pcil <- ggplot(dcil1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  geom_jitter(data=dcil2, aes(x=day, y=ciliate_per_ml, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  geom_vline(data=thresh, aes(xintercept=day, linetype=event)) +
  #ggtitle('Ciliate Biomass') +
  #labs(y="Individuals ml-1", x="Day") +
  #scale_x_continuous(limits=c(0,61)) +
  scale_y_log10() + 
  my_theme    


## NEMATODE

startingwrm <- tibble::tribble(
  ~day,  ~pred, ~std.error, ~pred.low, ~pred.high, ~treatment,
  0L, 10,         NA,    10,     10,   "HN",
  0L, 10,         NA,    10,     10,   "HNPanc",
  0L, 10,         NA,    10,     10,   "HNPevo"
)

dwrm1 <- ggemmeans(mwrm, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low, pred.high=conf.high, treatment=group) #%>%
#bind_rows(startingwrm,. )

dwrm2 <- dwrm %>%
  filter(wrm_interp == 0)

pwrm <- ggplot(dwrm1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  geom_jitter(data=dwrm2, aes(x=day, y=worm_per_ml, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  geom_vline(data=thresh, aes(xintercept=day, linetype=event)) +
  #ggtitle('Ciliate Biomass') +
  #labs(y="Individuals ml-1", x="Day") +
  #scale_x_continuous(limits=c(0,61)) +
  scale_y_log10() + 
  my_theme


## BACTERIA

dbac1 <- ggemmeans(mbac, terms = c("day", "treatment")) %>% 
  as_tibble() %>% 
  dplyr::rename(day=x, pred=predicted, pred.low=conf.low, pred.high=conf.high, treatment=group)

pbac <- ggplot(dbac1) + 
  geom_ribbon(aes(x=day, ymin=pred.low, ymax=pred.high, fill=treatment), alpha=0.2) + 
  geom_line(aes(x=day, y=pred, color=treatment)) + 
  geom_jitter(data=dbac, aes(x=day, y=OD, color=treatment), size=size, shape=shape, alpha=alpha, width=width) + 
  scale_fill_manual(values=my_colors) + 
  scale_color_manual(values=my_colors) + 
  #scale_x_continuous(limits=c(0,61)) +
  geom_vline(data=thresh, aes(xintercept=day, linetype=event)) +
  #ggtitle('Ciliate Biomass') +
  #labs(y="Individuals ml-1", x="Day") +
  #scale_y_log10() + 
  my_theme     


# COEXISTENCE LOTKA VOLTERRA

library(gauseR)


# Following approach from Box 1 in this paper:
#  [Using ecological coexistence theory to understand antibiotic resistance and microbial competition](https://www.nature.com/articles/s41559-020-01385-w)

# The following is verbatim from that paper:
  
#   "Although differences in competitive ability and niche overlap can hypothetically be obtained for any pairwise model of competition, the most convenient formulas, and the most readily used by empirical ecologists, are those derived for phenomenological Lotka–Volterra-type competition models of the general form."

# Thus the abundance of the ciliate would take the form:
  
#   $$ \frac{dN_{cil}}{dt} = N_{cil}(r_{cil} - \alpha_{cil,cil}N_{cil} - \alpha_{cil,wrm}N_{wrm}) $$
#  And the abundance of the nematode would take the form:
#  $$ \frac{dN_{wrm}}{dt} = N_{wrm}(r_{wrm} - \alpha_{wrm,wrm}N_{wrm} - \alpha_{wrm,cil}N_{cil}) $$
#  Here $N_{cil}$ is the abundance of the ciliate, $N_{wrm}$ is the abundance of species the nematode, $t$ is time in days, $r_{cil}$ is the per capita intrinsic rate of increase of the ciliate species, $\alpha_{cil,cil}$ is the linear effect of intraspecific competition (sometimes parameterized as carrying capacity, $K=\frac{1}{\alpha_{cil,cil}}$) and $\alpha_{cil,wrm}$ is the linear effect of interspecific competition of the nematode on the ciliate.

# The niche overlap is then simply the geometric mean ratio of the intraspecici and interspecific coefficients

# $$ \rho = \sqrt\frac{\alpha_{cil,wrm}\alpha_{wrm,cil}}{\alpha_{cil,cil}\alpha_{wrm,wrm}}$$
#   Therefore the niche difference between the two species is just:
  
  $$ 1 - \rho$$
  
#  The competitive ratio of the ciliate to the nematode is given by the ratio of intrinsic growth rates multiplied by the geometric mean ratio of each species competition sensitivity coefficients

#$$ \frac{k_{c}}{k_{n}} = \frac{r_{c}}{r_{n}} \times \sqrt\frac{\alpha_{nn}\alpha_{nc}}{\alpha_{cc}\alpha_{cn}}$$
#  You can then fit empirical data to these Lotka-Volterra equations to derive these coefficients. The `gauseR` package has a convenient wrapper for doing this all at once.

# Log transform abundance for better fits

dat_ltv <- predator %>% filter(treatment=="HNPanc") %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml) %>%
  mutate(ciliate=log(ciliate+1), worm=log(worm+1)) %>%
  mutate(ciliate=case_when(ciliate < 6 ~ ciliate-4,
                           TRUE ~ ciliate))

dat_htv <- predator %>% filter(treatment=="HNPevo") %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml) %>%
  mutate(ciliate=case_when(day==5 ~ 10^3,
                           TRUE ~ ciliate)) %>%
  mutate(ciliate=log(ciliate+1), worm=log(worm+1)) 



# Now we will run a big dumb loop over each individual replicate to calculate these competitive ratios and niche differences
df <- tibble()

with_preserve_seed(
  for (i in c("A", "B", "C", "D")) {
    dat_ltv_i <- dat_ltv %>% filter(replicate == i)
    time_ltv_i <- dat_ltv_i$day
    species_ltv_i <-
      select(dat_ltv_i, ciliate, worm) %>% as.data.frame()
    
    gause_ltv <-
      gause_wrapper(time = time_ltv_i,
                    species = species_ltv_i,
                    doplot = FALSE)
    
    ncltv <- gause_ltv$parameter_intervals %>%
      select(value = mu) %>%
      rownames_to_column(var = "parameter") %>%
      pivot_wider(names_from = "parameter", values_from = "value") %>%
      mutate(onemp = sqrt((a12 * a21) / (a11 * a22)),
             kckn = (r1 / r2) * sqrt((a22 * a21) / (a11 * a12))) %>%
      select(onemp, kckn) %>% mutate(type = "LTV")
    
    dat_htv_i <- dat_htv %>% filter(replicate == i)
    time_htv_i <- dat_htv_i$day
    species_htv_i <-
      select(dat_htv_i, ciliate, worm) %>% as.data.frame()
    
    gause_htv <-
      gause_wrapper(time = time_htv_i,
                    species = species_htv_i,
                    doplot = FALSE)
    
    nchtv <- gause_htv$parameter_intervals %>%
      select(value = mu) %>%
      rownames_to_column(var = "parameter") %>%
      pivot_wider(names_from = "parameter", values_from = "value") %>%
      mutate(onemp = sqrt((a12 * a21) / (a11 * a22)),
             kckn = (r1 / r2) * sqrt((a22 * a21) / (a11 * a12))) %>%
      select(onemp, kckn) %>% mutate(type = "HTV")
    
    df <- bind_rows(df, bind_rows(ncltv, nchtv))
  }
)


# Summarize across replicates

df_sum <- df %>% group_by(type) %>%
  filter(onemp > 0) %>%
  summarize(monemp =  mean((1-onemp), na.rm=T),
            mkckn  =  mean(kckn, na.rm=T),
            sdonemp = sd((1-onemp), na.rm=T),
            sdkckn  = sd(kckn, na.rm=T))


# Final coexistence plot

pcompetitive <- ggplot(df_sum) + 
  geom_hline(yintercept=1, linetype="dashed", color="grey", size=1) + 
  geom_function(fun = function(x) (1-x)) +
  geom_function(fun = function(x) (1-x)^-1) +
  geom_linerange(aes(y=mkckn, xmin = monemp-sdonemp, xmax = monemp+sdonemp, color=type)) + 
  geom_linerange(aes(x=monemp, ymin = mkckn-sdkckn, ymax = mkckn+sdkckn, color=type)) + 
  geom_point(aes(x=monemp, y=mkckn, color=type)) +
  scale_y_log10() +
  scale_x_continuous(breaks=seq(0, 1, 0.25)) +
  coord_fixed(xlim=c(0,1), ylim=c(0.1,3), expand=FALSE) +
  labs(x="niche difference", 
       y="competitive difference") +
  my_theme

pcompetitive


# FINAL PLOT

pfinal <- (pbac / pcil / pwrm + plot_layout(nrow=3, ncol=1, guides = 'collect')) | 
  (pcompetitive / guide_area() + plot_layout(heights=c(2,1))) + plot_annotation(tag_levels = 'A')
pfinal


ggsave(here::here("figs", "densities1.svg"), plot=pfinal, device="svg", units="cm", width=17.8, height=17.8)
