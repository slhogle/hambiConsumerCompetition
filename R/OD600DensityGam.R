library(tidyverse)
library(here)
library(mgcv)
library(itsadug)
library(performance)
library(gratia)
library(rcartocolor)


# LOAD DATA

predator <- read_tsv(here("data", "predator_density.tsv"), col_types="ccccdddd") %>%
  mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
         microcosmID=factor(microcosmID))


# FIRST LOOK

predator %>%
  filter(!is.na(OD)) %>%
  ggplot(aes(x=day, y=OD, color=treatment, group=microcosmID)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=0.5) + # starting concentration
  scale_fill_carto_d(palette = "Vivid") + 
  scale_color_carto_d(palette = "Vivid") + 
  labs(y="Ciliates per mL", x="Day") + 
  #scale_y_log10() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# GAMS

## FORMULAS
# GAM formula specification

# f1 - most complex model
# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations
# c) addition of a microcosm-specific smooth with different wiggliness for each microcosms
# d) random effect smooth for each of the 24 microcosms


f1 <- formula(OD ~ treatment + 
                s(day, k=12, m=2, bs="tp") +
                s(day, k=12, by=treatment, bs="tp", m=1) + 
                s(microcosmID, bs="re", k=12) +
                1)


# f2 - no random effect smooth
# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations
# c) addition of a microcosm-specific smooth with different wiggliness for each sample


f2 <- formula(OD ~ treatment + 
                s(day, k=12, m=2, bs="tp") +
                s(day, k=12, by=treatment, bs="tp", m=1) +
                1)


# f3 - no sample specific time smooth
# a) fixed effect of evolution/predation interaction 
# b) one smoother for time for all observations


f3 <- formula(OD ~ treatment + 
                s(day, k=12, m=2, bs="tp") +
                1)



## FITS
# Fit the models. Comparing multiple probability distributions


input <- dplyr::select(predator, -worm_per_ml, -ciliate_per_ml) %>%
  drop_na()


### GAUSSIAN ID-LINK

#control = gam.control(trace = TRUE)
m01 <- bam(f1, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)

m02 <- bam(f2, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)

m03 <- bam(f3, data=input, family=gaussian(link="identity"), method="fREML", 
           discrete=TRUE)


### GAUSSIAN LOG-LINK

#control = gam.control(trace = TRUE)
m04 <- bam(f1, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)

m05 <- bam(f2, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)

m06 <- bam(f3, data=input, family=gaussian(link="log"), method="fREML", 
           discrete=TRUE)


### GAMMA ID-LINK

#control = gam.control(trace = TRUE)
m07 <- bam(f1, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)

m08 <- bam(f2, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)

m09 <- bam(f3, data=input, family=Gamma(link="identity"), method="fREML", 
           discrete=TRUE)


### GAMMA INVERSE-LINK

#control = gam.control(trace = TRUE)
m10 <- bam(f1, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE)

m11 <- bam(f2, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE)

m12 <- bam(f3, data=input, family=Gamma(link="inverse"), method="fREML", 
           discrete=TRUE)


### GAMMA LOG-LINK

#control = gam.control(trace = TRUE)
m13 <- bam(f1, data=input, family=Gamma(link="log"), method="fREML", 
           discrete=TRUE)

m14 <- bam(f2, data=input, family=Gamma(link="log"), method="fREML",
           discrete=TRUE)

m15 <- bam(f3, data=input, family=Gamma(link="log"), method="fREML",
           discrete=TRUE)


### INVERSE GAUSSIAN INVERSE QUADRATIC-LINK

m16 <- bam(f1, data=input, family=inverse.gaussian(link="1/mu^2"), method="fREML", 
           discrete=TRUE)

m17 <- bam(f2, data=input, family=inverse.gaussian(link="1/mu^2"), method="fREML", 
           discrete=TRUE)


### INVERSE GAUSSIAN LOG-LINK

m19 <- bam(f1, data=input, family=inverse.gaussian(link="log"), method="fREML", 
           discrete=TRUE)

m20 <- bam(f2, data=input, family=inverse.gaussian(link="log"), method="fREML", 
           discrete=TRUE)

m21 <- bam(f3, data=input, family=inverse.gaussian(link="log"), method="fREML", 
           discrete=TRUE)


### INVERSE GAUSSIAN INVERSE-LINK

m22 <- bam(f1, data=input, family=inverse.gaussian(link="inverse"), method="fREML", 
           discrete=TRUE)

m23 <- bam(f2, data=input, family=inverse.gaussian(link="inverse"), method="fREML", 
           discrete=TRUE)

m24 <- bam(f3, data=input, family=inverse.gaussian(link="inverse"), method="fREML", 
           discrete=TRUE)


### INVERSE GAUSSIAN ID-LINK

m25 <- bam(f1, data=input, family=inverse.gaussian(link="identity"), method="fREML", 
           discrete=TRUE)

m26 <- bam(f2, data=input, family=inverse.gaussian(link="identity"), method="fREML", 
           discrete=TRUE)

m27 <- bam(f3, data=input, family=inverse.gaussian(link="identity"), method="fREML", 
           discrete=TRUE)



## COMPARE MODELS

compare_performance(m01, m04, m07, 
                    m10, m13, m16,
                    m19, m22, m25,
                    rank = TRUE)



compare_performance(m02, m05, m08, 
                    m11, m14, m17,
                    m20, m23, m26,
                    rank = TRUE)




compare_performance(m03, m06, m09, 
                    m12, m15,
                    m21, m24, m27,
                    rank = TRUE)

# Looks like guassian with log-link is best

# Which formula is best?
  
  
compare_performance(m04, m05, rank = TRUE)



test_performance(m04, m05)



# So including random effects term is not really any better BF = 1

# The microcosm-level random effects don't appear to be necessary.

# Proceed with formula f2 using the gaussian distribution with log-link.

# After doing a bunch of gam.checks with different knot values (k) it appears that k=9 is the best. Or at least it indicates there is no longer any additional nonlinearity or structure in the residuals that can be explained by a further smooth of time.


gam.check(m05)



# VISUALIZE MODEL SMOOTHS

draw(m05)


# VISUALIZE MODEL DIAGNOSTICS

appraise(m05)
