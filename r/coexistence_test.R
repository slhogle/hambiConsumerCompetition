# LTV
library(here)
library(tidyverse)

dat_ltv <- predator2 %>% filter(treatment=="HNPanc") %>%
  #mutate(ciliate_per_ml=ifelse(cil_interp==1, 0, ciliate_per_ml)) %>%
  #mutate(worm_per_ml=ifelse(wrm_interp==1, 0, worm_per_ml)) %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml)

dat_htv <- predator %>% filter(treatment=="HNPevo") %>%
  mutate(ciliate_per_ml=ifelse(cil_interp==1, 0, ciliate_per_ml)) %>%
  mutate(worm_per_ml=ifelse(wrm_interp==1, 0, worm_per_ml)) %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml) %>%
  mutate(ciliate=case_when(day==5 ~ 10^3,
                           TRUE ~ ciliate)) %>%
  mutate(ciliate=log(ciliate+1), worm=log(worm+1))

dat_ltv <- predator2 %>% filter(treatment=="HNPanc") %>%
  mutate(ciliate_per_ml=ifelse(cil_interp==1, 0, ciliate_per_ml)) %>%
  #arrange(replicate, day) %>%
  #group_by(replicate) %>%
  #mutate(ciliate_per_ml=ifelse(lag(ciliate_per_ml)==0 & lead(ciliate_per_ml)==0, 0, ciliate_per_ml)) %>%
  #replace_na(list(ciliate_per_ml=0)) %>%
  #mutate(worm_per_ml=ifelse(wrm_interp==1, 0, worm_per_ml)) %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml)

### LTV
#A 134660 j=-1.8 k=-1.8
#B 135981
#C 108634
#D 118441

rep="A"
ceiling=134660
LTV_df <- tibble(day=seq(0,4), worm=0, replicate=rep) %>%
  mutate(ciliate=ceiling/(1+exp(-1.8*(day-1.8)))) %>%
  mutate(worm=exp(day-2.5)) %>%
  bind_rows(., dat_ltv) %>%
  mutate(ciliate=log10(ciliate+1), worm=log10(worm+1))# %>%
#filter(replicate==rep)

ggplot(LTV_df) +
  geom_point(aes(x=day, y=ciliate)) + 
  geom_point(aes(x=day, y=worm, color="red")) 



dat_ltv_i <- LTV_df #%>% filter(replicate==rep) 

time_ltv_i <- dat_ltv_i$day
species_ltv_i <- select(dat_ltv_i, ciliate, worm) %>% as.data.frame()

gw = gause_wrapper(time=time_ltv_i, species=species_ltv_i)



gw$parameter_intervals %>%
  select(value=mu) %>%
  rownames_to_column(var="parameter") %>%
  pivot_wider(names_from="parameter", values_from="value") %>%
  mutate(onemp = 1-sqrt((a12*a21)/(a11*a22)),
         kckn = (r1/r2)*sqrt((a22*a21)/(a11*a12))) %>%
  select(onemp, kckn) %>%
  ggplot(aes(x=onemp, y=kckn)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey", size=1) + 
  geom_function(fun = f1) +
  geom_function(fun = f2) +
  geom_point() + 
  scale_y_log10(limits=c(0.1,10)) +
  scale_x_continuous(breaks=seq(0, 1, 0.25), limits=c(0,1)) +
  coord_fixed(xlim=c(0,0.5), ylim=c(0.5,2), expand=FALSE) +
  labs(x="niche difference", 
       y="competitive difference") +
  my_theme

# number of species
N<-ncol(gw$rawdata)-1
# parameters
pars_full<-c(gw$parameter_intervals$mu)
# data.frame for optimization
fittigdata<-data.frame(y=unlist(gw$rawdata[,-1]),
                       time=gw$rawdata$time,
                       N=N)

yest <- ode_prediction(pars_full, time=fittigdata$time, N=fittigdata$N)
plot(fittigdata$y, yest, xlab="observation", ylab="prediction")
abline(a=0, b=1, lty=2)

# HTV

dat_htv <- predator2 %>% filter(treatment=="HNPevo") %>%
  #mutate(ciliate_per_ml=ifelse(cil_interp==1, 0, ciliate_per_ml)) %>%
  mutate(worm_per_ml=ifelse(wrm_interp==1, 0, worm_per_ml)) %>%
  select(day, replicate, ciliate=ciliate_per_ml, worm=worm_per_ml) %>%
  mutate(ciliate=case_when(day==5 ~ 10^3,
                           TRUE ~ ciliate))

### LTV
#A 134660
#B 135981
#C 108634
#D 118441

rep="A"
ceiling=1346
HTV_df <- tibble(day=seq(0,4), worm=0, replicate=rep) %>%
  mutate(ciliate=ceiling/(1+exp(-1.5*(day-2)))) %>%
  #mutate(worm=exp(day-2.5)) %>%
  bind_rows(., dat_htv) %>%
  mutate(ciliate=log10(ciliate+1), worm=log10(worm+1)) #%>%
#filter(replicate==rep)

ggplot(HTV_df) +
  geom_point(aes(x=day, y=ciliate)) + 
  geom_point(aes(x=day, y=worm, color="red")) 

dat_htv_i <- HTV_df #%>% filter(replicate==rep)
time_htv_i <- dat_htv_i$day
species_htv_i <- select(dat_htv_i, ciliate, worm) %>% as.data.frame()

gw = gause_wrapper(time=time_htv_i, species=species_htv_i)

gw$parameter_intervals %>%
  select(value=mu) %>%
  rownames_to_column(var="parameter") %>%
  pivot_wider(names_from="parameter", values_from="value") %>%
  mutate(onemp = 1-sqrt((a12*a21)/(a11*a22)),
         kckn = sqrt((a22*a21)/(a11*a12))) %>%
  select(onemp, kckn) %>%
  ggplot(aes(x=onemp, y=kckn)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey", size=1) + 
  geom_vline(xintercept=0, linetype="dashed", color="grey", size=1) + 
  geom_function(fun = f1) +
  geom_function(fun = f2) +
  geom_point() +
  annotate(geom="point", x=0.3723788, y=0.5040107) +
  #scale_y_log10(limits=c(0.5,3)) +
  scale_x_continuous(breaks=seq(-1, 1, 0.25), limits=c(-1,1)) +
  scale_y_continuous(limits=c(0.5,3)) +
  #coord_fixed(xlim=c(0,0.5), ylim=c(0.5,2), expand=FALSE) +
  labs(x="niche difference", 
       y="competitive difference") +
  my_theme

# PLOT

tplot = tibble(p=c(0.2301603, 0.286703, 0.1632033, 0.1541913, 0.005403006, 0.6686968,
                   0.315572, 0.2806075, 0.2867988, 0.2584979), 
               k=c(0.6193553, 0.6635671, 0.7298156, 0.6740984, 1.03432, 3.600045,
                   0.8771774, 0.894811, 0.8874714, 0.8962965),
               type=c(rep("none", 6), rep("w_interp", 4)))

f1 = function(x){1-x}
f2 = function(x){(1-x)^-1}

ggplot() + 
  geom_hline(yintercept=1, linetype="dashed", color="grey", size=1) + 
  geom_function(fun = f1) +
  geom_function(fun = f2) +
  geom_point(data=tplot, aes(x=p, y=k, color=type)) + 
  scale_y_log10(limits=c(0.1,10)) +
  scale_x_continuous(breaks=seq(0, 1, 0.25), limits=c(0,1)) +
  coord_fixed(xlim=c(0,0.5), ylim=c(0.5,2), expand=FALSE) +
  labs(x="niche difference", 
       y="competitive difference") +
  my_theme