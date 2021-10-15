library(here)
library(tidyverse)
library(gauseR)
library(withr)

# Load data ---------------------------------------------------------------
predator <- read_tsv(here("data", "formattedPredatorPreyDensity.tsv"), col_types="ccccdddd") %>%
  mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
         microcosmID=factor(microcosmID))

dat_ltv <- predator %>% filter(treatment == "HNPanc") %>%
  select(day, replicate, ciliate = ciliate_per_ml, worm = worm_per_ml) %>%
  mutate(ciliate = log(ciliate + 1), worm = log(worm + 1)) %>%
  mutate(ciliate = case_when(ciliate < 6 ~ ciliate - 4,
                             TRUE ~ ciliate))

dat_htv <- predator %>% filter(treatment == "HNPevo") %>%
  select(day, replicate, ciliate = ciliate_per_ml, worm = worm_per_ml) %>%
  mutate(ciliate = case_when(day == 5 ~ 10 ^ 3,
                             TRUE ~ ciliate)) %>%
  mutate(ciliate = log(ciliate + 1), worm = log(worm + 1))


# Functions ---------------------------------------------------------------

subsetGW <- function(data, i){
  data_i    <- data %>% filter(replicate == i)
  time_i    <- data_i$day
  species_i <- select(data_i, ciliate, worm) %>% as.data.frame()
  
  mygause <- gause_wrapper(time = time_i,
                           species = species_i,
                           doplot = FALSE)
  
  return(mygause)
}

fitnessNicheDiff <- function(mygause){
  mygause$parameter_intervals %>%
    select(value = mu) %>%
    rownames_to_column(var = "parameter") %>%
    pivot_wider(names_from = "parameter", values_from = "value") %>%
    mutate(rho = sqrt((a12 * a21) / (a11 * a22)),
           kckn = (r1 / r2) * sqrt((a22 * a21) / (a11 * a12))) %>%
    select(rho, kckn)
}

GWloop <- function(i, df1=dat_ltv, df2=dat_htv){
  a <- subsetGW(df1, i)
  b <- fitnessNicheDiff(a) %>% mutate(type = "LTV")
  c <- subsetGW(df2, i)
  d <- fitnessNicheDiff(c) %>% mutate(type = "HTV")
  
  return(bind_rows(b, d))
}

# Run ---------------------------------------------------------------------

reps <- c("A", "B", "C", "D")

with_seed(453785,
          df <- map_df(reps, GWloop, df1=dat_ltv, df2=dat_htv))

df_sum <- df %>% group_by(type) %>%
  filter(rho > 0) %>%
  summarize(
    mrho   = mean(rho, na.rm = T),
    mkckn  = mean(kckn, na.rm = T),
    sdrho  = sd(rho, na.rm = T),
    sdkckn = sd(kckn, na.rm = T)
  )

saveRDS(df_sum, here::here("data", "competitiveLV.rds"))