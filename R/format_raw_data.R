library(here)
library(tidyverse)
library(withr)


# Predator prey densities -------------------------------------------------

df <- read_tsv(here::here("data_raw", "predator_prey_density.tsv"))  %>%
  filter(replicate != 1) %>%
  mutate(replicateABC=LETTERS[replicate-1]) %>%
  mutate(microcosmID=paste0(treatment, "_", replicateABC),
         sample=paste0(treatment, replicate, "-d", day )) %>%
  select(sample, microcosmID, treatment, replicate=replicateABC, day, OD, ciliate_per_ml, worm_per_ml)

with_seed(123784, 
          cil_d <- df %>%
            filter(day > 0) %>%
            filter(!(treatment %in% c("H", "HN"))) %>%
            select(sample, microcosmID, treatment, replicate, day, ciliate_per_ml) %>%
            distinct() %>%
            mutate(ciliate_per_ml_imp = ifelse(ciliate_per_ml == 0, 
                                               rlnorm(100, mean = 5, sd = 0.5), 
                                               ciliate_per_ml)) %>%
            mutate(ciliate_per_ml_imp = round(ciliate_per_ml_imp))
  )

with_seed(123784, 
          wrm_d <- df %>%
            filter(day > 0) %>%
            filter(!(treatment %in% c("H", "HPanc", "HPevo"))) %>%
            select(sample, microcosmID, treatment, replicate, day, worm_per_ml) %>%
            distinct() %>%
            mutate(worm_per_ml_imp = ifelse(worm_per_ml == 0, 
                                            rlnorm(100, mean = 1, sd = 0.5), 
                                            worm_per_ml)) %>%
            mutate(worm_per_ml_imp = round(worm_per_ml_imp))
  )

od_d <- df %>%
  dplyr::select(sample, microcosmID, treatment, replicate, day, OD)

df1 <- left_join(od_d, cil_d) %>%
  left_join(., wrm_d) %>%
  # for setting NAs to 0
  #mutate(across(c(ciliate_count_raw, ciliate_per_ml, 
  #                worm_count_raw, worm_per_ml), ~replace_na(.x, 0))) %>%
  mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
         microcosmID=factor(microcosmID))

write_rds(df1, here::here("data", "formatted_predator_prey_density.rds"))

