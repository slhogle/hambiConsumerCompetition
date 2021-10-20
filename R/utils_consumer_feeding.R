library(here)
library(tidyverse)
library(scales)

# Load data ---------------------------------------------------------------

traits   <- readr::read_tsv(here::here("data", "species_traits_11-11-2020.tsv"))
predator <- readr::read_rds(here::here("data", "formatted_predator_prey_density.rds"))
counts   <- readr::read_rds(here::here("data", "normalized_corrected_species_counts.rds")) 

# Prepare trait data ------------------------------------------------------

# https://www.chem.agilent.com/store/biocalculators/calcODBacterial.jsp?_requestid=564535
# For bacterial cell cultures OD600 of 1.0 = 8 x 10^8 cells/ml.

traits1 <- traits %>%
  mutate(ANC_BM=ifelse(ANC_BM == 0, 0.03, ANC_BM)) %>%
  dplyr::select(Species:Worm_BM) %>%
  pivot_longer(ANC:Worm, names_to="a", values_to="consumer_cell_ml") %>%
  pivot_longer(ANC_BM:Worm_BM, names_to="b", values_to="prey_OD600") %>% 
  mutate(prey_cell_ml=prey_OD600*(8*10^8)) %>% # converting OD600 to cell density. Based on e coli
  mutate(b=str_extract(b, "^[:alnum:]+")) %>%
  filter(a==b) %>%
  mutate(consumertype=case_when(a=="ANC" ~ "anc",
                                a=="EVO" ~ "evo",
                                a=="Worm" ~ "wrm")) %>% dplyr::select(-a, -b)

# Convert to rates
traits2 <- traits1 %>% 
  mutate(prey_clearance=prey_cell_ml/consumer_cell_ml,
         consumer_yield=prey_clearance^-1*(10^6),
  ) %>%
  filter(!(Species == "HAMBI-2159" & consumertype=="wrm")) # removing species that didnt grow well with worm

# prepares grazing specificity data
counts_f <- counts %>%
  group_by(sample) %>%
  mutate(count1=sum(count.correct)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq=count.correct/count1) %>% dplyr::select(-count1) %>%
  ungroup()

df <- counts_f %>%
  filter(day > 0) %>%
  filter(freq >= 0.01) %>%
  filter(treatment %in% c("H", "HN", "HPanc", "HPevo")) %>%
  mutate(abund=freq) %>%
  dplyr::select(strainID, abund, treatment, replicate, day) %>% 
  pivot_wider(id_cols=c(strainID, treatment, replicate, day), 
              names_from=treatment, values_from=abund)

df1 <- df %>% 
  dplyr::select(strainID, replicate, day, H, HN, HPanc, HPevo) %>%
  pivot_longer(cols = c(HN, HPanc, HPevo), names_to="treatment", values_to="freq") %>%
  left_join(., filter(predator, treatment %in% c("HPanc", "HPevo", "HN"))) %>% 
  filter(day<=21) %>%
  filter(!(treatment=="HN" & day==21)) %>%
  mutate(treatment=factor(treatment, levels=c("HPanc", "HPevo", "HN"), 
                          labels=c("Canc", "Cevo", "N")))

df2 <- df1 %>%
  group_by(replicate, treatment) %>%
  mutate(pred_per_ml=case_when(!is.na(ciliate_per_ml) ~ ciliate_per_ml/max(ciliate_per_ml),
                               is.na(ciliate_per_ml) ~ worm_per_ml/max(worm_per_ml),
                               TRUE ~ NA_real_)) %>%
  mutate(freq=ifelse(is.na(freq), 0.001, freq),
         H=ifelse(is.na(H), 0.001, H)) %>%
  mutate(prey_removed=H-freq) %>%
  mutate(prey_consumed_pc=case_when(treatment=="N" ~ prey_removed/worm_per_ml,
                                    TRUE ~ (prey_removed/ciliate_per_ml)*1)) 
