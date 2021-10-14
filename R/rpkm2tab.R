library(here)
library(tidyverse)

##### for reading in directory full of counts files
##### https://github.com/STAT545-UBC/Discussion/issues/398

tax <- tibble::tribble(
  ~strainID,             ~genus,        ~species,
  "HAMBI-0097",    "Acinetobacter",       "lwoffii",
  "HAMBI-1972",        "Aeromonas",        "caviae",
  "HAMBI-0105",    "Agrobacterium",   "tumefaciens",
  "HAMBI-2160",       "Bordetella",         "avium",
  "HAMBI-0262",    "Brevundimonas",       "bullata",
  "HAMBI-1988",     "Chitinophaga",        "sancti",
  "HAMBI-1287",      "Citrobacter",        "koseri",
  "HAMBI-0403",        "Comamonas",  "testosteroni",
  "HAMBI-2164",      "Cupriavidus",       "necator",
  "HAMBI-1279",           "Hafnia",         "alvei",
  "HAMBI-1299",         "Kluyvera",    "intermedia",
  "HAMBI-3237",       "Microvirga",   "lotononidis",
  "HAMBI-2792",        "Moraxella",         "canis",
  "HAMBI-1292",       "Morganella",      "morganii",
  "HAMBI-1923",         "Myroides",      "odoratus",
  "HAMBI-3031",         "Niabella",  "yanshanensis",
  "HAMBI-2159", "Paraburkholderia",   "caryophylli",
  "HAMBI-2494", "Paraburkholderia",   "kururiensis",
  "HAMBI-2443",       "Paracoccus", "denitrificans",
  "HAMBI-1977",      "Pseudomonas",  "chlororaphis",
  "HAMBI-0006",      "Pseudomonas",        "putida",
  "HAMBI-1896", "Sphingobacterium",  "spiritivorum",
  "HAMBI-1842",      "Sphingobium",    "yanoikuyae",
  "HAMBI-2659", "Stenotrophomonas",   "maltophilia"
)

files <- set_names(list.files(here("bbmap_rpkm"), full.names = TRUE), 
                   str_extract(list.files( here("bbmap_rpkm"), full.names = TRUE),
                               regex("(?<=[/])([^/]+)(?=\\.[^.]+)")))

counts <- map_df(files, read_tsv, 
                 comment = "#", 
                 col_names= c("strainID", "Length", "Bases","Coverage", "count", "RPKM","Frags", "FPKM"),
                 .id = "sample") %>%
  left_join(., tax) %>%
  select(sample, strainID, genus, species, count)

#### accounting for multiple 16S copies in genomes.

## HAMBI-1842 has 4 copies for sure
## HAMBI-2659 has 2 copies for sure

counts1 <- counts #%>%
  # mutate(
  #   count = case_when(
  #     .$strainID == "HAMBI-1842" ~ trunc(.$count/4), # ensures we get integers. important for some alpha div measures
  #     .$strainID == "HAMBI-2659" ~ trunc(.$count/2),
  #     TRUE ~ as.numeric(.$count)
  #   ))

counts1_wide <- counts1 %>% 
  group_by_at(vars(-count)) %>%
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
  spread(key=sample, value=count) %>% # spread
  select(-row_id) %>% # drop the index 
  drop_na()# %>%
  column_to_rownames(var = "strainID")

write_tsv(counts1, here("summary_sheets", "counts.tsv"))
write_tsv(counts1_wide, here("summary_sheets", "counts_wide.tsv"))
