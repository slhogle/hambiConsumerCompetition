library(here)
library(tidyverse)
source(here::here("R", "genericUtils.R"))

mapdir <- here::here("dataRaw", "16SAmplicon", "mapping", "bbmapRPKM")

files <- set_names(list.files(mapdir, full.names = TRUE),
                   str_extract(
                     list.files(mapdir, full.names = TRUE),
                     regex("(?<=[/])([^/]+)(?=\\.[^.]+)")
                   ))

counts <- map_df(
  files,
  read_tsv,
  comment = "#",
  col_names = c(
    "strainID",
    "Length",
    "Bases",
    "Coverage",
    "count",
    "RPKM",
    "Frags",
    "FPKM"
  ),
  .id = "sample"
) %>%
  left_join(., tax) %>%
  select(sample, strainID, genus, species, count)

counts_wide <- counts %>%
  group_by_at(vars(-count)) %>%
  mutate(row_id = 1:n()) %>% ungroup() %>%  # build group index
  spread(key = sample, value = count) %>% # spread
  select(-row_id) %>% # drop the index
  drop_na()# %>%
column_to_rownames(var = "strainID")

write_tsv(counts, here("data", "speciesCounts.tsv"))
write_tsv(counts_wide, here("data", "speciesCountsWide.tsv"))
