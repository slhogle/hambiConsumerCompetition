library(here)
library(tidyverse)
source(here::here("r", "utils_generic.R"))

untar(tarfile="data_raw/16S_amplicon/mapping.tar.gz", exdir="data_raw/16S_amplicon")
untar(tarfile="data_raw/16S_amplicon/mapping/bbmap_rpkm.tar.gz", exdir="data_raw/16S_amplicon/mapping")

mapdir <- here::here("data_raw", "16S_amplicon", "mapping", "bbmap_rpkm")

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
  drop_na()

write_rds(counts, here("data", "species_counts.rds"))
write_rds(counts_wide, here("data", "species_counts_wide.rds"))

# Clean up ----------------------------------------------------------------

unlink("data_raw/16S_amplicon/mapping", recursive = T)
