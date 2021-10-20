predator <- readr::read_rds(here::here("data", "formatted_predator_prey_density.rds")) #%>%
  # mutate(replicate=factor(replicate, levels=c("A", "B", "C", "D"))) %>%
  # mutate(treatment=factor(treatment, levels=c("H", "HN", "HPanc", "HPevo", "HNPanc", "HNPevo")),
  #        microcosmID=factor(microcosmID))

counts <- readr::read_rds(here::here("data", "normalized_corrected_species_counts.rds"))

# This function implements equation (16) from Chromosomal barcoding of E. coli 
# populations reveals lineage diversity dynamics at high resolution

Dq1 <- function(.data){
  pseudocount <- 1e-7
  M <- dim(.data)[1]
  species_freq <- apply(.data, 1, function(x) as.numeric(x)/sum(x))
  species_freq[species_freq == 0] <- pseudocount
  species_freq_pool <- as.numeric(colSums(.data)/sum(.data))
  species_freq_pool[species_freq_pool == 0] <- pseudocount
  Meff <- exp(sum(apply(species_freq, 2, 
                        function(x) sum(x*log(x/species_freq_pool), na.rm=TRUE)))/M)
  Diss <- (Meff - 1)/(M-1)
  return(Diss)
}
