library(tidyverse)
library(here)

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


# Plotting ----------------------------------------------------------------
library(rcartocolor)
library(scales)

my_colors <- carto_pal(6, "Vivid")
names(my_colors) <- c("HPanc", "HPevo", "HN", "HNPanc", "HNPevo", "H")

my_theme <- theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # e+00 becomes 1
  l <- gsub("e\\+00", "", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}

scientific_10_exp_labels <- trans_format("log10", math_format(10^.x) )
scientific_10_exp_breaks <- trans_format("log10", function(x) 10^x )

# mytheme = function() {
#   theme(
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.line.x = element_line(color = "black"),
#     axis.line.y = element_line(color = "black"),
#     legend.title = element_blank(),
#     legend.background = element_blank(),
#     legend.key = element_blank()
#   )
# }


# Convenience functions ---------------------------------------------------

# opposite of %in% fuction
`%nin%` = Negate(`%in%`)

# creates 95% quantile
quibble95 = function(x, q = c(0.025, 0.5, 0.975)) {
  tibble(x = quantile(x, q), quantile = c("q2.5", "q50", "q97.5"))
}

# scales mean by sd
scale2 <- function(x, na.rm = FALSE) {
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
}

# calculates the coefficient of variation
cv <- function(x, na.rm=TRUE){
  sd(x, na.rm=na.rm)/mean(x, na.rm=na.rm)
}

# for scaling plot axis by arcsin 
arcsinsqrt_trans <- trans_new(
  name = "arcsinsqrt",
  trans = function(x) asin(x^0.5),
  inverse = function(x) (sin(x))^2,
  breaks = breaks_extended(n = 6)
)
