source(here::here("r", "utils_generic.R"))
source(here::here("r", "utils_community_composition.R"))

library(withr)
library(phyloseq)
library(breakaway)
library(DivNet)

# Not necessary for main analysis. Just exploratory

# Read data ---------------------------------------------------------------

dv.cov1 <- read_rds(here::here("data", "divnet_out.rds"))
summary <- read_rds(here::here("data", "shannon_summary.rds"))

# Format data -------------------------------------------------------------

dvbc <- data.frame(dv.cov1$`bray-curtis`) %>% 
  rownames_to_column(var="sampleA") %>%
  pivot_longer(-sampleA, names_to="sampleB", values_to="mean") %>%
  mutate(sampleA=str_replace_all(sampleA, "_|-", "."),
         sampleB=str_replace_all(sampleB, "_|-", ".")) %>%
  mutate(day=case_when(str_detect(sampleA, "d0") & str_detect(sampleB, "d0") ~ 0,
                       str_detect(sampleA, "d5") & str_detect(sampleB, "d5") ~ 5,
                       str_detect(sampleA, "d9") & str_detect(sampleB, "d9") ~ 9,
                       str_detect(sampleA, "d13") & str_detect(sampleB, "d13") ~ 13,
                       str_detect(sampleA, "d17") & str_detect(sampleB, "d17") ~ 17,
                       str_detect(sampleA, "d21") & str_detect(sampleB, "d21") ~ 21,
                       str_detect(sampleA, "d29") & str_detect(sampleB, "d29") ~ 29,
                       str_detect(sampleA, "d45") & str_detect(sampleB, "d45") ~ 45,
                       str_detect(sampleA, "d61") & str_detect(sampleB, "d61") ~ 61,
                       TRUE ~ NA_real_)) %>% drop_na() %>% 
  mutate(sampleA=str_replace_all(sampleA, "\\d\\.d\\d+", "")) %>%
  mutate(sampleB=str_replace_all(sampleB, "\\d\\.d\\d+", "")) %>%
  distinct() %>%
  mutate(sampleA=factor(sampleA, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  mutate(sampleB=factor(sampleB, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  arrange(sampleA, sampleB) %>%
  group_by(grp = paste0(pmin(as.character(sampleA), as.character(sampleB)),
                        pmax(as.character(sampleA), as.character(sampleB)), 
                        as.character(day))) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-grp) %>%
  mutate(sampleA=factor(sampleA, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  mutate(sampleB=factor(sampleB, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo"))))

dvbc_var <- data.frame(dv.cov1$`bray-curtis-variance`) %>% 
  rownames_to_column(var="sampleA") %>%
  pivot_longer(-sampleA, names_to="sampleB", values_to="variance") %>%
  mutate(sampleA=str_replace_all(sampleA, "_|-", "."),
         sampleB=str_replace_all(sampleB, "_|-", ".")) %>%
  mutate(day=case_when(str_detect(sampleA, "d0") & str_detect(sampleB, "d0") ~ 0,
                       str_detect(sampleA, "d5") & str_detect(sampleB, "d5") ~ 5,
                       str_detect(sampleA, "d9") & str_detect(sampleB, "d9") ~ 9,
                       str_detect(sampleA, "d13") & str_detect(sampleB, "d13") ~ 13,
                       str_detect(sampleA, "d17") & str_detect(sampleB, "d17") ~ 17,
                       str_detect(sampleA, "d21") & str_detect(sampleB, "d21") ~ 21,
                       str_detect(sampleA, "d29") & str_detect(sampleB, "d29") ~ 29,
                       str_detect(sampleA, "d45") & str_detect(sampleB, "d45") ~ 45,
                       str_detect(sampleA, "d61") & str_detect(sampleB, "d61") ~ 61,
                       TRUE ~ NA_real_)) %>% drop_na() %>% 
  mutate(sampleA=str_replace_all(sampleA, "\\d\\.d\\d+", "")) %>%
  mutate(sampleB=str_replace_all(sampleB, "\\d\\.d\\d+", "")) %>%
  distinct() %>%
  mutate(sampleA=factor(sampleA, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  mutate(sampleB=factor(sampleB, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  arrange(sampleA, sampleB) %>%
  group_by(grp = paste0(pmin(as.character(sampleA), as.character(sampleB)),
                        pmax(as.character(sampleA), as.character(sampleB)), 
                        as.character(day))) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-grp) %>%
  mutate(sampleA=factor(sampleA, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))) %>%
  mutate(sampleB=factor(sampleB, levels=rev(c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo"))))


dvbc.fmt <- dvbc %>% 
  arrange(sampleA, sampleB) %>%
  group_by(grp = paste0(pmin(as.character(sampleA), as.character(sampleB)),
                        pmax(as.character(sampleA), as.character(sampleB)), 
                        as.character(day))) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-grp) #%>% 
  # filter(sampleA != sampleB) 

dvbc_var.fmt <- dvbc_var %>% 
  arrange(sampleA, sampleB) %>%
  group_by(grp = paste0(pmin(as.character(sampleA), as.character(sampleB)),
                        pmax(as.character(sampleA), as.character(sampleB)), 
                        as.character(day))) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-grp) %>% 
  filter(sampleA != sampleB) 

# Plots -------------------------------------------------------------------

# plugin shannon diversity

ggplot(summary) + 
  geom_line(aes(x=day, y=estimate, color=treatment, group=microcosmID)) + 
  geom_pointrange(aes(x=day, y=estimate, ymin=lower, ymax=upper, 
                      color=treatment, group=microcosmID), fatten=0.5) + 
  facet_grid(inference~treatment) +
  scale_fill_viridis_c(option="B", name=NULL) + 
  labs(y="Shannon Diversity", x="Day") +
  theme_bw() +
  theme(legend.position = "none")

# bray curtis matrix
dvbc %>% 
  filter(day != 0) %>%
  ggplot(aes(x=sampleA, y=sampleB, fill=mean)) + 
  geom_tile() + 
  facet_grid(~day) +
  scale_fill_viridis_c() + 
  coord_fixed()

# bray curtis matrix w/ time
left_join(dvbc.fmt, dvbc_var.fmt) %>%
  ggplot(aes(x=day, y=mean, color=sampleB)) +
  geom_linerange(aes(ymin=mean-variance*10, ymax=mean+variance*10)) +
  geom_point() + 
  geom_line() +
  facet_grid(sampleA~sampleB) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
