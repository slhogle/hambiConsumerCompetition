source(here::here("r", "utils_generic.R"))
source(here::here("r", "utils_community_composition.R"))

library(patchwork)
library(scales)
library(ggforce)
library(withr)

# Read plot from ordination.R ---------------------------------------------

p2 <- readr::read_rds(here::here("data", "fig3c.rds"))


# Species abdundance distribution -----------------------------------------

# order strains by abundance
strainrecode <- counts %>%
  filter(day > 0) %>%
  group_by(treatment, day, replicate, strainID) %>%
  summarize(spsum=sum(count.correct)) %>%
  ungroup() %>%
  group_by(treatment, day, replicate) %>%
  mutate(nreads=sum(spsum)) %>%
  ungroup() %>%
  mutate(RA=round(spsum/nreads, 2)) %>%
  group_by(strainID) %>%
  summarize(mRA=mean(RA), sdRA=sd(RA)) %>%
  ungroup() %>%
  arrange(mRA) %>%
  distinct(strainID) %>%
  rowid_to_column("strain.num")

strainrecode <- counts %>% 
  group_by(strainID) %>%
  summarize(spsum=sum(count.correct)) %>%
  ungroup() %>%
  mutate(RA=spsum/sum(spsum)) %>%
  arrange(desc(RA)) %>%
  distinct(strainID) %>%
  rowid_to_column("strain.num")


# Setup tibble
sad <- counts %>%
  mutate(count.correct=ifelse(count.correct==0, 1, count.correct)) %>%
  filter(day > 0) %>%
  group_by(treatment, day, replicate, strainID) %>%
  summarize(spsum=sum(count.correct)) %>%
  ungroup() %>%
  group_by(treatment, day, replicate) %>%
  mutate(nreads=sum(spsum)) %>%
  ungroup() %>%
  mutate(RA=round(spsum/nreads, 2)) %>%
  group_by(treatment, strainID) %>%
  summarize(mRA=mean(RA), sdRA=sd(RA)) %>%
  mutate(mRA=ifelse(mRA < 1e-4, 1e-4, mRA)) %>%
  ungroup() %>%
  mutate(ymin=mRA-sdRA, ymax=mRA+sdRA) %>%
  mutate(ymin=ifelse(ymin<=0, 1e-4, ymin)) %>%
  left_join(., strainrecode) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))

# plot
p0 <- ggplot(sad) + 
  geom_linerange(aes(x=strain.num, ymin = ymin, ymax = ymax), color="grey80") + 
  geom_point(aes(x=strain.num, y=mRA), size=0.5) + 
  facet_grid(~ treatment) +
  labs(y="", x="") +
  scale_y_continuous(trans = 'log10',
                     breaks = c(1, 0.5, 0.1, 0.01, 0.001),
                     labels = label_percent(accuracy = 0.1)) + 
  scale_x_continuous(n.breaks=24, limits = c(0,25), expand=c(0,0)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Heat map ----------------------------------------------------------------



# Calculate means over replicates
hmapdf.mn <- counts %>% 
  mutate(count=exp(rlogcount)) %>%
  dplyr::select(-lograrecount, -log10rab, -sflogcount, -rlogcount, -PA) %>%
  arrange(strainID, microcosmID, day) %>%
  group_by(treatment, strainID, day) %>%
  mutate(count.mean=mean(count)) %>% ungroup() %>%
  dplyr::select(day, treatment, cil, wrm, strainID, count.mean) %>%
  distinct() %>%
  group_by(treatment, strainID) %>%
  mutate(count.mean.f=count.mean/max(count.mean)) %>%
  mutate(count.mean.f1=log(count.mean/max(count.mean))/sqrt(2)) %>%
  mutate(xmin=day, xmax=lead(day)) %>%
  mutate(xmax=ifelse(is.na(xmax), 65, xmax)) %>%
  ungroup() %>%
  left_join(., strainrecode) %>%
  mutate(treatment=factor(treatment, levels=c("H", "HPanc", "HPevo", "HN", "HNPanc", "HNPevo")))


p1 <- ggplot(hmapdf.mn) +
  geom_rect(aes(ymin=xmin, ymax=xmax, xmin = strain.num-0.5, 
                xmax = strain.num+0.5, fill = count.mean.f)) +
  geom_hline(yintercept=32) +
  geom_hline(yintercept=11) +
  scale_fill_viridis_c(option = "C", begin = 0.05, end=0.95, direction =-1,
                       trans = scales::pseudo_log_trans(sigma = 1e-3, base = exp(1)),
                       breaks = breaks_log(n = 5, base = 10),
                       labels = label_percent(accuracy = 0.1)) +
  facet_grid(~ treatment) +
  labs(y="", x="") +
  scale_x_continuous(n.breaks=24, limits = c(0,25), expand=c(0,0)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


# Combined plot -----------------------------------------------------------

pfinal <- (p0 / p1 / p2) + plot_layout(guides = 'collect', 
                                       widths = c(1,1,1),
                                       heights= c(2,4,3),
                                       ncol=1)

ggsave(here("figs", "fig3.svg"), pfinal, width=20.8, height=12, units="cm",
       device="svg")
