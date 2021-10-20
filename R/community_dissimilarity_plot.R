source(here::here("R", "utils_generic.R"))
source(here::here("R", "utils_community_composition.R"))

# Format data -------------------------------------------------------------

# Get a list of data frames that contain only species 
# counts in columns. only integers

reps.list <- counts %>%
  select(sample, treatment, day, replicate, strainID, count.correct) %>%
  pivot_wider(names_from="strainID", values_from="count.correct") %>%
  group_by(day, treatment) %>%
  group_split() %>%
  map(., . %>% select(5:28))

reps.list.names <- counts %>%
  select(sample, treatment, day, replicate, strainID, count.correct) %>%
  pivot_wider(names_from="strainID", values_from="count.correct") %>%
  group_by(day, treatment) %>%
  group_keys() %>%
  unite("ID", day:treatment, sep="_") %>% pull(ID)

names(reps.list) <- reps.list.names


# Calculate ---------------------------------------------------------------

# Calculate dissimilarities
Dq1.res <- map_df(reps.list, Dq1) %>%
  pivot_longer(cols=everything(), names_to="sample", values_to = "shannon_div")


# Randomly shuffle species occurrences across treatments and replicates
df1 <- counts %>%
  select(sample, treatment, day, replicate, strainID, count.correct) %>%
  pivot_wider(names_from="strainID", values_from="count.correct") %>%
  select(5:28)

reps.r.list <- list()

for (i in 1:54){
  reps.r.list[[i]] <- sample_n(df1, 4)
}

names(reps.r.list) <- reps.list.names


# Calculate dissimilarities expected if you randomly shuffled the species distributions

Dq1.res.rand <- map_df(reps.r.list, Dq1) %>%
  pivot_longer(cols=everything(), names_to="sample", values_to = "shannon_div_rand")

dq1.f <- left_join(Dq1.res.rand, Dq1.res) %>%
  pivot_longer(cols=-sample, names_to="index") %>%
  mutate(index=factor(index, levels=c("shannon_div_rand", "shannon_div")))


write_rds(dq1.f, here::here("data", "community_dissimilarities.rds"))

# Plot data ---------------------------------------------------------------

p <- ggplot(dq1.f, aes(x=index)) + 
  geom_boxplot(aes(y=value)) +
  labs(x="", y="D'") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(here::here("figs", "figS1b.svg"), p, width=17.8, height=11.8, units="cm",
       device="svg")

print(dq1.f %>%
  group_by(index) %>%
  summarize(m=mean(value),
            sd=sd(value)))