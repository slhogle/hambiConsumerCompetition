source(here::here("r", "utils_consumer_feeding.R"))
source(here::here("r", "utils_generic.R"))

library(patchwork)

# Fig 5 A -----------------------------------------------------------------

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

p5a <- ggplot(df2, aes(x=H, y=freq)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method=lm, se=F, size=0.5, color="red", fullrange=TRUE) +
  geom_point(alpha=0.75, shape=16, aes(color=pred_per_ml)) + 
  labs(x="Prey relative abund. without consumer", 
       y="Prey relative abund. with consumer",
       color="Consumer relative\nabundance") + 
  coord_fixed(xlim=c(0,0.96), ylim=c(0,0.96)) + 
  scale_x_continuous(trans=arcsinsqrt_trans) + scale_y_continuous(trans=arcsinsqrt_trans) +
  scale_color_viridis_c(begin=0, end=0.9) +
  facet_wrap(~treatment, ncol=2) + 
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# Fig 5B ------------------------------------------------------------------

p5b <- traits2 %>%
  mutate(consumertype=factor(consumertype, levels=c("evo", "anc", "wrm"))) %>%
  group_by(consumertype) %>%
  mutate(speciesFct=fct_reorder(Species, prey_clearance, max)) %>%
  ungroup() %>%
  ggplot() +
  geom_hline(yintercept=c(2676.469, 7138.521, 75865.864)) +
  geom_smooth(aes(x=speciesFct, y=prey_clearance, 
                  group=consumertype, color=consumertype), 
              method="loess", se=F, show.legend = FALSE) +
  geom_point(aes(x=speciesFct, y=prey_clearance, 
                 color=consumertype), show.legend = FALSE) +
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                 limits = c(10^3, 10^5)) + 
  scale_color_carto_d(palette="Vivid") +
  labs(x="Prey species", y="Prey clearance [prey/consumer]") +
  annotation_logticks(sides = "l", outside=T) +
  coord_cartesian(clip = "off") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

# Final plot --------------------------------------------------------------

p5 <- p5a + p5b +
  plot_layout(widths=c(1.5, 1), heights=c(1.5, 2), ncol=2, guides = 'collect') + 
  plot_annotation(tag_levels = 'A')

p5

ggsave(here("figs", "fig5.svg"), p5, width=25, height=17.8, units="cm",
       device="svg")
