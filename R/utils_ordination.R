`%nin%` <- Negate(`%in%`)

logit <- function(x){
  log(x/(1-x))
}

# get dataframe for ordination from tibble
makemydf <- function(tib, d, tr){
  tib %>%
    filter(day==d) %>%
    group_by(sample) %>%
    pivot_wider(id_cols="sample", 
                names_from = "strainID", 
                values_from = tr) %>%
    column_to_rownames(var="sample") %>%
    as.data.frame()
}

# fortify ordination to dataframe
ordfort <- function(ord, metadata){
  fortify(ord) %>%
    rename(sample=Label) %>%
    mutate(sample=as.character(sample)) %>%
    left_join(., metadata)
}

# compute group centroids groupwise, which are the mean coordinate on each axis
groupcentroids <- function(fort_ord){
  fort_ord %>%
    filter(Score=="sites") %>%
    group_by(treatment) %>%
    summarize(PC1=mean(PC1), PC2=mean(PC2))
}

# get the segments for the spider plot
spiderize <- function(fort_ord, centroids){
  fort_ord %>% 
    filter(Score=="sites") %>%
    select(treatment, fromx=PC1, fromy=PC2) %>%
    left_join(., centroids) %>%
    rename(tox=PC1, toy=PC2)
}

# combine it all to make a spider plot
spiderplot <- function(rda_ordination, metadata){
  a <- ordfort(rda_ordination, metadata)
  b <- groupcentroids(a)
  c <- spiderize(a, b)
  
  ggplot(filter(a, Score=="sites")) +
    geom_point(aes(x=PC1, y=PC2, color=treatment)) +
    geom_segment(data = c, aes(x=fromx, xend = tox, 
                               y=fromy, yend = toy,
                               color=treatment)) +
    geom_point(data=b, aes(x=PC1, y=PC2, fill=treatment), shape=21, size=3) +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
}