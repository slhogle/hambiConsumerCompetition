parsemodname <- function(data){
  data %>%
    mutate(thin=str_extract(model1, "(?<=thin_)\\d+"),
           samples=str_extract(model1, "(?<=samples_)\\d+"),
           chains=str_extract(model1, "(?<=chains_)\\d+"),
           phase=ifelse(str_detect(model1, "qeq"), "quasi-equilibrium", "sorting"),
           distribution=ifelse(str_detect(model1, "ma"), "normal", "probit"),
           model=case_when(str_detect(model1, "envi") ~ "no_random_effects",
                            str_detect(model1, "time") ~ "no_fixed_effects",
                            str_detect(model1, "full") ~ "full")) %>%
    select(-model1)
}

loadhmsc <- function(myname){
  return(readr::read_rds(here::here("data", "JSDM_fit", myname)))
}

loadhmsccoda <- function(myname){
  mymodel <- loadhmsc(myname)
  return(convertToCodaObject(mymodel, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F)))
}

psrf_ess <- function(myname, myparam){
  mpost      <- loadhmsccoda(myname)
  if (myparam == "beta") {
    psrf <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
    ess  <- effectiveSize(mpost$Beta)
    return(list(tibble(parameter=myparam,
                       term=row.names(psrf), 
                       p.est=psrf[,1],
                       ess=ess,
                       model1=myname) %>% parsemodname()))
  } else if (myparam == "gamma") {
    psrf  <- gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
    ess  <- effectiveSize(mpost$Gamma)
    return(list(tibble(parameter=myparam,
                       term=row.names(psrf), 
                       p.est=psrf[,1],
                       ess=ess,
                       model1=myname) %>% parsemodname()))
  } else if (myparam == "rho") {
    psrf  <- gelman.diag(mpost$Rho, multivariate=FALSE)$psrf
    ess  <- effectiveSize(mpost$Rho)
    return(list(tibble(parameter=myparam,
                       term=NA,
                       p.est=psrf[,1],
                       ess=ess,
                       model1=myname) %>% parsemodname()))
  } else if(myparam == "omega" & str_detect(myname, "full")) {
    tmp  <- mpost$Omega[[1]]
    psrf <- gelman.diag(tmp, multivariate=FALSE)$psrf
    ess  <- effectiveSize(tmp)
    return(list(tibble(parameter=myparam,
                       term=row.names(psrf),
                       p.est=psrf[,1],
                       ess=ess,
                       model1=myname) %>% parsemodname()))
    }
}



R2fit <- function(myname){
  mymodel   <- loadhmsc(myname)
  preds <- computePredictedValues(mymodel)
  MF    <- evaluateModelFit(hM=mymodel, predY=preds)
  if(str_detect(myname, "_mp")){
    R2   <- MF$TjurR2
    AUC  <- MF$AUC
    RMSE <- MF$RMSE
    type <- "probit"
  } else {
    R2   <- MF$R2
    RMSE <- MF$RMSE
    type <- "gaussian"
    AUC <- NA
  }
  return(list(tibble(species=mymodel$spNames,
                         model1=myname,
                         type=type,
                         R2=as.numeric(R2),
                         RMSE=as.numeric(RMSE),
                         AUC=AUC) %>% parsemodname()))
}

WAICfit <- function(myname, mythin){
  mymodel      <- loadhmsc(myname)
  WAIC_species <- computeWAIC(mymodel, byColumn = T)
  WAIC_full    <- computeWAIC(mymodel, byColumn = F)
  
  return(list(tibble(strainID=mymodel$spNames,
                         model1=myname,
                         WAIC_species=WAIC_species,
                         WAIC_full=WAIC_full) %>% parsemodname()))
}

cv5fit <- function(myname){

  mymodel <- readr::read_rds(here::here("data", "JSDM_fit", myname))
  mypreds <- readr::read_rds(here::here("data", "JSDM_assess_fit", myname))
  
  MF <- evaluateModelFit(hM=mymodel, predY=mypreds)
  
  if(str_detect(myname, "_mp")){
    R2   <- MF$TjurR2
    AUC  <- MF$AUC
    RMSE <- MF$RMSE
    type <- "probit"
  } else {
    R2   <- MF$R2
    RMSE <- MF$RMSE
    type <- "gaussian"
    AUC <- NA
  }
  return(list(tibble(species=mymodel$spNames,
                         model1=myname,
                         type=type,
                         R2=as.numeric(R2),
                         RMSE=as.numeric(RMSE),
                         AUC=AUC) %>% parsemodname()))
}


vpformat <- function(name) {
  #m <- readr::read_rds(here::here("JSDM_fit", name))
  m  <- loadhmsc(name)
  vp <- computeVariancePartitioning(m)
  df <- data.frame(vp$vals) %>% as_tibble(rownames = "variable") %>%
    pivot_longer(starts_with("HAMBI"), names_to="species", values_to="fvar") %>%
    mutate(species=factor(str_replace(string = species, pattern = "\\.", replacement = "-"))) %>%
    mutate(variable=factor(variable, levels=c("Random: day", "Random: microcosmID", "effort", "day",
                                              "cil", "day:cil", "wrm", "day:wrm", "cil:wrm", "day:cil:wrm")),
           model1=name) %>% parsemodname()
  return(list(df))
}

vptraitformat <- function(name) {
  m  <- loadhmsc(name)
  #m <- readRDS(here::here("models_fit", name))
  vp <- computeVariancePartitioning(m)
  df <- data.frame(R2T_beta=vp$R2T$Beta) %>% as_tibble(rownames = "variable") %>%
    mutate(variable=factor(variable, levels=c("(Intercept)", "effort", "day",
                                              "day:cil1", "day:cil2", "day:wrm1",
                                              "day:cil1:wrm1", "day:cil2:wrm1",
                                              "cil1", "cil2", "wrm1",
                                              "cil1:wrm1", "cil2:wrm1")),
           model1=name,
           R2T_y=vp$R2T$Y) %>% parsemodname()
  return(list(df))
}

betaformat <- function(name, level) {
  #m <- readRDS(here::here("models_fit", name))
  m  <- loadhmsc(name)
  pb <- getPostEstimate(m, parName="Beta")
  
  pos <- pb$support
  neg <- pb$supportNeg
  
  df.s <- data.frame(as.matrix(pos))
  df.s[pos < neg] <- -neg[pos < neg]
  rownames(df.s) <- m$covNames
  
  betas <- as_tibble(df.s, rownames = "variable") %>%
    pivot_longer(-variable, names_to="species", values_to="beta") %>%
    mutate(variable=str_replace_all(variable, "present", "")) %>%
    mutate(species=factor(str_replace(string = species, pattern = "\\.", replacement = "-"))) %>%
    mutate(variable=factor(variable, levels=c("(Intercept)", "effort", "day",
                                              "day:cil1", "day:cil2", "day:wrm1",
                                              "day:cil1:wrm1", "day:cil2:wrm1",
                                              "cil1", "cil2", "wrm1",
                                              "cil1:wrm1", "cil2:wrm1"))) %>%
    mutate(sign=ifelse(abs(beta) >= level, round(beta), 0),
           model1=name) %>% parsemodname()
  
  return(list(betas))
}


gammaformat <- function(name, level) {
  #m <- readRDS(here::here("models_fit", name))
  m  <- loadhmsc(name)
  pg <- getPostEstimate(m, parName="Gamma")
  
  pos <- pg$support
  neg <- pg$supportNeg
  
  df.s <- data.frame(as.matrix(pos))
  df.s[pos < neg] <- -neg[pos < neg]
  rownames(df.s) <- m$covNames
  
  gammas <- as_tibble(df.s, rownames = "variable") %>%
    pivot_longer(-variable, names_to="trait", values_to="gamma") %>%
    mutate(variable=str_replace_all(variable, "present", "")) %>%
    mutate(trait=case_when(trait=="X1" ~ "Intercept",
                              trait=="X2" ~ "D",
                              trait=="X3" ~ "N_carbon_ecoplate",
                              trait=="X4" ~ "r_specific_growth_rate",
                              trait=="X5" ~ "biofilm_formation")) %>%
    mutate(variable=factor(variable, levels=c("(Intercept)", "effort", "day",
                                              "day:cil1", "day:cil2", "day:wrm1",
                                              "day:cil1:wrm1", "day:cil2:wrm1",
                                              "cil1", "cil2", "wrm1",
                                              "cil1:wrm1", "cil2:wrm1")),
           trait=factor(trait, levels=c("Intercept", "D", "N_carbon_ecoplate", 
                                              "r_specific_growth_rate", "biofilm_formation"))) %>%
    mutate(sign=ifelse(abs(gamma) >= level, round(gamma), 0),
           model1=name) %>% parsemodname()
    
    return(list(gammas))
}

traitgradient <- function(name){
  m  <- loadhmsc(name)
  #m <- readRDS(here::here("models_fit", name))
  
  grad.wrm0 <- constructGradient(m, focalVariable = "cil", non.focalVariables = list("wrm"=list(3,0)))
  grad.wrm1 <- constructGradient(m, focalVariable = "cil", non.focalVariables = list("wrm"=list(3,1)))
  
  py.wrm0 <- predict(m, XData=grad.wrm0$XDataNew, studyDesign=grad.wrm0$studyDesignNew,
                     ranLevels=grad.wrm0$rLNew, expected=T)
  py.wrm1 <- predict(m, XData=grad.wrm1$XDataNew, studyDesign=grad.wrm1$studyDesignNew,
                     ranLevels=grad.wrm1$rLNew, expected=T)
  
  traitwrm0 <- list()
  for (i in 2:5) {
    t <- plotGradient(m, grad.wrm0, pred=py.wrm0, measure="T", index=i, showData = F)$data %>% 
      as_tibble() %>%
      rename(cil=xx)
    traitwrm0[[i]] <- t
  }
  
  traitwrm0 <- bind_rows(traitwrm0, .id = "id") %>%
    mutate(trait=case_when(id==1 ~ "D",
                           id==2 ~ "N_carbon_ecoplate",
                           id==3 ~ "r_specific_growth_rate",
                           id==4 ~ "biofilm_formation")) %>%
    mutate(wrm=0)
  
  traitwrm1 <- list()
  for (i in 2:5) {
    t <- plotGradient(m, grad.wrm1, pred=py.wrm1, measure="T", index=i, showData = F)$data %>% 
      as_tibble() %>%
      rename(cil=xx)
    traitwrm1[[i]] <- t
  }
  
  traitwrm1 <- bind_rows(traitwrm1, .id = "id") %>%
    mutate(trait=case_when(id==1 ~ "D",
                           id==2 ~ "N_carbon_ecoplate",
                           id==3 ~ "r_specific_growth_rate",
                           id==4 ~ "biofilm_formation")) %>%
    mutate(wrm=1)
  
  bind_rows(traitwrm0, traitwrm1) %>% mutate(model1=name) %>% parsemodname()
}

myvarplot <- function(d, cols){
  ggplot(data=d) +
    geom_bar(aes(x=species, y=tvar, fill=variable), stat="identity") + 
    labs(x="", y="", fill="") + 
    coord_flip() + 
    scale_y_reverse(limits=c(1,0)) + 
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

mybetaplot <- function(d){
  ggplot(d) +
    geom_tile(aes(x=variable, y=species, fill=sign, alpha=abs(sign))) +
    labs(x="", y="", fill="") +
    scale_fill_carto_c(palette = "Fall", guide=FALSE) + 
    scale_alpha_continuous(range=c(0.35,1), guide=FALSE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

mygammaplot <- function(d) {
  ggplot(d) +
    geom_tile(aes(x=variable, y=trait, fill=sign, alpha=abs(sign))) +
    labs(x="", y="", fill="") +
    scale_fill_carto_c(palette = "Fall", guide=FALSE) + 
    scale_alpha_continuous(range=c(0.35,1), guide=FALSE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank())
}