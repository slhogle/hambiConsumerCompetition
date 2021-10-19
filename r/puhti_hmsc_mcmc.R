#!/usr/bin/env Rscript

# set custom library path for Hmsc
.libPaths(c("/projappl/project_2001175/rpackages-r-env-singul-4.0.2", .libPaths()))
libpath <- .libPaths()[1]

library(Hmsc)
library(ape)

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)!=3) {
  stop("Wrong number of arguments", call.=FALSE)
} 

model_name <- as.character(args[1])
samples    <- as.numeric(args[2])
thin       <- as.numeric(args[3])
nChains    <- 4

model <- readRDS(paste0("models_unfit/", as.character(model_name), 
                        ".rds"))

model <- sampleMcmc(model,
                    samples = samples,
                    thin=thin,
                    adaptNf=rep(ceiling(0.4*samples*thin),model$nr),
                    transient = ceiling(0.5*samples*thin),
                    nChains = nChains, 
                    nParallel = nChains)

filename <- paste("models_fit/", model_name, 
                  "_thin_",      as.character(thin),
                  "_samples_",   as.character(samples),
                  "_chains_",    as.character(nChains),
                  ".rds",        sep = "")

saveRDS(model, filename)
