#!/usr/bin/env Rscript

# set custom library path for Hmsc
.libPaths(c("/projappl/project_2001175/rpackages-r-env-singul-4.0.2", .libPaths()))
libpath <- .libPaths()[1]

library(Hmsc)
library(ape)

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)!=1) {
  stop("Wrong number of arguments", call.=FALSE)
} 

model_name <- as.character(args[1])
model <- readRDS(paste0("models_fit/", as.character(model_name)))

# number of cores cannot be greater than number of chains
partition <- createPartition(model, nfolds = 5)
CFpreds <- computePredictedValues(model, partition=partition, nParallel=4) 

filename <- paste0("models_assess_fit/", as.character(model_name))
saveRDS(CFpreds, filename)