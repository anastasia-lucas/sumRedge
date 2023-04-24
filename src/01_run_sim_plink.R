#!/usr/bin/env Rscript

####################################################
#                      readme                      #  
####################################################

# run simulations and get average for 1K SNPs

####################################################
#                        env                       #  
####################################################

rm(list = ls())

library(dplyr)

set.seed(42)

####################################################
#                     variables                    #  
####################################################

args <- commandArgs(trailingOnly = TRUE)

args <- c("simu-maf-0.25-0.35", "1-0.5.param")

dir <- args[1]

pattern <- paste0("*", args[2], "*raw")

all.files <- list.files(dir, pattern = glob2rx(pattern))

out.dir <- file.path("snapshots", "summaries")

####################################################
#                     functions                    #  
####################################################

source(file.path("src", "functions", "get_allele_counts.R"))
source(file.path("src", "functions", "get_summaries.R"))

####################################################
#                       main                       #  
####################################################

for(filepath in all.files){
  print(filepath)
  
  result <- get_summary_comparison(file.path(dir, filepath))
  
  write.table(result, 
              file = file.path(out.dir, 
                               sub(".raw", "-summary.txt", filepath)),
              row.names = FALSE, 
              quote = FALSE,
              sep = "\t")
}


