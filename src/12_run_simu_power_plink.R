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

####################################################
#                     variables                    #  
####################################################

args <- commandArgs(trailingOnly = TRUE)

dir <- args[1]

filepath <- args[2]

out.dir <- file.path("snapshots", "2pop-power")

####################################################
#                     functions                    #  
####################################################

source(file.path("src", "functions", "sumRedge.R"))

####################################################
#                       main                       #  
####################################################

result <- evaluate_sumRedge(file.path(dir, filepath))

write.table(result,
            file = file.path(out.dir, sub(".raw", "-summary.txt", filepath)),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
