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

out.dir <- file.path("snapshots", "2pop")

####################################################
#                     functions                    #  
####################################################

source(file.path("src", "functions", "sumRedge.R"))

####################################################
#                       main                       #  
####################################################

for(seed in 1:100){
  
    print(seed)
    
    result <- evaluate_sumRedge_null(file.path(dir, filepath), seed)
    
    write.table(result, 
                file = file.path(out.dir, 
                                 paste0(sub(".raw", "-summary", filepath),
                                        "-", seed, ".txt")),
                row.names = FALSE, 
                quote = FALSE,
                sep = "\t")
    
}
  



