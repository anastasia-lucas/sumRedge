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

dir <- file.path("simulated_data", "2pop")

pattern <- paste0("*raw")

all.files <- list.files(dir, pattern = glob2rx(pattern))

out.dir <- file.path("snapshots", "2pop")

####################################################
#                     functions                    #  
####################################################

source(file.path("src", "functions", "sumRedge.R"))

####################################################
#                       main                       #  
####################################################

for(filepath in all.files){
  print(filepath)
  
  result <- evaluate_sumRedge(file.path(dir, filepath))
  
  write.table(result, 
              file = file.path(out.dir, 
                               sub(".raw", "-summary.txt", filepath)),
              row.names = FALSE, 
              quote = FALSE,
              sep = "\t")
}


