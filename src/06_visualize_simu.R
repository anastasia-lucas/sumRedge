####################################################
#                      readme                      #  
####################################################

# visualize results

####################################################
#                        env                       #  
####################################################

rm(list = ls())

library(dplyr)
library(ggplot2)

set.seed(42)

####################################################
#                     variables                    #  
####################################################

# inputs
data.path <- file.path("snapshots", "summary-disease-all.txt")

# outputs


####################################################
#                     functions                    #  
####################################################

####################################################
#                       main                       #  
####################################################

data <- as.data.frame(data.table::fread(data.path)) %>%
  tidyr::pivot_longer(-c(dataset, edges.med), 
                      names_to = "metric", values_to = "value")

data$maf <- sub("simu-maf-", "", gsub("/.*", "", data$dataset))
data$cc <- sub("_.*", "", sub(".*param_", "", data$dataset))
data$or <- sub("\\.param.*", "", sub(".*or_", "", data$dataset))
data$prev <- sub("\\.raw", "", sub(".*_", "", data$dataset))

ggplot(data = data[data$metric %in% c("cor.freq.case", "cor.freq.ctrl") &
                     data$or != "1-10.0",], 
       aes(x = edges.med,
           y = value,
           shape = cc,
           color = prev)) +
  geom_point() +
  facet_grid(maf ~ or) +
  theme_bw()


