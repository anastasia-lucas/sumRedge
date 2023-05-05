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

# system(' cat <(head -1 snapshots/2pop/disease_dom_maf_0.85-0.9_or_1-0.5.param_2000-20000_0.01-0-hwe-summary.txt) <(cat snapshots/2pop/disease_dom_maf_0.85-0.9_or_1-* |grep -v "dataset") > snapshots/sumRedge-disease-all.txt')

####################################################
#                     variables                    #  
####################################################

# inputs
data.path <- file.path("snapshots", "sumRedge-disease-all.txt")

# outputs
out.path <- file.path("snapshots", "sumRedge_comparison.pdf")

# variables 
pal <- c("#0072B2", "#009E73", "#E69F00", "#D55E00")

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

p <- ggplot(data = data[data$metric %in% c("cor.freq.case") &
                     data$or != "1-10.0",], 
       aes(x = edges.med,
           y = value,
           shape = cc,
           color = prev)) +
  geom_point(size = 2) +
  facet_grid(maf ~ or) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = pal)

ggsave(p, filename = out.path, heigh = 5.5, width = 9)



