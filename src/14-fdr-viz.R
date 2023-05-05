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

# system('cat <(head -1 snapshots/2pop-power/disease_dom_maf_0.85-0.9_or_1-0.5.param_2000-20000_0.01-0-hwe-summary.txt) <(cat snapshots/2pop-power/disease_dom_maf_0.85-0.9_or_1-* |grep -v "dataset") > snapshots/2pop-power/sumRedge-disease-all.txt')

####################################################
#                     variables                    #  
####################################################

# inputs
data.dir <- file.path("snapshots", "2pop-fdr")

# outputs
out.path <- file.path("snapshots", "sumRedge_comparison-fpr")

# variables 
pal <- c("#0072B2", "#009E73", "#E69F00", "#D55E00",
         "#999999", "#56B4E9", "#F0E442", "#CC79A7")

####################################################
#                     functions                    #  
####################################################

get_fpr <- function(filepath){

  fpr <- data.frame()
  
  for(i in 1:100){
    data <- data.table::fread(paste0(filepath, "-", i, ".txt"))
    fpr.i <- data.frame(
      add.05 = nrow(data[data$add.p < 0.05, ])/nrow(data),
      edge.05 = nrow(data[data$edge.p < 0.05, ])/nrow(data),
      sumRedge.05 = nrow(data[data$sumRedge.p < 0.05, ])/nrow(data),
      # n.tests = nrow(data),
      dataset = unique(data$dataset),
      seed = i
    )

   fpr <- rbind(fpr, fpr.i)
    
  }
  
  return(fpr)
}

####################################################
#                       main                       #  
####################################################

all.files <- unique(lapply(strsplit(list.files(data.dir), "-"), function(x) return(x[-8])))

filepaths <- list()
for(i in 1:length(all.files)){
  filepaths[i] <- file.path(data.dir, paste(all.files[[i]], collapse = "-"))
}

all.fpr <- lapply(filepaths, get_fpr)

fpr.df <- do.call("rbind", all.fpr)

fpr.df$maf <- "0.85-0.90"
fpr.df$cc <- sub("_.*", "", sub(".*param_", "", fpr.df$dataset))
fpr.df$or <- sub("\\.param.*", "", sub(".*or_", "", fpr.df$dataset))
fpr.df$prev <- sub("-0-hwe\\.raw", "", sub(".*_", "", fpr.df$dataset))

fpr.df <- fpr.df[fpr.df$or != "1-10.0", ]

# false positive rates
p <- fpr.df %>% 
  tidyr::pivot_longer(c(add.05, edge.05, sumRedge.05),
                      names_to = 'method',
                      values_to = 'FPR') %>%
     ggplot(aes(x = FPR,
                y = method,
                color = method)) +
  geom_vline(xintercept = 0.025, linetype = 2, color = "#999999") +
  geom_vline(xintercept = 0.075, linetype = 2, color = "#999999") +
  geom_boxplot() +
  facet_grid(or ~ prev) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = pal, guide = "none") +
  geom_vline(xintercept = 0.05, color = "#999999") +
  ggtitle("False positive rate of null models in 100 permuted datasets (p < 0.05)")

ggsave(p, filename = paste(out.path, "box.pdf", sep = "-"), height = 8, width = 12)


# edge values
all.data <- list()
j <- 1
for(i in 1:length(all.files)){
  filepath <- file.path(data.dir, paste(all.files[[i]], collapse = "-"))
  for(k in 1:100){
    data <- data.table::fread(paste0(filepath, "-", k, ".txt"))
    data$seed <- k
    all.data[[j]] <- data
    j <- j + 1
  }
}

all.df <- do.call("rbind", all.data)
all.df$maf <- "0.85-0.90"
all.df$cc <- sub("_.*", "", sub(".*param_", "", all.df$dataset))
all.df$or <- sub("\\.param.*", "", sub(".*or_", "", all.df$dataset))
all.df$prev <- sub("-0-hwe\\.raw", "", sub(".*_", "", all.df$dataset))
all.df <- all.df[all.df$or != "1-10.0", ]

p <- all.df %>% 
  tidyr::pivot_longer(c(sumRedge.edge, edge.edge),
                      names_to = 'method',
                      values_to = 'alpha') %>%
  filter(alpha < Inf & alpha > -Inf) %>%
  ggplot(aes(x = alpha,
             y = method,
             color = method)) +
  geom_vline(xintercept = 0.5, linetype = 2, color = "#999999") +
  geom_vline(xintercept = 1, color = "#999999") +
  geom_boxplot(outlier.size = .75, outlier.alpha = 0.5) +
  facet_grid(or ~ prev) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_log10() +
  ggtitle("Computed heterozygous alpha values across 100 null 1K SNP datasets")

ggsave(p, filename = paste(out.path, "edges.pdf", sep = "-"), height = 8, width = 12)





