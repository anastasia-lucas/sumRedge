####################################################
#                      readme                      #  
####################################################

# check that genotype reconstruction is working

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

  
test.dir <- file.path("simulated_data/add_maf_0.01-0.02_or_0.5-0.5.param_1000-1000_0.01-0.raw")

out.dir <- file.path("simulated_data-additive_results/")

####################################################
#                     functions                    #  
####################################################

source(file.path("src", "functions", "get_allele_counts.R"))

get_cc <- function(filepath){
  sub <- strsplit(strsplit(strsplit(test.dir, "param_")[[1]][2], "_")[[1]][1], "-")[[1]]
  cc <- list()
  cc$case <- as.integer(sub[1])
  cc$ctrl <- as.integer(sub[1])
  return(cc)
}

#' get frequency of allele coded as 2
#' @param x vector of snps coded 0,1,2
get_freq2 <- function(x){
  a2f <- (length(x[x == 2])*2 + length(x[x == 1]))/(2*length(x))
  return(a2f)
}

####################################################
#                       main                       #  
####################################################

raw <- data.table::fread(file.path(test.dir))
raw$PHENOTYPE <- raw$PHENOTYPE - 1 
genos <- raw[, 7:ncol(raw)]

######### original regression 

regression <- lapply(genos, function(x){
  res <- summary(glm(raw$PHENOTYPE ~ x))$coefficients
  # beta, or, se, p
  res.list <- c(res[2,1], exp(res[2,1]), res[2,2], res[2,4])
  return(res.list)
})

beta.org <- unlist(lapply(regression, function(x) x[1]))
or.org <- unlist(lapply(regression, function(x) x[2]))
se.org <- unlist(lapply(regression, function(x) x[3]))
p.org <- unlist(lapply(regression, function(x) x[4]))

data.frame(or = or.org, p = p.org) %>%
  ggplot(aes(y = -log10(p), x = or)) +
  geom_point() +
  theme_classic()

######### reconstruct from summary

n.cases <- rep(get_cc(temp.dir)$case, length(or.org)) 
n.ctrls <- rep(get_cc(temp.dir)$ctrl, length(or.org))
freq.org <- unlist(lapply(genos, get_freq2))

all.res <- list()
for(i in 1:length(or.org)){
  all.res[[i]] <- get_gcount(n.case = n.cases[i],
                             n.ctrl = n.ctrls[i],
                             or = or.org[i],
                             se = se.org[i],
                             freq = freq.org[i])
}

######### rerun regression

genos.new <- lapply(all.res, construct_genotypes, expand = FALSE)

regression <- lapply(genos.new, function(x){
  res <- tryCatch(summary(glm(pheno ~ snp, data = x))$coefficients, error = function(e) c(0,0,0,0))
  # beta, or, se, p
  res.list <- tryCatch(c(res[2,1], exp(res[2,1]), res[2,2], res[2,4]), error = function(e) c(0,0,0,0))
  return(res.list)
})

beta.new <- unlist(lapply(regression, function(x) x[1]))
or.new <- unlist(lapply(regression, function(x) x[2]))
se.new <- unlist(lapply(regression, function(x) x[3]))
p.new <- unlist(lapply(regression, function(x) x[4]))

data.frame(or.org = or.org, or.new = or.new) %>% 
  ggplot(aes(x = or.org, y = or.new)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(min(min(or.org), min(or.new)), 
       max(max(or.org), max(or.new))) + 
  xlim(min(min(or.org), min(or.new)), 
       max(max(or.org), max(or.new))) + 
  theme_classic()


data.frame(beta.org = beta.org, beta.new = beta.new) %>% 
  ggplot(aes(x = beta.org, y = beta.new)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(min(min(beta.org), min(beta.new)), 
       max(max(beta.org), max(beta.new))) + 
  xlim(min(min(beta.org), min(beta.new)), 
       max(max(beta.org), max(beta.new))) + 
  theme_classic()


data.frame(p.org = -log10(p.org), p.new = -log10(p.new)) %>% 
  ggplot(aes(x = p.org, y = p.new)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() 

