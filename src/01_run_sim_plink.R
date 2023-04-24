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

  
test.dir <- file.path("simulated_data")

all.files <- 

out.dir <- file.path("snapshots/simulated_data_results/dominant.txt")

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
#' @param vector of phenotypes
#' @param group 'case', 'control', 'population'
get_freq2 <- function(x, y, group = "population"){
  if(group == 'case') {x <- x[y == 1]}
  if(group == 'control') {x <- x[y == 0]}
  a2f <- (length(x[x == 2])*2 + length(x[x == 1]))/(2*length(x))
  return(a2f)
}

####################################################
#                       main                       #  
####################################################

raw <- as.data.frame(data.table::fread(file.path(test.dir)))
raw$PHENOTYPE <- raw$PHENOTYPE - 1 
genos <- raw[, grep("disease", names(raw))]

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
case.org <- unlist(lapply(genos, get_freq2, group = "case", y = raw$PHENOTYPE))
ctrl.org <- unlist(lapply(genos, get_freq2, group = "control", y = raw$PHENOTYPE))

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


######### make figures

# grid columns = metric
# grid rows = maf
# x-axis = OR with penetrance grouped
# color CC ratio
# difference

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

######## compare actual het and homozygous counts
table(genos.new[[1]]$pheno, genos.new[[1]]$snp)
tab <- table(raw$PHENOTYPE, genos$disease_0_D)

# p(case | ref)
p.ref <- (tab[2,1]/(tab[1,1] + tab[2,1])); or.ref <- p.ref / (1 - p.ref)
# p(case | ref)
p.alt <- (tab[2,3]/(tab[1,3] + tab[2,3])); or.alt <- p.alt/ (1 - p.alt)

p.ref/p.alt




