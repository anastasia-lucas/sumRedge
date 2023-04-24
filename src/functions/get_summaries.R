#' get case control ratio from string
#' @param filepath
get_cc <- function(filepath){
  sub <- strsplit(strsplit(strsplit(filepath, "param_")[[1]][2], "_")[[1]][1], "-")[[1]]
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

#' get penetrance from string
#' @param filepath
get_penetrance <- function(filepath){
  sub <- strsplit(strsplit(strsplit(test.dir, "param_")[[1]][2], "_")[[1]][2], "-")[[1]][1]
  return(sub)
}

#' load & format data
#' @param filepath
#' @return data.frame of genotype matrix and list of phenotypes
format_data <- function(filepath){
  raw <- as.data.frame(data.table::fread(file.path(filepath)))
  raw$PHENOTYPE <- raw$PHENOTYPE - 1 
  genos <- raw[, grep("disease", names(raw))]
  
  x <- list(geno = genos, pheno = raw$PHENOTYPE)
  
  return(x)
}

#' calculate edge weights
#' @param x genotype vector
#' @param y phenotypes
#' @return edge alpha values
calc_edge <- function(x, y){
  x <- data.frame(raw = x)
  x$het <- ifelse(x$raw == 1, 1, 0)
  x$homo.alt <- ifelse(x$raw == 2, 1, 0)

  mod <- summary(glm(y ~ x$het + x$homo.alt - 1))
  
  alpha <- mod$coefficients[1,1]/mod$coefficients[2,1]
  
  return(alpha)
  
}
 
#' get summary stats for simulated data
#' @param geno genotype as returned from format_data()
#' @param pheno phenotype as returned from format_data()
#' @return list of summary stats
get_summary_stats <- function(geno, pheno){
  
  models <- paste0("pheno ~ geno$", names(genos))
  
  regression <- lapply(models, function(x){
    res <- summary(glm(formula = x, family = "binomial"))$coefficients
    # beta, or, se, p
    res.list <- c(res[2,1], exp(res[2,1]), res[2,2], res[2,4])
    return(res.list)
  })
  
  x <- list()
  
  x$beta <- unlist(lapply(regression, function(x) x[1]))
  x$or <- unlist(lapply(regression, function(x) x[2]))
  x$se <- unlist(lapply(regression, function(x) x[3]))
  x$p <- unlist(lapply(regression, function(x) x[4]))
  
  return(x)
  
}
  
  
#' get summary stats for simulated data
#' @param beta genotype as returned from format_data()
#' @param pheno phenotype as returned from format_data()
#' @return list of summary stats
reconstruct_from_summary <- function(n.cases, n.ctrls, freq, or, se){
  
  gcounts <- list()
  for(i in 1:length(or)){
    gcounts[[i]] <- get_gcount(n.case = n.cases,
                               n.ctrl = n.ctrls,
                               or = or[i],
                               se = se[i],
                               freq = freq[i])
  }
  
  return(gcounts)
}

#' run regression on genotype matrix reconstruction
#' @param reconstructed.geno reconstructed genotypes as returned from reconstruct_from_summary()
#' @return list, regression results
rerun_regression <- function(reconstructed.geno){
  
  genos.new <- lapply(reconstructed.geno, construct_genotypes, expand = FALSE)
  
  regression <- lapply(genos.new, function(x){
    res <- tryCatch(summary(glm(pheno ~ snp, data = x))$coefficients, error = function(e) c(0,0,0,0))
    # beta, or, se, p
    res.list <- tryCatch(c(res[2,1], exp(res[2,1]), res[2,2], res[2,4]), error = function(e) c(0,0,0,0))
    return(res.list)
  })
  
  x <- list()
  
  x$beta <- unlist(lapply(regression, function(x) x[1]))
  x$or <- unlist(lapply(regression, function(x) x[2]))
  x$se <- unlist(lapply(regression, function(x) x[3]))
  x$p <- unlist(lapply(regression, function(x) x[4]))
  
  return(x)
  
}

#' get summary for 10K snp set
#' @param filepath simulated data location
#' @return summary comparison between original and reconstructed genotypes
get_summary_comparison <- function(filepath){
  
  print("Loading data...")
  
  data <- format_data(filepath)
  
  print("Calculating edges...")
  
  edges <- lapply(data$geno, calc_edge, y = data$pheno)
  
  print("Computing summary statistics...")
  
  reg.org <- get_summary_stats(geno = data$geno, pheno = data$pheno)
  
  print("Calculating allele frequencies...")
  
  freq.pop <- lapply(data$geno, get_freq2, data$pheno)
  
  freq.case <- lapply(data$geno, get_freq2, y = data$pheno, group = "case")
  
  freq.ctrl <- lapply(data$geno, get_freq2, y = data$pheno, group = "control")
  
  print("Reconstructing genotype vectors...")
  
  cc <- get_cc(filepath)
  
  gcounts <- reconstruct_from_summary(n.cases = cc$case, 
                                      n.ctrl = cc$ctrl,
                                      or = reg.org$or,
                                      se = reg.org$se,
                                      freq = unlist(freq.pop))
  
  print("Rerunning regressions....")
  
  reg.new <- rerun_regression(gcounts)
  
  result <- list(original = reg.org, reconstruction = reg.new, 
                 case.freq = freq.case, ctrl.freq = freq.ctrl,
                 gcounts = gcounts, edges = edges)
  
  
  print("Finished")
  
  return(result)
  
}


