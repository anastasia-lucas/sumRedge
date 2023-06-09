source(file.path("src", "functions", "get_allele_counts.R"))
source(file.path("src", "functions", "get_summaries.R"))

#' load & format data
#' @param filepath
#' @return data.frame of genotype matrix and list of phenotypes
format_subsets <- function(filepath){
  set.seed(42)
  
  raw <- as.data.frame(data.table::fread(file.path(filepath)))
  raw <- raw[sample(1:nrow(raw)), ]
  raw$PHENOTYPE <- raw$PHENOTYPE - 1 
  
  genos <- raw[, 7:ncol(raw)]
  
  x1 <- list(geno = genos[1:11000,], pheno = raw$PHENOTYPE[1:11000])
  x2 <- list(geno = genos[11001:22000, ], pheno = raw$PHENOTYPE[11001:22000])
  
  x <- list(train = x1, test = x2)
  
  return(x)
}

#' run EDGE regression on genotype matrix reconstruction
#' @param geno genotype matrix
#' @param pheno phenotype vector
#' @param edges edge-derived alpha values as returned from calc_edges()
#' @return list, regression results
sumRedge <- function(geno, pheno, edges){
  
  genos.new <- lapply(geno, function(x){
    y <- data.frame(pheno = pheno,
                    snp = x)
  })
  
  genos.edge <- mapply(function(x, e){
    x$snp <- x$snp/2
    x$snp[x$snp == 0.5] <- e
    return(x)
  }, x = genos.new, e = edges, SIMPLIFY = FALSE)
  
  regression <- lapply(genos.edge, function(x){
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
  x$edge <- unname(unlist(edges))
  
  return(x)
  
}

#' get summary for 10K snp set
#' @param filepath simulated data location
#' @return summary comparison between original and reconstructed genotypes
evaluate_sumRedge <- function(filepath){
  
  print("Loading data...")
  
  set.seed(0)
  
  data <- format_subsets(filepath)
  
  # to speed up
  set.seed(42)

  print("Computing summary statistics...")
  
  reg.org <- get_summary_stats(geno = data$train$geno, 
                               pheno = data$train$pheno)
  
  freq.pop <- lapply(data$train$geno, get_freq2, data$train$pheno)
  
  print("Reconstructing genotype vectors...")
  
  gcounts <- reconstruct_from_summary(n.cases = length(data$train$pheno[data$train$pheno == 1]), 
                                      n.ctrl = length(data$train$pheno[data$train$pheno == 0]),
                                      or = reg.org$or,
                                      se = reg.org$se,
                                      freq = unname(unlist(freq.pop)))
  
  print("Calculating edges")
  
  edges.train <- lapply(data$train$geno, calc_edge, y = data$train$pheno)
  edges.test <- lapply(data$test$geno, calc_edge, y = data$test$pheno)
  
  print("Rerunning regressions...")
  
  reg.new <- sumRedge(data$test$geno, data$test$pheno, edges.train)
  
  reg.edge <- sumRedge(data$test$geno, data$test$pheno, edges.test)
  
  reg.add <- get_summary_stats(geno = data$test$geno, 
                               pheno = data$test$pheno)
  
  print("Concatenating results...")
  
  result <- data.frame(dataset = filepath,
                       sumRedge.edge = reg.new$edge,
                       sumRedge.beta = reg.new$beta,
                       sumRedge.or = reg.new$or,
                       sumRedge.se = reg.new$se,
                       sumRedge.p = reg.new$p,
                       edge.edge = reg.edge$edge,
                       edge.beta = reg.edge$beta,
                       edge.or = reg.edge$or,
                       edge.se = reg.edge$se,
                       edge.p = reg.edge$p,
                       add.beta = reg.add$beta,
                       add.or = reg.add$or,
                       add.se = reg.add$se,
                       add.p = reg.add$p)
  
  print("Finished")
  
  return(result)
  
}

#' get summary for 1K snp set
#' @param filepath simulated data location
#' @param seed seed for permuting phenotypes
#' @return summary comparison between original and reconstructed genotypes
evaluate_sumRedge_null <- function(filepath, seed){
  
  print("Loading data...")
  
  set.seed(42)
  
  data <- format_subsets(filepath)
  
  print("Computing summary statistics...")
  
  reg.org <- get_summary_stats(geno = data$train$geno, 
                               pheno = data$train$pheno)
  
  freq.pop <- lapply(data$train$geno, get_freq2, data$train$pheno)
  
  print("Reconstructing genotype vectors...")
  
  gcounts <- reconstruct_from_summary(n.cases = length(data$train$pheno[data$train$pheno == 1]), 
                                      n.ctrl = length(data$train$pheno[data$train$pheno == 0]),
                                      or = reg.org$or,
                                      se = reg.org$se,
                                      freq = unname(unlist(freq.pop)))
  
  print("Calculating edges")
  
  edges.train <- lapply(data$train$geno, calc_edge, y = data$train$pheno)
  
  print("Rerunning regressions...")
  
  set.seed(seed)
  
  y <- sample(data$test$pheno, size = length(data$test$pheno))
  
  reg.new <- sumRedge(data$test$geno, pheno = y, edges.train)
  
  edges.test <- lapply(data$test$geno, calc_edge, y = y)
  reg.edge <- sumRedge(data$test$geno, pheno = y, edges.test)
  
  reg.add <- get_summary_stats(geno = data$test$geno, pheno = y)
  
  # print("Running permutations...")
  
  print("Concatenating results...")
  
  result <- data.frame(dataset = filepath,
                       sumRedge.edge = reg.new$edge,
                       sumRedge.beta = reg.new$beta,
                       sumRedge.or = reg.new$or,
                       sumRedge.se = reg.new$se,
                       sumRedge.p = reg.new$p,
                       edge.edge = reg.edge$edge,
                       edge.beta = reg.edge$beta,
                       edge.or = reg.edge$or,
                       edge.se = reg.edge$se,
                       edge.p = reg.edge$p,
                       add.beta = reg.add$beta,
                       add.or = reg.add$or,
                       add.se = reg.add$se,
                       add.p = reg.add$p)
  
  print("Finished")
  
  return(result)
  
}
