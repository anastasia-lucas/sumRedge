#' Simulate genotypes based off of MAF, penetrance, & encoding
#' @param maf MAF of SNP
#' @param eff biological action of SNP
#' @param cases number of cases
#' @param controls number of controls
#' @param model genetic model to generate, list of length 3
#' @param penetrance penetrance table, list of length 9
#' @param penbase baseline to use for penetrance tables, default = 0.25
#' @param pendiff difference between min and max probabilities in the penetrance table (default = 1-2*penbase)
#' @param seed seed

simulate_genoytypes <- function(maf, eff, cases, controls, model, penetrance, 
                                penbase = 0.25, pendiff, seed = 42){
  
  set.seed(seed)
  n <- cases + controls
  
  penbl <- penbase
  if(missing(pendiff)) pendiff <- 1-2*penbase
  
  pentable <- matrix(0, nrow = 3, ncol = 3)
  
  CC <- TRUE
  pen <- FALSE
  
  if(length(model) == 3){
    # auto-scale the coefficients
    coeff <- (model-min(model))/sum(model)
    coeff <- c(penetrance, pendiff * model)
  } else {
    stop("Coefficients in model must be 3 elements")
  }
  
  # Specify the penetrance table directly
  if(!missing(penetrance)){
    pen <- TRUE
    pentable <- t(matrix(as.numeric(penetrance), nrow = 3, ncol = 3))
    
    if(min(pentable) != 0 || max(pentable) != 1){
      warning("Penetrance table not scaled from 0 to 1, scaling")
      pentable = (pentable - min(pentable)) / (max(pentable) - min(pentable))
    }
    
    if(CC){
      if(penbase < 0 || penbase >= 1){
        stop("Penetrance baseline not in range [0,1), exiting")
      }
      
      if(pendf < 0 || penbl + pendf > 1){
        stop("Penetrance baseline and difference have nonsenical values (diff < 0 or baseline + diff > 1)")
      }
      
    }
    
    pentable <- penbl + pendf * pentable;
    
  }
  
  if(min(eff) != 0 & max(eff) != 1){
    warning("Effect for SNP not scaled from 0 to 1, scaling")
    eff = (eff - min(eff)) / (max(eff) - min(eff))
  }

  id_num <- 1
  
  store.snp <- data.frame()
  
  while (n_case + n_control > 0 && n_samp > 0){
    snp1 = sum(runif(2) < maf1)
    snp2 = sum(runif(2) < maf2)
    
    snpvec = as.numeric(c(snp1 > 1, snp1 > 0, snp2 > 1, snp2 > 0)) + 1
    
    if(pen){
      prob = pentable[snp1+1, snp2+1]
    } else {
      prob = coeff[1] + coeff[2]*eff1[snp1+1] + coeff[3]*eff2[snp2+1] + coeff[4]*eff1[snp1+1]*eff2[snp2+1]
    }
    
    print <- FALSE
    if(CC){
      status <- runif(1) < prob
      if (status == 1 & n_case > 0){
        n_case <- n_case - 1
        print <- TRUE
      } else if( status == 0 & n_control > 0){
        n_control <- n_control - 1
        print <- TRUE
      }
    } else {
      status <- prob + rnorm(1)
      n_samp <- n_samp - 1
      print <- TRUE
    }
    
    store.snp <- rbind(store.snp, cbind(id_num, status, t(as.data.frame(snpvec))))
    id_num <- id_num + 1
    
    
  }
  
  names(store.snp) <- c("IID", "Pheno", "SNP1_1", "SNP1_2", "SNP2_1", "SNP2_2")
  
  return(store.snp)
  
}





