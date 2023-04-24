###############################################################
#                          functions                          #
###############################################################

#' function to calculate group frequencies from summary stats
#' @param ncase number of cases
#' @param nctrl number of controls
#' @param or odds ratio
#' @param se standard error
#' @param freq frequency of alt allele
calc_group_freq <- function(ncase, nctrl, or, se, freq){
  
  x <- 2*ncase
  y <- 2*nctrl
  w <- se^2
  z <- or
  resFrq <- list()

  for(i in 0:49){
    w <- w*(1.001^i)
    tmp2 <- ((2*z*y*(1-z)-w*x*y*z)^2) - (4*y*z*(y*z+x)*((z-1)^2 + w*x*z))
    if(tmp2 >= 0) break
  }
  
  # no solution due to discriminant < 0
  # if this happens, it means the missing rate of this SNP is more than 20%
  # need to check the input quality
  if (tmp2 < 0) {
    resFrq$pCa <- 0.0;
    resFrq$pCon <- 0.0;
    resFrq$pPop <- 0.0;
  } else {
    tmp1 <- w*x*y*z - 2*y*z*(1-z)
    tmp2 <- tmp2^0.5
    tmp3 <- 2*((z-1)^2 + w*x*z)
    d1 <- (tmp1-tmp2)/tmp3
    c1 <- y-d1
    b1 <- x*d1/(z*y - z*d1 + d1)
    a1 <- x-b1
    frq1 <- c1/(c1+d1)
    d2 <- (tmp1 + tmp2)/tmp3
    c2 <- y-d2
    b2 <- (x*d2)/(z*y - z*d2 + d2)
    a2 <- x-b2
    frq2 <- c2/(c2+d2)
    flag1 <- 0
    flag2 <- 0
    
    if((a1 > 0) & (b1 > 0) & (c1 > 0) & (d1 > 0) ){
      flag1 <- 1; TmpRes1 <- list()
      TmpRes1$res11 <- a1
      TmpRes1$res12 <- b1
      TmpRes1$res21 <- b1
      TmpRes1$res22 <- d1
    }
    
    if((a2 > 0) & (b2 > 0) & (c2 > 0) & (d2 > 0)){
      flag2 <- 1; TmpRes2 <- list()
      TmpRes2$res11 <- a2
      TmpRes2$res12 <- b2
      TmpRes2$res21 <- c2
      TmpRes2$res22 <- d2
    }
    
    if(flag1 == 1 & flag2 == 0){
      res <- TmpRes1
    } else if(flag1 == 0 & flag2 == 1) {
      res <- TmpRes2
    } else if(flag1 == 1 & flag2 == 1){
      if(abs(freq - frq1) < abs(freq - frq2)){
        res <- TmpRes1
      } else {
        res <- TmpRes2
      }
    } else {
      res <- list()
      res$res11 <- 0.0
      res$res12 <- 0.0
      res$res21 <- 0.0
      res$res22 <- 0.0
    }
    
    resFrq$pCa <- res$res11/(res$res11 + res$res12)
    resFrq$pCon <- res$res21/(res$res21 + res$res22)
    resFrq$pPop <- (res$res11 + res$res21)/(res$res11 + res$res12 + res$res22)
    if(is.nan(resFrq$pCa)) resFrq$pCa <- 0
    if(is.nan(resFrq$pCon)) resFrq$pCon <- 0
    if(is.nan(resFrq$pPop)) resFrq$pPop <- 0
  }
  
  return(resFrq)  

}

#' function to calculate genotype matrix from group frequencies
#' 1:6 Case AA0/Aa1/aa2, Control AA0/Aa1/aa2
#' @param n.case number cases
#' @param n.ctrl number controls
#' @param or odds ratio
#' @param se standard error
#' @param freq frequency of alt allele
#' @param epsilon if maf < epsilon, no matrix will be constructed
get_gcount <- function(n.case, n.ctrl, or, se, freq, epsilon = 0.01){
  
  group.freq <- calc_group_freq(ncase = n.case, 
                                nctrl = n.ctrl,
                                or = or, se = se,
                                freq = freq)
  
  gcount <- list()
  if(abs(group.freq$pPop)  > epsilon){
    gcount$case[1] <- n.case*((1-group.freq$pCa)^2)
    gcount$case[2] <- n.case*2*group.freq$pCa*(1-group.freq$pCa)
    gcount$case[3] <- n.case*(group.freq$pCa^2)
    gcount$ctrl[1] <- n.ctrl*((1-group.freq$pCon)^2)
    gcount$ctrl[2] <- n.ctrl*2*group.freq$pCon*(1-group.freq$pCon)
    gcount$ctrl[3] <- n.ctrl*(group.freq$pCon^2)
  } else {
    gcount$case[1] <- 0; gcount$case[2] <- 0; gcount$case[3] <- 0
    gcount$ctrl[1] <- 0; gcount$ctrl[2] <- 0; gcount$ctrl[3] <- 0
  }
  # pop, case, ctrl
  gcount$freqs[1] <- group.freq$pPop
  gcount$freqs[2] <- group.freq$pCa
  gcount$freqs[3] <- group.freq$pCo
  
  return(gcount)
}

#' internal function to construct genotype vectors
#' @param g.count output from get_gcount()
#' @param expand expand into dummy (vs. additive), default TRUE
construct_genotypes <- function(g.count, expand = TRUE){
  # just make sure 0,1,2 refers to minor
  if(((g.count$case[1] + g.count$ctrl[1])*2) > ((g.count$case[3] + g.count$ctrl[3])*2)){
    minor <- 3; major <- 1
  } else{
    minor <- 1; major <- 3
  }
  
  if(expand){
    mat <- data.frame(pheno = c(rep(1, sum(round(g.count$case))),
                                rep(0, sum(round(g.count$ctrl)))),
                      het = c(rep(0, round(g.count$case[major])),
                              rep(1, round(g.count$case[2])),
                              rep(0, round(g.count$case[minor])),
                              rep(0, round(g.count$ctrl[major])),
                              rep(1, round(g.count$ctrl[2])),
                              rep(0, round(g.count$ctrl[minor]))),
                      homo.alt = c(rep(0, round(g.count$case[major])),
                                   rep(0, round(g.count$case[2])),
                                   rep(1, round(g.count$case[minor])),
                                   rep(0, round(g.count$ctrl[major])),
                                   rep(0, round(g.count$ctrl[2])),
                                   rep(1, round(g.count$ctrl[minor]))))
  } else {
    mat <- data.frame(pheno = c(rep(1, sum(round(g.count$case))),
                                rep(0, sum(round(g.count$ctrl)))),
                      snp = c(rep(0, round(g.count$case[major])),
                              rep(1, round(g.count$case[2])),
                              rep(2, round(g.count$case[minor])),
                              rep(0, round(g.count$ctrl[major])),
                              rep(1, round(g.count$ctrl[2])),
                              rep(2, round(g.count$ctrl[minor]))))
  }

  
  return(mat)
  
}

#' internal function to calculate EDGE weights
#' @param g.count output from get_gcount()
calc_edge <- function(g.count){
  
  mat <- construct_genotypes(g.count)
  
  mod <- summary(glm(pheno ~ het + homo.alt - 1, data = mat))
  
  alpha <- mod$coefficients[1,1]/mod$coefficients[2,1]
  
  return(alpha)
}

#' run edge regression from output of summary stats
#' @param n.case number cases
#' @param n.ctrl number controls
#' @param or odds ratio
#' @param se standard error
#' @param freq frequency of alt allele
#' @param epsilon if maf < epsilon, no matrix will be constructed
run_edge_regression <- function(n.case, n.ctrl, or, se, freq, epsilon){
  
  
  g.count <- get_gcount(n.case, n.ctrl, or, se, freq, epsilon)
  alpha <- calc_edge(g.count)
  
  
  # just make sure 0,alpha,2 refers to minor
  if(((g.count$case[1] + g.count$ctrl[1])*2) > ((g.count$case[3] + g.count$ctrl[3])*2)){
    minor <- 3; major <- 1
  } else{
    minor <- 1; major <- 3
  }
  
  mat <- data.frame(pheno = c(rep(1, sum(round(g.count$case))),
                              rep(0, sum(round(g.count$ctrl)))),
                    het = c(rep(0, round(g.count$case[major])),
                            rep(1, round(g.count$case[2])),
                            rep(0, round(g.count$case[minor])),
                            rep(0, round(g.count$ctrl[major])),
                            rep(1, round(g.count$ctrl[2])),
                            rep(0, round(g.count$ctrl[minor]))),
                    homo.alt = c(rep(0, round(g.count$case[major])),
                                 rep(0, round(g.count$case[2])),
                                 rep(1, round(g.count$case[minor])),
                                 rep(0, round(g.count$ctrl[major])),
                                 rep(0, round(g.count$ctrl[2])),
                                 rep(1, round(g.count$ctrl[minor]))))
  
  return(mat)
  
  
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


