library(rdrobust)
library(rdd)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)
library(MASS)

source("frdfunctions.R")
try(source("slackrISE.R"), T)

#### simulation parameters ####
Nci <- 999
Nbc <- 500
Nsimu <- 5000
Ncore <- 38
a <- 0.05
commonseed <- 798

#### a wrapper for FRD simulations ####
frdsimu <- function(model.id, c, rho, cluster = NULL, 
                    heteroskedasticity = F, continuous = F, kernel = "tri") {
  
  # the RDD package uses different kernel name
  kernel2 <- ifelse(kernel == "uni", "rectangular",
                    ifelse(kernel == "tri", "triangular", "epanechnikov"))
  
  ## a function for each random draw
  eachdraw <- function() {
    
    dta <- gen.data(model.id, c, rho, cluster, heteroskedasticity, continuous)
    
    # for clustered data
    if (!is.null(cluster)) {
      bws <- rdbwselect(dta$y, dta$x, cluster = dta$g, kernel = kernel)$bws
      hIK <- IKbandwidth(dta$x, dta$y, kernel = kernel2)
      # wild bootstrap
      wildboot <- frdboot(dta$y, dta$t, dta$x, h = bws[1:2], b = bws[3:4], cluster = dta$g, a = a, 
                          Nbc = Nbc, Nci = Nci, kernel = kernel, residual = "HC3")[2:4]
      # CCT robust
      cctrobust <- rdrobust(dta$y, dta$x, fuzzy = dta$t, cluster = dta$g, h = bws[1:2], b = bws[3:4],
                            kernel = kernel, level = 100 - a*100)
      cctrobust <- c(cctrobust$coef[2], cctrobust$ci[3, ])
      # naive
      naive <- RDestimate(y ~ x + t, dta, bw = hIK, kernel = kernel2, cluster = dta$g)
      naive <- c(naive$est[1], naive$ci[1, ])
    } else {
      
      # for non-clustered data
      if (heteroskedasticity == F) {
        bws <- rdbwselect(dta$y, dta$x, kernel = kernel)$bws
        hIK <- IKbandwidth(dta$x, dta$y, kernel = kernel2)
      } else {
        bws <- rdbwselect(dta$y0, dta$x, kernel = kernel)$bws
        hIK <- IKbandwidth(dta$x, dta$y0, kernel = kernel2)
      }
      
      # wild bootstrap
      wildboot <- frdboot(dta$y, dta$t, dta$x, h = bws[1:2], b = bws[3:4], a = a, 
                          Nbc = Nbc, Nci = Nci, kernel = kernel, residual = "HC3")[2:4]
      # CCT robust
      cctrobust <- rdrobust(dta$y, dta$x, fuzzy = dta$t, h = bws[1:2], b = bws[3:4],
                            kernel = kernel, level = 100 - a*100)
      cctrobust <- c(cctrobust$coef[2], cctrobust$ci[3, ])
      # naive
      naive <- RDestimate(y ~ x + t, dta, bw = hIK, kernel = kernel2)
      naive <- c(naive$est[1], naive$ci[1, ])
    }
    
    # summary
    out <- c(wildboot, cctrobust, naive, bws[c(1, 3)], hIK)
    names(out) <- c("wild_t", "wild_l", "wild_u",
                    "cct_t", "cct_l", "cct_u",
                    "naive_t", "naive_l", "naive_u",
                    "h", "b", "hIK")
    
    return(out)
  }
  
  ## parallel computing
  cl <- makeCluster(Ncore)
  registerDoParallel(cl)
  export.obj <- c("Nci", "Nbc", "wild_values", "wild_weights", "a",
                  "frd_bc", "frd_dist", "frd_estimator", "frd_parameter",
                  "frdboot", "gen.data", "inputlist", "kweight", "ysd")
  result <- foreach(i=1:Nsimu, .combine="rbind", .packages=c("rdrobust", "MASS", "rdd"),
                    .export=export.obj, .inorder=F) %dorng% eachdraw()
  stopCluster(cl)
  
  ## clean the results
  tau    <- ifelse(model.id == 2, -3.45, 0.04)
  Nestimator <- 3
  table <- matrix(nrow = Nestimator, ncol = 8)
  colnames(table) <- c("Bias", "SD", "RMSE","EC(%)", "IL", "hCCT", "bCCT", "hIK")
  rownames(table) <- c("Wild bootstrap", "CCT robust", "Conventional")
  
  table[c(1, 2), 6] <- mean(result[ , 10])
  table[c(1, 2), 7] <- mean(result[ , 11])
  table[3, 8] <- mean(result[ , 12])
  
  for (i in 1:Nestimator) {
    table[i, 1] <- mean(result[ , (3*(i - 1) + 1)]) - tau
    table[i, 2] <- sd(result[ , (3*(i - 1) + 1)])
    table[i, 3] <- sqrt(mean((result[ , (3*(i - 1) + 1)] - tau)^2))
    table[i, 4] <- mean(result[ , (3*(i - 1) + 2)] <= tau &
                          result[ , (3*(i - 1) + 3)] >= tau)*100
    table[i, 5] <- mean(result[ , (3*(i - 1) + 3)] - result[ , (3*(i - 1) + 2)])
  }
  return(list(summary = table, raw = result))
}

#### organize results ####
jumps <- c(0.9)
kernels <- c("tri")
conti <- c(F)

auto.simu <- function(jumps, rhos, clusters, hetero, conti, kernels) {
  for (j in jumps) {
    for (r in rhos) {
      for (cl in clusters) {
        for (h in hetero) {
          for (con in conti) {
            for (k in kernels) {
              
              # naming
              if (cl == 0) {
                cl.n <- NULL
                cl.type <- "_nocl"
              } else {
                cl.n <- cl
                cl.type <- paste("_cl", cl, sep = "")
              }
              rho.type <- ifelse(r >= 0, paste("_rho", r*100, sep = ""),
                                 paste("_rhoN", -r*100, sep = ""))
              sd.type <- ifelse(h, "_hetero", "_homo")
              t.type <- ifelse(con, "_conti", "_binary")
              
              # simulations for m = 1, 2, 3
              set.seed(commonseed)
              simu.name <- paste(k, "_m", c(1, 2, 3), "_jump", j*100, rho.type,
                                 cl.type, sd.type, t.type, sep = "")
              assign(simu.name[1], frdsimu(1, j, r, cl.n, h, con, k), envir = .GlobalEnv)
              assign(simu.name[2], frdsimu(2, j, r, cl.n, h, con, k), envir = .GlobalEnv)
              assign(simu.name[3], frdsimu(3, j, r, cl.n, h, con, k), envir = .GlobalEnv)
              
              # combine results for m = 1, 2, 3
              all.name <- paste(k, "_jump", j*100, rho.type,
                                cl.type, sd.type, t.type, sep = "")
              all <- rbind(eval(parse(text = paste(simu.name[1], "$summary", sep = ""))),
                           eval(parse(text = paste(simu.name[2], "$summary", sep = ""))),
                           eval(parse(text = paste(simu.name[3], "$summary", sep = ""))))
              rownames(all) <- NULL
              all <- as.data.frame(all)
              all$DGP <- c("1", NA, NA,
                           "2", NA, NA,
                           "3", NA, NA)
              all$Method <- rep(c("Wild bootstrap", "CCT robust", "Conventional"), 3)
              all <- all[ , c("DGP", "Method", "Bias", "SD", "RMSE", 
                              "EC(%)", "IL", "hCCT", "bCCT", "hIK")]
            
              # send to slack and save as tex
              try(slackr(print(all.name)), T)
              try(slackr(print(all)), T)
              print(xtable(all, digits = c(0, 0, 3, 3, 3, 3, 1, 3, 3, 3, 3)), 
                    include.rownames = F, type = "latex", 
                    file = paste("output/",all.name, ".tex", sep = ""))
            }
          }
        }
      }
    }
  }
}

# basic simulation
auto.simu(jumps, rhos = 0, clusters = 0, hetero = F, conti, kernels)

# heteroskedasticity
auto.simu(jumps, rhos = 0, clusters = 0, hetero = T, conti, kernels)

# endogeneity
auto.simu(jumps, rhos = c(0.9, -0.9), clusters = 0, hetero = F, conti, kernels)

# clustered data
auto.simu(jumps, rhos = 0, clusters = c(5, 10, 25), hetero = F, conti, kernels)



#### additional speficication ####

#### warning: running the following code will override previous results

# in the function "frdsimu", when the data is heteroskedastic, y0 instead of y
# is used to calculate optimal bandwidth. y0 is outcome with homoskedastic error
# and y is outcome with heteroskedastic error. The intention is to make this two
# cases comparable (difference not driven by the bandwidth choice).

# the following function "frdsimu_2" is a copy of "frdsimu" except that y is used
# to calculate the optimal bandwidth in the heteroskedastic case.

frdsimu <- function(model.id, c, rho, cluster = NULL, 
                    heteroskedasticity = F, continuous = F, kernel = "tri") {
  
  # the RDD package uses different kernel name
  kernel2 <- ifelse(kernel == "uni", "rectangular",
                    ifelse(kernel == "tri", "triangular", "epanechnikov"))
  
  ## a function for each random draw
  eachdraw <- function() {
    
    dta <- gen.data(model.id, c, rho, cluster, heteroskedasticity, continuous)
    
    # for clustered data
    if (!is.null(cluster)) {
      bws <- rdbwselect(dta$y, dta$x, cluster = dta$g, kernel = kernel)$bws
      hIK <- IKbandwidth(dta$x, dta$y, kernel = kernel2)
      # wild bootstrap
      wildboot <- frdboot(dta$y, dta$t, dta$x, h = bws[1:2], b = bws[3:4], cluster = dta$g, a = a, 
                          Nbc = Nbc, Nci = Nci, kernel = kernel, residual = "HC3")[2:4]
      # CCT robust
      cctrobust <- rdrobust(dta$y, dta$x, fuzzy = dta$t, cluster = dta$g, h = bws[1:2], b = bws[3:4],
                            kernel = kernel, level = 100 - a*100)
      cctrobust <- c(cctrobust$coef[2], cctrobust$ci[3, ])
      # naive
      naive <- RDestimate(y ~ x + t, dta, bw = hIK, kernel = kernel2, cluster = dta$g)
      naive <- c(naive$est[1], naive$ci[1, ])
    } else {
      
      # for non-clustered data
      if (heteroskedasticity == F) {
        bws <- rdbwselect(dta$y, dta$x, kernel = kernel)$bws
        hIK <- IKbandwidth(dta$x, dta$y, kernel = kernel2)
      } else {
        bws <- rdbwselect(dta$y, dta$x, kernel = kernel)$bws
        hIK <- IKbandwidth(dta$x, dta$y, kernel = kernel2)
      }
      
      # wild bootstrap
      wildboot <- frdboot(dta$y, dta$t, dta$x, h = bws[1:2], b = bws[3:4], a = a, 
                          Nbc = Nbc, Nci = Nci, kernel = kernel, residual = "HC3")[2:4]
      # CCT robust
      cctrobust <- rdrobust(dta$y, dta$x, fuzzy = dta$t, h = bws[1:2], b = bws[3:4],
                            kernel = kernel, level = 100 - a*100)
      cctrobust <- c(cctrobust$coef[2], cctrobust$ci[3, ])
      # naive
      naive <- RDestimate(y ~ x + t, dta, bw = hIK, kernel = kernel2)
      naive <- c(naive$est[1], naive$ci[1, ])
    }
    
    # summary
    out <- c(wildboot, cctrobust, naive, bws[c(1, 3)], hIK)
    names(out) <- c("wild_t", "wild_l", "wild_u",
                    "cct_t", "cct_l", "cct_u",
                    "naive_t", "naive_l", "naive_u",
                    "h", "b", "hIK")
    
    return(out)
  }
  
  ## parallel computing
  cl <- makeCluster(Ncore)
  registerDoParallel(cl)
  export.obj <- c("Nci", "Nbc", "wild_values", "wild_weights", "a",
                  "frd_bc", "frd_dist", "frd_estimator", "frd_parameter",
                  "frdboot", "gen.data", "inputlist", "kweight", "ysd")
  result <- foreach(i=1:Nsimu, .combine="rbind", .packages=c("rdrobust", "MASS", "rdd"),
                    .export=export.obj, .inorder=F) %dorng% eachdraw()
  stopCluster(cl)
  
  ## clean the results
  tau    <- ifelse(model.id == 2, -3.45, 0.04)
  Nestimator <- 3
  table <- matrix(nrow = Nestimator, ncol = 8)
  colnames(table) <- c("Bias", "SD", "RMSE","EC(%)", "IL", "hCCT", "bCCT", "hIK")
  rownames(table) <- c("Wild bootstrap", "CCT robust", "Conventional")
  
  table[c(1, 2), 6] <- mean(result[ , 10])
  table[c(1, 2), 7] <- mean(result[ , 11])
  table[3, 8] <- mean(result[ , 12])
  
  for (i in 1:Nestimator) {
    table[i, 1] <- mean(result[ , (3*(i - 1) + 1)]) - tau
    table[i, 2] <- sd(result[ , (3*(i - 1) + 1)])
    table[i, 3] <- sqrt(mean((result[ , (3*(i - 1) + 1)] - tau)^2))
    table[i, 4] <- mean(result[ , (3*(i - 1) + 2)] <= tau &
                          result[ , (3*(i - 1) + 3)] >= tau)*100
    table[i, 5] <- mean(result[ , (3*(i - 1) + 3)] - result[ , (3*(i - 1) + 2)])
  }
  return(list(summary = table, raw = result))
}

# re-run the heteroskedastic case
auto.simu(jumps, rhos = 0, clusters = 0, hetero = T, conti, kernels)




