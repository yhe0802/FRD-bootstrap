kweight <- rdrobust:::rdrobust_kweight

wild_values <- c(1 - sqrt(5), 1 + sqrt(5)) / 2
wild_weights <- c(sqrt(5) + 1, sqrt(5) - 1) / (2 * sqrt(5))

# this function applies to both outcome Y and treatment T.
# it create a list of matrix, which will be repeatedly used.
# notes on notation: L(left), R(right), i(0-1 indicator).
inputlist <- function(y, x, h, b, p, q, kernel, residual) {
  
  # units within bandwidth max(h, b)
  iL <- x < 0 & x > -max(h[1], b[1])
  iR <- x >=0 & x <  max(h[2], b[2])
  ihL <- x[iL] > - h[1]
  ihR <- x[iR] < h[2]
  
  # if b < h, the following vectors will include units beyond b.
  # it is necessary because we need fitted values and residuals beyond b but
  # within h. It will not affect estimation because their weights are zero.
  yL <- y[iL]
  yR <- y[iR]
  xL <- x[iL]
  xR <- x[iR]
  
  # matrix for the outer bootstrap
  # orthogonal polynomials
  xL.poly <- poly(xL, q) 
  xR.poly <- poly(xR, q)
  # the design matrix
  XL <- cbind(1, xL.poly) 
  XR <- cbind(1, xR.poly)
  KXL <- kweight(xL, 0, b[1], kernel) * XL
  KXR <- kweight(xR, 0, b[2], kernel) * XR
  # coefficient maker
  coefL <- solve(crossprod(XL, KXL), t(KXL)) 
  coefR <- solve(crossprod(XR, KXR), t(KXR))
  # the original intercept maker
  b0L <- c(1, predict(xL.poly, 0)) %*% coefL
  b0R <- c(1, predict(xR.poly, 0)) %*% coefR
  
  # residual adjustment
  KXL.sqrt <- sqrt(kweight(xL, 0, b[1], kernel)) * XL
  KXR.sqrt <- sqrt(kweight(xR, 0, b[2], kernel)) * XR
  hL <- diag(KXL.sqrt %*% solve(crossprod(XL, KXL), t(KXL.sqrt)))
  hR <- diag(KXR.sqrt %*% solve(crossprod(XR, KXR), t(KXR.sqrt)))
  
  adjustL <- switch (residual,
                      HC0 = 1,
                      # HC1 = sqrt(length(yL)/(length(yL) - q - 1)), wrong
                      HC2 = 1/sqrt(1 - hL),
                      HC3 = 1/(1 - hL)
  )
  
  adjustR <- switch (residual,
                      HC0 = 1,
                      # HC1 = sqrt(length(yR)/(length(yR) - q - 1)), wrong
                      HC2 = 1/sqrt(1 - hR),
                      HC3 = 1/(1 - hR)
  )
  
  # matrix for the inner bootstrap (with "p" in notation)
  # orthogonal polynomials
  xpL.poly <- poly(xL[ihL], p) 
  xpR.poly <- poly(xR[ihR], p)
  # the design matrix
  XpL <- cbind(1, xpL.poly) 
  XpR <- cbind(1, xpR.poly)
  KXpL <- kweight(xL[ihL], 0, h[1], kernel) * XpL
  KXpR <- kweight(xR[ihR], 0, h[2], kernel) * XpR
  # coefficient maker
  coefpL <- solve(crossprod(XpL, KXpL), t(KXpL)) 
  coefpR <- solve(crossprod(XpR, KXpR), t(KXpR))
  # the original intercept maker
  b0pL <- c(1, predict(xpL.poly, 0)) %*% coefpL
  b0pR <- c(1, predict(xpR.poly, 0)) %*% coefpR
  
  return(list(L=yL, R=yR, ihL=ihL, ihR=ihR, XL=XL, XR=XR, 
              coefL=coefL, coefR=coefR, b0L=b0L, b0R=b0R, 
              b0pL=b0pL, b0pR=b0pR, adjL=adjustL, adjR=adjustR))
}

# this function finds the bootstrap parameter.
# i.e., intercept from fitting a qth order polynomial.
frd_parameter <- function(Yinput, Tinput) {
  (Yinput$b0R %*% Yinput$R - Yinput$b0L %*% Yinput$L)/
  (Tinput$b0R %*% Tinput$R - Tinput$b0L %*% Tinput$L)
}

# this function finds the FRD estimate.
# i.e., intercept from fitting a pth order polynomial.
frd_estimator <- function(Yinput, Tinput) {
  (Yinput$b0pR %*% Yinput$R[Yinput$ihR] - Yinput$b0pL %*% Yinput$L[Yinput$ihL])/
  (Tinput$b0pR %*% Tinput$R[Tinput$ihR] - Tinput$b0pL %*% Tinput$L[Tinput$ihL])
}

# this function finds the bias-corrected estimate.
frd_bc <- function(Yinput, Tinput, gen.wild, Nbc) {
  parameter <- frd_parameter(Yinput, Tinput)
  estimate <- frd_estimator(Yinput, Tinput)
  
  yL.fit <- Yinput$XL %*% (Yinput$coefL %*% Yinput$L)
  yR.fit <- Yinput$XR %*% (Yinput$coefR %*% Yinput$R)
  
  yL.res <- (Yinput$L - yL.fit) * Yinput$adjL
  yR.res <- (Yinput$R - yR.fit) * Yinput$adjR
  
  tL.fit <- Tinput$XL %*% (Tinput$coefL %*% Tinput$L)
  tR.fit <- Tinput$XR %*% (Tinput$coefR %*% Tinput$R)
  
  tL.res <- (Tinput$L - tL.fit) * Tinput$adjL
  tR.res <- (Tinput$R - tR.fit) * Tinput$adjR

  boot <- replicate(Nbc, {
    wild.e <- gen.wild()
    Yinput$L <- yL.fit + yL.res * wild.e$wild.e.L
    Yinput$R <- yR.fit + yR.res * wild.e$wild.e.R
    Tinput$L <- tL.fit + tL.res * wild.e$wild.e.L
    Tinput$R <- tR.fit + tR.res * wild.e$wild.e.R
    frd_estimator(Yinput, Tinput)
  })
  
  return(estimate - mean(boot) + parameter)
}

# this function finds the bootstrap distribution of bias-corrected estimate.
frd_dist <- function(Yinput, Tinput, gen.wild, Nbc, Nci) {
  yL.fit <- Yinput$XL %*% (Yinput$coefL %*% Yinput$L)
  yR.fit <- Yinput$XR %*% (Yinput$coefR %*% Yinput$R)
  
  yL.res <- (Yinput$L - yL.fit) * Yinput$adjL
  yR.res <- (Yinput$R - yR.fit) * Yinput$adjR
  
  tL.fit <- Tinput$XL %*% (Tinput$coefL %*% Tinput$L)
  tR.fit <- Tinput$XR %*% (Tinput$coefR %*% Tinput$R)
  
  tL.res <- (Tinput$L - tL.fit) * Tinput$adjL
  tR.res <- (Tinput$R - tR.fit) * Tinput$adjR
  
  boot <- replicate(Nci, {
    wild.e <- gen.wild()
    Yinput$L <- yL.fit + yL.res * wild.e$wild.e.L
    Yinput$R <- yR.fit + yR.res * wild.e$wild.e.R
    Tinput$L <- tL.fit + tL.res * wild.e$wild.e.L
    Tinput$R <- tR.fit + tR.res * wild.e$wild.e.R
    frd_bc(Yinput, Tinput, gen.wild, Nbc)
  })
  
  return(boot)
}

# a wrapper for FRD bootstrap
frdboot <- function(y, t, x, h, b, cluster = NULL, 
                    a = 0.05, Nbc = 500, Nci = 999, p = 1, q = 2, 
                    kernel = "tri", residual = "HC3"){
  
  # if the input bandwidth is a scalar, use it for both sides
  if (length(h) == 1) h <- c(h, h)
  if (length(b) == 1) b <- c(b, b)

  # lists of input objects
  Yinput <- inputlist(y, x, h, b, p, q, kernel, residual)
  Tinput <- inputlist(t, x, h, b, p, q, kernel, residual)

  # a function generating wild bootstrap values
  i.hb <- (x > -max(h[1], b[1]) & x < max(h[1], b[1]))
  i.hb.L <- x[i.hb] < 0
  i.hb.R <- !i.hb.L
  
  if (!is.null(cluster)) {
    cluster <- cluster[i.hb]
    
    gen.wild <- function() {
      e <- vector(length = length(cluster))
      for (i in unique(cluster)) {
        e[cluster == i] <- sample(wild_values, 1, prob = wild_weights)
      }
      return(list(wild.e.L = e[i.hb.L], wild.e.R = e[i.hb.R]))
    }
  } else {
    n.L <- sum(i.hb.L)
    n.R <- sum(i.hb.R)
    
    gen.wild <- function() {
      return(list(wild.e.L = sample(wild_values, n.L, T, wild_weights),
                  wild.e.R = sample(wild_values, n.R, T, wild_weights)))
    }
  }
  
  # estimation
  parameter <- frd_parameter(Yinput, Tinput)
  tau <- frd_estimator(Yinput, Tinput)
  taubc <- frd_bc(Yinput, Tinput, gen.wild, Nbc)
  taubcboot <- frd_dist(Yinput, Tinput, gen.wild, Nbc, Nci)
  ci <- taubc - quantile(taubcboot, c(1 - a/2, a/2)) + parameter
  
  return(c(tau, taubc, ci))
}

# a function to link the running variable to standard error in outcome
ysd <- function(x) 0.1295 + 9*x^2

# c is the jump in mean treatment.
# rho is the correlation coefficient between ut and uy.
gen.data <- function(model.id, c, rho = 0, cluster = NULL,
                     heteroskedasticity = F, continuous = F, Nobs = 1000) {
  
  x <- 2*rbeta(Nobs, 2, 4) - 1
  
  # firstly create joint standard normal ut and uy, with correlation rho.
  # ut is used to create treatment with a jump equal to c.
  # uy is rescaled and added to expected y.
  uty <- mvrnorm(Nobs, c(0, 0), matrix(c(1, rho, rho, 1), 2), empirical = T)
  
  # the treatment variable
  if (continuous == F) {
    t <- as.numeric(ifelse(x < 0, uty[ , 1] < qnorm(0.5 - c/2), 
                           uty[ , 1] < qnorm(0.5 + c/2)))
  } else {
    t <- ifelse(x < 0, uty[ , 1], uty[ , 1] + c)
  }
  
  # expected outcome
  Ey <- switch(model.id,
               ifelse(x < 0, 0.04*t + 1.27*x +  7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5,
                      0.04*t + 0.84*x -  3.00*x^2 +  7.99*x^3 -  9.01*x^4 + 3.56*x^5),
               ifelse(x < 0, -3.45*t + 2.30*x +  3.28*x^2 +  1.45*x^3 +  0.23*x^4 + 0.03*x^5,
                      -3.45*t + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5),
               ifelse(x < 0, 0.04*t + 1.27*x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5,
                      0.04*t + 0.84*x - 0.1*3.00*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5))
  
  # in the case of clustered data, the option "heteroskedasticity" will be ignored
  # notice that "cluster" is the number of clusters in each side.
  if (!is.null(cluster)) {
    g <- sample.int(cluster, Nobs, T)
    g[x >= 0] <- g[x >= 0] + cluster
    e.g <- vector(length = Nobs)
    for (i in 1:(2*cluster)) e.g[g == i] <- rnorm(1, 0, 0.1295/sqrt(2))
    return(data.frame(y = Ey + e.g + uty[ , 2]*0.1295/sqrt(2),
                      t = t, x = x, g = g))
    
  }
  
  # if it's not clustered data
  if (heteroskedasticity == F) {
    return(data.frame(y = Ey + uty[ , 2]*0.1295, t = t, x = x))
  } else {
    # in order to compare the cases with and without heteroskedasticity,
    # the same bandwidth is required (otherwise both different bandwidth and
    # different variance structure will affect the results).
    # so, both y0 (no heteroskedasticity) and y (with heteroskedasticity) will be 
    # generated. y0 will be used to obtain bandwidth.
    
    return(data.frame(y0 = Ey + uty[ , 2]*0.1295, 
                      y = Ey + uty[ , 2]*ysd(x),
                      t = t, x = x))
  }
}