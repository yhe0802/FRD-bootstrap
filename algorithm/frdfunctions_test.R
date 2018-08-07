library(testthat)
library(rdrobust)
source("frdfunctions.R")

# a simple DGP for testing
generate.data <- function(model.id, c) {
  x <- 2*rbeta(500, 2, 4) - 1
  uty <- matrix(rnorm(1000), ncol = 2)
  t <- as.numeric(ifelse(x < 0, uty[ , 1] < qnorm(0.5 - c/2), uty[ , 1] < qnorm(0.5 + c/2)))
  
  y <- uty[ , 2]*0.1295^2 +
    switch(model.id,
           ifelse(x < 0, 0.04*t + 1.27*x +  7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5,
                  0.04*t + 0.84*x -  3.00*x^2 +  7.99*x^3 -  9.01*x^4 + 3.56*x^5),
           ifelse(x < 0, -3.45*t + 2.30*x +  3.28*x^2 +  1.45*x^3 +  0.23*x^4 + 0.03*x^5,
                  -3.45*t + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5),
           ifelse(x < 0, 0.04*t + 1.27*x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5,
                  0.04*t + 0.84*x - 0.1*3.00*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5))
  return(data.frame(y = y, t = t, x = x))
}

test_that("Function `inputlist' works correctly.", {
  D <- generate.data(1, 0.8)
  
  # a simple version of function `inputlist'
  inputlist_simple <- function(y, x, h, b, p, q, kernel, residual) {
    
    # units within bandwidth max(h, b)
    iL <- x < 0 & x > -max(h[1], b[1])
    iR <- x >=0 & x <  max(h[2], b[2])
    ihL <- x[iL] > - h[1]
    ihR <- x[iR] < h[2]
    
    yL <- y[iL]
    yR <- y[iR]
    xL <- x[iL]
    xR <- x[iR]
    
    # matrix for the outer bootstrap
    # polynomials
    xL.poly <- poly(xL, q, raw = T) 
    xR.poly <- poly(xR, q, raw = T)
    # the design matrix
    XL <- cbind(1, xL.poly) 
    XR <- cbind(1, xR.poly)
    # coefficient maker
    WL <- diag(kweight(xL, 0, b[1], kernel))
    WR <- diag(kweight(xR, 0, b[2], kernel))
    coefL <- solve(t(XL) %*% WL %*% XL) %*% t(XL) %*% WL
    coefR <- solve(t(XR) %*% WR %*% XR) %*% t(XR) %*% WR
    
    # residual adjustment
    WXL <- sqrt(kweight(xL, 0, b[1], kernel)) * cbind(1, poly(xL, q, raw = T))
    WXR <- sqrt(kweight(xR, 0, b[2], kernel)) * cbind(1, poly(xR, q, raw = T))
    
    hL <- diag(WXL %*% solve(t(WXL) %*% WXL) %*% t(WXL))
    hR <- diag(WXR %*% solve(t(WXR) %*% WXR) %*% t(WXR))
    
    adjustL <- switch (residual,
                       HC0 = 1,
                       HC1 = sqrt(length(yL)/(length(yL) - q - 1)),
                       HC2 = 1/sqrt(1 - hL),
                       HC3 = 1/(1 - hL)
    )
    
    adjustR <- switch (residual,
                       HC0 = 1,
                       HC1 = sqrt(length(yR)/(length(yR) - q - 1)),
                       HC2 = 1/sqrt(1 - hR),
                       HC3 = 1/(1 - hR)
    )
    
    # matrix for the inner bootstrap (with "p" in notation)
    # polynomials
    xpL.poly <- poly(xL[ihL], p, raw = T) 
    xpR.poly <- poly(xR[ihR], p, raw = T)
    # the design matrix
    XpL <- cbind(1, xpL.poly) 
    XpR <- cbind(1, xpR.poly)
    # coefficient maker
    WLp <- diag(kweight(xL[ihL], 0, h[1], kernel))
    WRp <- diag(kweight(xR[ihR], 0, h[2], kernel))
    coefpL <- solve(t(XpL) %*% WLp %*% XpL) %*% t(XpL) %*% WLp
    coefpR <- solve(t(XpR) %*% WRp %*% XpR) %*% t(XpR) %*% WRp
    
    return(list(L=yL, R=yR, ihL=ihL, ihR=ihR, XL=XL, XR=XR, 
                coefL=coefL, coefR=coefR, 
                coefpL=coefpL, coefpR=coefpR,
                adjL=adjustL, adjR=adjustR))
  }
  
  # results from original version
  A <- inputlist(D$y, D$x, c(0.2, 0.2), c(0.3, 0.3), 1, 2, "tri", "HC3")
  
  # results from simple version
  B <- inputlist_simple(D$y, D$x, c(0.2, 0.2), c(0.3, 0.3), 1, 2, "tri", "HC3")

  # test fitted values
  expect_equal(A$XL %*% A$coefL %*% A$L, B$XL %*% B$coefL %*% B$L)
  expect_equal(A$XR %*% A$coefR %*% A$R, B$XR %*% B$coefR %*% B$R)
  
  # test intercept
  expect_equal(A$b0L %*% A$L, c(1, 0, 0) %*% B$coefL %*% B$L)
  expect_equal(A$b0R %*% A$R, c(1, 0, 0) %*% B$coefR %*% B$R)
  expect_equal(A$b0pL %*% A$L[A$ihL], c(1, 0) %*% B$coefpL %*% B$L[B$ihL])
  expect_equal(A$b0pR %*% A$R[A$ihR], c(1, 0) %*% B$coefpR %*% B$R[B$ihR])
  
  # test adjust factor
  expect_equal(A$adjL, B$adjL)
  expect_equal(A$adjR, B$adjR)
})

test_that("Function `frdboot' produces the same bias-uncorrected estimate as
          `rdrobust'.", {
  D <- generate.data(1, 0.8)
  h <- 0.2
  b <- 0.3
  
  # result from frdboot
  A <- frdboot(D$y, D$t, D$x, h = h, b = b, cluster = NULL,
               Nbc = 5, Nci = 5, kernel = "tri")
  
  # result from rdrobust
  B <- rdrobust(D$y, D$x, fuzzy = D$t, h = h, b = b, kernel = "tri")
  
  # test bias-uncorrected estimate
  expect_equal(A[1], B$coef[1])
})