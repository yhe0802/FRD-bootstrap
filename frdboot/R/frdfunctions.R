# calculate weight
kweight <- function (X, c, h, kernel) {
  u = (X - c)/h
  if (kernel == "epanechnikov" | kernel == "epa") {
    w = (0.75 * (1 - u^2) * (abs(u) <= 1))/h
  }
  else if (kernel == "uniform" | kernel == "uni") {
    w = (0.5 * (abs(u) <= 1))/h
  }
  else {
    w = ((1 - abs(u)) * (abs(u) <= 1))/h
  }
  return(w)
}

# distribution of the random variable used in wild bootstrap
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
