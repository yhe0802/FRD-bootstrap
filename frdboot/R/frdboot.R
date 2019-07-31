#' frdboot
#'
#' A function to calculate bias-corrected estimator and its confidence interval in fuzzy regression discontinuity design.
#'
#' @param y the outcome variable.
#' @param t the treatment variable.
#' @param x the running variable.
#' @param h the bandwidth for estimation.
#' @param b the bandwidth for bias correction.
#' @param cluster the cluster variable. Default value is NULL.
#' @param a the significance level. Default value is 0.05.
#' @param Nbc the number of bootstrap replications for bias correction. Default value is 500.
#' @param Nci the number of bootstrap replications for confidence interval. Default value is 999.
#' @param p the polynomial order for estimation. Default value is 1.
#' @param q the polynomial order for bias correction. Default value is 2.
#' @param kernel the kernel function: epa, uni or tri. Default value is "tri".
#' @param residual the residual adjustment used in wild bootstrap: HC0, HC2 or HC3. Default value is "HC3".
#' @return The original estimator, the bias-corrected estimator and the robust confidence interval.
#' @export
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
