library(rdrobust)
library(rdd)
rm(list = ls())

source("../algorithm/frdfunctions.R")

#### original data ####

final4 <- read.csv("final4.csv")
final4 <- final4[ , c("classize", "c_size", "avgverb", "avgmath")]
final4 <- final4[complete.cases(final4), ]

final5 <- read.csv("final5.csv")
final5 <- final5[ , c("classize", "c_size", "avgverb", "avgmath")]
final5 <- final5[complete.cases(final5), ]

#### Replicate Table 1 ####

summary(final4)
summary(final5)

#### RD plot ####

final4_sub <- final4[final4$c_size <= 80, ] # only use the first threshold at 40.
final5_sub <- final5[final5$c_size <= 80, ]

visualize_rd <- function(dta) {
  
  q1 <- ggplot(dta, aes(c_size, classize)) + 
    geom_point(size = 0.5) + 
    geom_smooth(data = subset(dta, c_size <= 40), method='lm',formula=y~poly(x,4)) +
    geom_smooth(data = subset(dta, c_size > 40), method='lm',formula=y~poly(x, 4)) +
    xlab("Enrollment") + ylab("Class size")
  
  q2 <- ggplot(dta, aes(c_size, avgverb)) + 
    geom_point(size = 0.5) + 
    geom_smooth(data = subset(dta, c_size <= 40), method='lm',formula=y~poly(x,4)) +
    geom_smooth(data = subset(dta, c_size > 40), method='lm',formula=y~poly(x, 4)) +
    xlab("Enrollment") + ylab("Verbal score")
  
  q3 <- ggplot(dta, aes(c_size, avgmath)) + 
    geom_point(size = 0.5) + 
    geom_smooth(data = subset(dta, c_size <= 40), method='lm',formula=y~poly(x,4)) +
    geom_smooth(data = subset(dta, c_size > 40), method='lm',formula=y~poly(x, 4)) +
    xlab("Enrollment") + ylab("Math score")
  
  pdf(paste(deparse(substitute(dta)), ".pdf", sep = ""), width=15, height=5)
  grid.arrange(q1, q2, q3, ncol = 3)
  dev.off()
}

visualize_rd(final4_sub)
visualize_rd(final5_sub)


#### RD estimate ####

estimate_rd <- function(dta, outcome, kernel = "tri") {
  
  # set the cutoff point to zero
  dta$c_size <- dta$c_size - 40.01
  
  # bandwidth choice
  bws <- rdbwselect(dta[ , outcome], dta$c_size, fuzzy = dta$classize,
                    kernel = kernel)$bws
  kernel2 <- ifelse(kernel == "uni", "rectangular",
                    ifelse(kernel == "tri", "triangular", "epanechnikov"))
  hIK <- IKbandwidth(dta$c_size, dta[ , outcome], kernel = kernel2)
  
  # CCT
  cct <- rdrobust(dta[ , outcome], dta$c_size, fuzzy = dta$classize,
                  kernel = kernel, h = bws[1:2], b = bws[3:4])
  cct <- c(cct$coef[1:2], cct$ci[3, ])
  
  # conventional
  if (outcome == "avgverb") {
    naive <- RDestimate(avgverb ~ c_size + classize, dta, bw = hIK, kernel = kernel2)
  } else {
    naive <- RDestimate(avgmath ~ c_size + classize, dta, bw = hIK, kernel = kernel2)
  }
  naive <- c(naive$est[1], naive$ci[1, ])
  
  # bootstrap
  boot <- frdboot(dta[ , outcome], dta$classize, dta$c_size,
                  kernel = kernel, h = bws[1:2], b = bws[3:4])
  return(list(bws = bws, hIK = hIK, boot = boot, cct = cct, naive = naive))
}

set.seed(798)
verb_4 <- estimate_rd(final4_sub, "avgverb")
math_4 <- estimate_rd(final4_sub, "avgmath")
verb_5 <- estimate_rd(final5_sub, "avgverb")
math_5 <- estimate_rd(final5_sub, "avgmath")