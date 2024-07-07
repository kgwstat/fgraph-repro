# CovRegEstmtr = function(X) {
#   m <- nrow(X)
#   n <- ncol(X)
#   CovInit <- t(X) %*% X / m
#   CovAvg <- matrix(0, n-1, n-1)
#   
#   for(i in 1:(n-1)) {
#     for(j in 1:(n-1)) {
#       if(abs(i-j) > 1) {
#         CovAvg[i,j] <- (CovInit[i,j] + CovInit[i+1,j] + CovInit[i,j+1] + CovInit[i+1,j+1]) / 4
#       }
#       if(i-j == 1) {
#         CovAvg[i,j] <- (CovInit[i,j] + CovInit[i+1,j] + CovInit[i+1,j+1]) / 3
#       }
#       if(i-j == -1) {
#         CovAvg[i,j] <- (CovInit[i,j] + CovInit[i,j+1] + CovInit[i+1,j+1]) / 3
#       }
#       if(abs(i-j) == 0) {
#         CovAvg[i,j] <- CovInit[i,j+1]
#       }
#     }
#   }
#   CovAvg
# }
CovSparseEstmtr <- function(X, M, r) {
  N <- nrow(X)
  R <- ncol(X)
  d <- R %/% M
  covM <- matrix(0, M, M)
  sumOP <- matrix(0, R, R)

  for(j in 1:N) {
    sumOP <- sumOP + X[j, ] %*% t(X[j, ])
  }
  diag(sumOP) <- 0

  for(i in 1:M) {
    for(j in 1:M) {
      idx <- 1:d + (i-1)*d
      jdx <- 1:d + (j-1)*d
      covM[i, j] <- sum(sumOP[idx, jdx])
    }
  }
  covM <- (M^2 / (N*r*(r-1))) * covM
  s <- cut(1:R, breaks=M, labels=FALSE, include.lowest = TRUE)
  f <- function(i, j) {covM[i, j]}
  outer(s, s, Vectorize(f))
}
# CovSparseEstmtrM <- function(X, M, r) {
#   N <- nrow(X)
#   R <- ncol(X)
#   d <- R %/% M
#   covM <- matrix(0, M, M)
#   sumOP <- matrix(0, R, R)
#   
#   for(j in 1:N) {
#     sumOP <- sumOP + X[j, ] %*% t(X[j, ])
#   }
#   diag(sumOP) <- 0
#   
#   for(i in 1:M) {
#     for(j in 1:M) {
#       idx <- 1:d + (i-1)*d
#       jdx <- 1:d + (j-1)*d
#       covM[i, j] <- sum(sumOP[idx, jdx])
#     }
#   }
#   covM <- (M^2 / (N*r*(r-1))) * covM
#   covM
# }
# CrossValidationError <- function(X, M, r) {
#   N <- nrow(X)
#   R <- ncol(X)
#   d <- R %/% M
#   cov <- CovSparseEstmtr(X, M, r)
#   sum <- 0
#   for(l in 1:N) {
#     OPM <- matrix(0, M, M)
#     OP <- X[l, ] %*% t(X[l, ])
#     diag(OP) <- 0
#     for(i in 1:M) {
#       for(j in 1:M) {
#         idx <- 1:d + (i-1)*d
#         jdx <- 1:d + (j-1)*d
#         OPM[i, j] <- sum(OP[idx, jdx])
#       }
#     }
#     sum <- sum + sum(CovSparseEstmtr(X[-l, ], M, r) * OPM)
#   }
#   sum <- sum / (N*r*(r-1))
#   sum(cov^2)/(M^{2}) - 2*sum
# }
# CrossValidationError2 <- function(X, M, r, q) {
#   N <- nrow(X)
#   R <- ncol(X)
#   d <- R %/% M
#   cov <- CovSparseEstmtr(X, M, r)
#   sum <- 0
#   s <- sample(1:N, q)
#   
#   for(l in s) {
#     OPM <- matrix(0, M, M)
#     OP <- X[l, ] %*% t(X[l, ])
#     diag(OP) <- 0
#     for(i in 1:M) {
#       for(j in 1:M) {
#         idx <- 1:d + (i-1)*d
#         jdx <- 1:d + (j-1)*d
#         OPM[i, j] <- sum(OP[idx, jdx])
#       }
#     }
#     sum <- sum + sum(CovSparseEstmtr(X[-l, ], M, r) * OPM) 
#   }
#   sum <- sum / (q*r*(r-1))
#   sum(cov^2)/(M^{2}) - 2*sum
# }

# KCrossValidationError <- function(X, M, r, q = 5) {
#   N <- nrow(X)
#   R <- ncol(X)
#   d <- R %/% M
#   cov <- CovSparseEstmtrM(X, M, r)
#   total <- 0
# 
#   # Split into q parts
#   s <- split(1:N, cut(1:N, breaks=q, labels=FALSE, include.lowest = TRUE))
# 
#   for(l in s) {
#     OPM <- matrix(0, M, M)
#     OP <- X[l, ] %*% t(X[l, ])
#     diag(OP) <- 0
#     for(i in 1:M) {
#       for(j in 1:M) {
#         idx <- 1:d + (i-1)*d
#         jdx <- 1:d + (j-1)*d
#         OPM[i, j] <- sum(OP[idx, jdx])
#       }
#     }
#     total <- total + sum(CovSparseEstmtrM(X[-l, ], M, r) * OPM)
#   }
#   total <- total / (q*r*(r-1))
#   sum(cov^2)/(M^{2}) - 2*total
# }

#library(Rcpp)
#library(RcppArmadillo)
#sourceCpp("./EquivalentRcppFunctions.cpp")

# R <- 600
# grid <- seq(from = 0, to = 1, length = R)
# CovFn <- function(x, y) { min(x, y)}
# Cov <- outer(grid, grid, Vectorize(CovFn))
# N <- 10000
# r <- 10
# X <- MASS::mvrnorm(n = N, mu = rep(0, R), Sigma = Cov)
# for(j in 1:N) {
#  s <- sample(1:R, r)
#  X[j, -s] <- 0
# }
# Mat <- CovSparseEstmtr(X, 20, r) 
# norm(Mat - Cov, type = "F")/norm(Cov, type = "F")
# 
# #cut(1:100, breaks=7, labels=FALSE, include.lowest = TRUE)
# library(doParallel)
# cl <- makeCluster(6)
# registerDoParallel(cl)
# foreach(k = 7:15) %dopar% {
#   Mat <- CovSparseEstmtr(X, k, r) 
#   norm(Mat - Cov, type = "F")/norm(Cov, type = "F")
# }
# 
# 
# 
# CovSparseEstmtr(X, 15, r) %>% lattice::levelplot()
# 
# 
# KCrossValidationError(X, 10, 10, 5)
#CrossValidationError2Rcpp(X, 2, 15, 100)

