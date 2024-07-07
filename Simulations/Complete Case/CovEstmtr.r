CovRegEstmtr = function(X) {
  n <- nrow(X)
  R <- ncol(X)
  CovInit <- t(X) %*% X / n
  CovAvg <- matrix(0, R-1, R-1)
  
  for(i in 1:(R-1)) {
    for(j in 1:(R-1)) {
      if(abs(i-j) > 1) {
        CovAvg[i,j] <- (CovInit[i,j] + CovInit[i+1,j] + CovInit[i,j+1] + CovInit[i+1,j+1]) / 4
      }
      if(i-j == 1) {
        CovAvg[i,j] <- (CovInit[i,j] + CovInit[i+1,j] + CovInit[i+1,j+1]) / 3
      }
      if(i-j == -1) {
        CovAvg[i,j] <- (CovInit[i,j] + CovInit[i,j+1] + CovInit[i+1,j+1]) / 3
      }
      if(abs(i-j) == 0) {
        CovAvg[i,j] <- CovInit[i,j+1]
      }
    }
  }
  CovAvg
}
CovSparseEstmtr <- function(X, M, r) {
  N <- nrow(X)
  R <- ncol(X)
  d <- R %/% M
  covM <- matrix(0, M, M)
  sumOP <- matrix(0, R, R)
  
  sumOP <- t(X) %*% X
  diag(sumOP) <- 0
  
  for(i in 1:M) {
    for(j in 1:M) {
      idx <- 1:d + (i-1)*d
      jdx <- 1:d + (j-1)*d
      covM[i, j] <- sum(sumOP[idx, jdx])
    }
  }
  covM <- (M^2 / (N*r*(r-1))) * covM
  kronecker(covM, matrix(rep(1, (R/M)^2), R/M, R/M))
}
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

#library(Rcpp)
#library(RcppArmadillo)
#sourceCpp("./EquivalentRcppFunctions.cpp")

#R <- 300
#grid <- seq(from = 0, to = 1, length = R)
#CovFn <- function(x, y) { min(x, y)}
#Cov <- outer(grid, grid, Vectorize(CovFn))
#N <- 10000
#r <- 10
#X <- MASS::mvrnorm(n = N, mu = rep(0, R), Sigma = Cov)
#for(j in 1:N) {
#  s <- sample(1:R, r)
#  X[j, -s] <- 0
#}

#CrossValidationError2(X, 2, 10, 100)
#CrossValidationError2Rcpp(X, 2, 15, 100)



