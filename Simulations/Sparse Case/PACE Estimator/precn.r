######################################################################
library(Matrix)
library(magrittr)
######################################################################
Diag <- function(Mat, p) {
    n <- nrow(Mat)
    ones <- function(k) {
        Matrix::Matrix(rep(1, k^2), k, k)
    } # returns k times k matrix of ones
    DiagofMat <- c(rep(n %/% p, p), n %% p) %>%
                 lapply(ones) %>% Matrix::bdiag()
    return(Mat * DiagofMat)
} # returns diagonal part of operator matrix
######################################################################
Param <- function(Mat) {
    n <- nrow(Mat)
    error <- function(x) {
        Diff <- Mat %*% solve(x * Matrix::Diagonal(n) + Mat, Mat) - Mat
        norm(Diff, type = "2") / norm(Mat, type = "2")
    }
    error <- purrr::possibly(error, otherwise = NULL, quiet = TRUE)

    Params <- 10^{-15:0} 
    Param  <- Params %>% 
              sapply(error, simplify = TRUE) %>% which.min() %>% 
              Params[argmin = .] # tuning
# fine tuning
#   Params <- c((1:9) / 10, 1:9) * Param 
#   Param  <- Params %>% 
#             sapply(error, simplify = TRUE) %>% which.min() %>% 
#             Params[argmin = _] # fine tuning
    return(Param)
} # returns regularization parameter for matrix inversion
######################################################################
ParamFn <- function(Mat) {
    n <- nrow(Mat)
    error <- function(x) {
        Diff <- Mat %*% solve(x * Matrix::Diagonal(n) + Mat, Mat) - Mat
        norm(Diff, type = "2") / norm(Mat, type = "2")
    }
    error <- purrr::possibly(error, otherwise = NULL, quiet = TRUE)
    
    Params <- 10^{-15:0} 
    Param  <- Params %>% 
        sapply(error, simplify = TRUE) %>% which.min() %>% 
        Params[argmin = .] # tuning
    # fine tuning
    #   Params <- c((1:9) / 10, 1:9) * Param 
    #   Param  <- Params %>% 
    #             sapply(error, simplify = TRUE) %>% which.min() %>% 
    #             Params[argmin = _] # fine tuning
    return(Param)
} # returns regularization parameter for matrix inversion
######################################################################
Corr <- function(Mat, p) {
    n <- nrow(Mat)
    # parameter choice
    Param <- Param(Diag(Mat, p) / n)
    # calculating eigendecomposition
    # for calculating inverse square root
    EigenDecompn <- eigen(Param * Matrix::Diagonal(n) + Diag(Mat, p) / n)
    PosEigenvalues <- EigenDecompn$values - min(EigenDecompn$values) + Param
    # additional regularization
    D <- Matrix::Diagonal(x = PosEigenvalues^{1 / 2})
    V <- EigenDecompn$vectors
    B <- V %*% D %*% t(V)
    C <- solve(B, solve(B, Mat  - Diag(Mat, p)) %>% t()) %>% t()
    R <- Matrix::Diagonal(n) + (C + t(C)) / 2
    # symmetrize C and construct R
    return(R)
}
#####################################################################
Prec <- function(Mat, p) {
    n  <- nrow(Mat)
    R0 <- Corr(Mat, p) - Matrix::Diagonal(n)
    P  <- Matrix::Diagonal(n) - solve(Matrix::Diagonal(n) + R0 / n, R0) / n
    P  <- (P + t(P)) / 2
}
######################################################################
MatofOpNorms <- function(Mat, p) {
    n <- nrow(Mat)
    d <- nrow(Mat) %/% p
    spltr <- {
        I <- (row(Mat) - 1) %/% d
        J <- (col(Mat) - 1) %/% d
        1 + I + J * (max(I) + 1)
    } # to split Mat into blocks
    conv2matrix <- function(x) { matrix(x, sqrt(length(x))) } 
    # constructs square matrices from vectors
    norm2       <- function(x) { norm(x, type = "2") / n } 
    # calculates operator norm
    MatofOps <- Mat %>% split(spltr) %>% lapply(conv2matrix) %>% 
                sapply(norm2, simplify = TRUE) %>% conv2matrix()
    return(MatofOps)
}
######################################################################
PrecReg <- function(Mat, p, param = 1) {
    n <- nrow(Mat)
    # parameter choice
    Param <- Param(Diag(Mat, p) / n)
    # calculating eigendecomposition
    # for calculating inverse square root
    EigenDecompn <- eigen(Param * Matrix::Diagonal(n) + Diag(Mat, p) / n) 
    PosEigenvalues <- EigenDecompn$values + 10*param*abs(min(EigenDecompn$values))
    # additional regularization
    D <- Matrix::Diagonal(x = PosEigenvalues^{1 / 2})
    V <- EigenDecompn$vectors
    B <- V %*% D %*% t(V)
    C <- solve(B, solve(B, Mat  - Diag(Mat, p)) %>% t()) %>% t()
    R <- Matrix::Diagonal(n) + (C + t(C)) / 2
    n  <- nrow(Mat)
    R0 <- R - Matrix::Diagonal(n)
    P  <- Matrix::Diagonal(n) - solve(Matrix::Diagonal(n) + R0 / n, R0) / n
    P  <- (P + t(P)) / 2
}
######################################################################
PrecMu <- function(Mat, p, mu) {
    n  <- nrow(Mat)
    Param <- mu 
    
    EigenDecompn <- eigen(Diag(Mat, p) / n) 
    PosEigenvalues <- EigenDecompn$values
    PosEigenvalues <- Param + (PosEigenvalues + abs(PosEigenvalues)) / 2
    
    D <- diag(x = PosEigenvalues^{1 / 2})
    V <- EigenDecompn$vectors
    B <- V %*% D %*% t(V)
    C <- solve(B, solve(B, Mat  - Diag(Mat, p)) %>% t()) %>% t()
    R <- Matrix::Diagonal(n) + (C + t(C)) / 2
    
    R0 <- R - Matrix::Diagonal(n)
    P  <- Matrix::Diagonal(n) - solve(Matrix::Diagonal(n) + R0 / n, R0) / n
    P  <- (P + t(P)) / 2
}
######################################################################
CV <- function(lambda, X, p, q) {
    S <- split(1:nrow(X), 1:nrow(X) %% q)
    m <- ncol(X)
    total <- 0
    for (s in S) {
        A   <- t(X[-s,]) %*% X[-s,] / length(-s)
        B   <- t(X[ s,]) %*% X[ s,] / length( s)
        
        DA  <- Diag(A, p) / m
        DB  <- Diag(B, p) / m
        DAi <- lambda * diag(m) + DA
        
        residual <- DB - DB %*% solve(DAi, DB)
        total <- total + norm(residual, type = "F")^2
    }
    return(total/length(S))
}
######################################################################
CVReg <- function(lambda, X, p, q) {
    S <- split(1:nrow(X), 1:nrow(X) %% q)
    m <- ncol(X)-1
    total <- 0
    for (s in S) {
        A   <- CovRegEstmtr(X[-s,])
        B   <- CovRegEstmtr(X[ s,])
        
        DA  <- Diag(A, p) / m
        DB  <- Diag(B, p) / m
        DAi <- lambda * diag(m) + DA
        
        residual <- DB - DB %*% solve(DAi, DB)
        total <- total + norm(residual, type = "F")^2
    }
    return(total/length(S))
}
######################################################################
