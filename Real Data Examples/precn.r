######################################################################
library(Matrix)
######################################################################
Diag <- function(Mat, p) {
    n <- nrow(Mat)
    ones <- function(k) {
        Matrix::Matrix(rep(1, k^2), k, k)
    } # returns k times k matrix of ones
    DiagofMat <- c(rep(n %/% p, p), n %% p) |>
                 lapply(ones) |> Matrix::bdiag()
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
    Param  <- Params |> 
              sapply(error, simplify = TRUE) |> which.min() |> 
              Params[argmin = _] # tuning
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
    D <- Matrix::Diagonal(x = EigenDecompn$values^{1 / 2})
    V <- EigenDecompn$vectors
    B <- V %*% D %*% t(V)
    C <- solve(B, solve(B, Mat  - Diag(Mat, p)) |> t()) |> t()
    R <- Matrix::Diagonal(n) + (C + t(C)) / 2 
    # symmetrize C and construct R
    return(R)
}
######################################################################
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
    conv2matrix <- \(x) { matrix(x, sqrt(length(x))) } 
    # constructs square matrices from vectors
    norm2       <- \(x) { norm(x, type = "2") / n } 
    # calculates operator norm
    MatofOps <- Mat |> split(spltr) |> lapply(conv2matrix) |> 
                sapply(norm2, simplify = TRUE) |> conv2matrix()
    return(MatofOps)
}
######################################################################

