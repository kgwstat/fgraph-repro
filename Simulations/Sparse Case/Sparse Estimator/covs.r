######################################################################
R <- 600
grid <- seq(from = 0, to = 1, length = R)
######################################################################
CovFn1 <- function(x, y) { exp(-abs(x - y)^2) }
CovFn2 <- function(x, y) { min(x, y) }
w <- 0.5
SmtFn <- function(x, y) {
    if ((x - y > 0) & (x - y < w)) {1} else {0} 
}
Smt <- outer(grid, grid, Vectorize(SmtFn))
TriFn <- function(a, x) {
    (abs(x) / a < 1) * (1 - abs(x) / a)
}
CovFn4 <- function(x, y) {
    0.8*TriFn(0.8, x-y) + 0.2*TriFn(1, x-y)
}
Q <- 10
a <- 0.3
CovFn5 <- function(x, y) {
    i <- 1 + floor(x * Q)
    j <- 1 + floor(y * Q)
    
    x <- x - i / Q
    y <- y - j / Q
    
    val <- (1 - x) * (1 - y) * a^{abs(i - j)} +
        (1 - x) * (  y  ) * a^{abs(i - j - 1)} +
        (  x  ) * (1 - y) * a^{abs(i + 1 - j)} +
        (  x  ) * (  y  ) * a^{abs(i - j)}
}

Cov1 <- outer(grid, grid, Vectorize(CovFn1))
Cov2 <- outer(grid, grid, Vectorize(CovFn2))
Cov3 <- t(Smt) %*% Cov2 %*% Smt / R^{2}
Cov4 <- outer(grid, grid, Vectorize(CovFn4))
Cov5 <- outer(grid, grid, Vectorize(CovFn5))
######################################################################
######################################################################
Covs <- list(Cov1, Cov2, Cov3, Cov4, Cov5)
rm(Cov1, Cov2, Cov3, Cov4, Cov5, Smt)
######################################################################

