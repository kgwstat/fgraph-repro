Partition <- function(n, p) {
  q <- floor(n / p)  
  r <- n %% p        

  start <- 1 
  partition <- list() 
  for (i in 1:p) {
    size <- q + as.integer(i <= r)
    partition[[i]] <- start:(start + size - 1)
    start <- start + size
  }
  
  return(partition)
}

DiagonalPart <- function(M, partition) {
    n <- nrow(M)
    if (is.numeric(partition) && partition %% 1 == 0 && partition > 0) {
        partition <- Partition(n, partition)
    }
    
    PartitionOnes <- matrix(0, n, n)
    for (ix in partition) {
        PartitionOnes[ix, ix] <- 1
    }    
    return(PartitionOnes * M)
}

OpNormMatrix <- function(M, partition) {
    n <- nrow(M)
    if (is.numeric(partition) && partition %% 1 == 0 && partition > 0) {
        partition <- Partition(n, partition)
    }

    p <- length(partition)
    OpNormMat <- matrix(0, p, p)
    for (i in 1:p) {
        for (j in i:p) {
            OpNormMat[i, j] <- norm(M[partition[[i]], partition[[j]]], type = "2")
            OpNormMat[j, i] <- OpNormMat[i, j]
        }
    }
    return(OpNormMat / n)
}

CorrelationOpMatrix <- function(M, partition, epsilon = 1e-15, type = "norms") {
    n <- nrow(M)
    if (is.numeric(partition) && partition > 0 && partition %% 1 == 0) {
        partition <- Partition(n, partition)
    }

    EigenDecomposition <- eigen(DiagonalPart(M, partition) / n) 
    Eigenvalues <- EigenDecomposition$values
    Eigenvalues <- epsilon + (Eigenvalues + abs(Eigenvalues)) / 2
    
    D <- diag(x = Eigenvalues^{1 / 2})
    V <- EigenDecomposition$vectors      
    B <- V %*% D %*% t(V)
    C <- t(solve(B, t(solve(B, M  - DiagonalPart(M, partition)))))
    R <- diag(n) + (C + t(C)) / 2

    switch(type, norms = OpNormMatrix(R, partition), 
                 matrices = R)
}

PrecisionOpMatrix <- function(M, partition, epsilon = 1e-15, type = "norms") {
    n <- nrow(M)
    if (is.numeric(partition) && partition > 0 && partition %% 1 == 0) {
        partition <- Partition(n, partition)
    }

    R0 <- CorrelationOpMatrix(M, partition, epsilon, type = "matrices") - diag(n)
    P  <- Diagonal(n) - solve(diag(n) + R0 / n, R0) / n
    P  <- (P + t(P)) / 2
    
    switch(type, norms = OpNormMatrix(P, partition), 
                 matrices = P)
}






