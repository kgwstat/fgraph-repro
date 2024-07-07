######################################################################
library(dplyr)
library(tidyr)
library(Matrix)
library(purrr)
library(MASS)
library(DescTools)
library(doParallel)
library("fdapace")
######################################################################
source("./covs.r")
source("./metrics.r")
source("./precn.r")
source("./CovEstmtr.r")
######################################################################
# R = 600
# p = 10, 20, 30, 40
# N = 500
# noise = 0, 0.01, 0.1
######################################################################
load("./CovDetail.RData")
AUCs <- function(j, p, N, noise, Covs, OmegaTrue) {
    R <- 600
    r <- 5
    
    # Starting cluster
    cl <- makeCluster(25)
    registerDoParallel(cl)
    clusterEvalQ(cl, {library(Matrix); source("./metrics.r"); source("./precn.r"); source("./CovEstmtr.r")})
    
    Param <- 10^{-9} * norm(Diag(Covs[[j]]/R, p), type = "2")
    
    # Calculating AUCs
    AUCs <- foreach(k = 1:100, .combine = 'c', .packages = c("DescTools")) %dopar% {
        set.seed(1238487 + k^2)
        X <- MASS::mvrnorm(n  = N,
                           mu = rep(0, R),
                           Sigma = Covs[[j]] + noise*diag(R))
        for(i in 1:N) {
            s <- sample(1:R, r)
            X[i, -s] <- NA
        }
        Cov.Estimate  <- CovPACE(X)
        Prec.Estimate <- PrecMu(Cov.Estimate, p, Param)
        Omega.True    <- OmegaTrue[[j]][[(p/10)]]
        AreaUnderCurveTrapezoid(Prec.Estimate, Omega.True) #Using trapezoidal rule
    }
    
    # Stopping cluster
    stopCluster(cl)
    
    message(c("AUC-Done: (j, p, N, noise) = ( ", as.character(j), 
              ", ", as.character(p),
              ", ", as.character(N),", ", as.character(noise), ")"))
    
    return(AUCs)
}
######################################################################
######################################################################
CovID       <- data.frame(j = c(1, 2, 3, 4, 5))
Resolution  <- data.frame(p = c(10, 20, 30, 40))
NoofSamples <- data.frame(N = c(500))
NoiseLevel  <- data.frame(noise = c(0, 0.01, 0.1))
Specs       <- tidyr::crossing(CovID, Resolution, NoofSamples, NoiseLevel)
######################################################################
######################################################################
AUCData <- Specs %>% 
    dplyr::mutate(AUCs = purrr::pmap(list(x = j, y = p, z = N, w = noise), function(x, y, z, w) {AUCs(x, y, z, w, Covs, OmegaTrue)}))
save(AUCData, file="AUCDataSparse.RData")

AUCSummary <- AUCData %>% 
    dplyr::mutate(Median = sapply(AUCs, median)) %>%
    dplyr::mutate(Median = round(Median, 2)) %>%
    dplyr::mutate(MAD = sapply(AUCs, mad)) %>%
    dplyr::mutate(MAD = round(MAD, 2)) %>%
    dplyr::mutate(Summary = paste(Median, MAD, sep = "\u00B1")) %>%
    subset(select = -c(AUCs, Median, MAD)) %>%
    tidyr::pivot_wider(names_from = N, values_from = Summary)
save(AUCSummary, file="AUCSummarySparse.RData")
######################################################################
######################################################################
