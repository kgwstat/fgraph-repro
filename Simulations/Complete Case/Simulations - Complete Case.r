######################################################################
library(dplyr)
library(tidyr)
library(Matrix)
library(purrr)
library(MASS)
library(DescTools)
library(doParallel)
######################################################################
source("./covs.r")
source("./metrics.r")
source("./precn.r")
######################################################################
# R = 600
# p = 20, 30, 40
# N = 50, 100, 200
######################################################################
load("./CovDetail.RData")
AUCs <- function(j, p, N, Covs, OmegaTrue) {
    R <- 600
    
    # Starting cluster
    cl <- makeCluster(25)
    registerDoParallel(cl)
    clusterEvalQ(cl, {library(Matrix); source("./metrics.r"); source("./precn.r")})

    # Choosing Tuning Parameter
    set.seed(1238487)
    X     <- MASS::mvrnorm(n = N, mu = rep(0, R), Sigma = Covs[[j]])
    vals  <- 10^{-seq(0, 15, length.out = 15)} * norm(Diag(Covs[[j]]/R, p), type = "2")
    CVErr <- parSapply(cl, vals, function(y) {CV(y, X, p, 5)})
    Param <- vals[which.min(CVErr)]
    
    message(c("Param-Done: (j, p, N) = ( ", as.character(j), 
              ", ", as.character(p),
              ", ", as.character(N), ")"))
    
    # Calculating AUCs
    AUCs <- foreach(k = 1:100, .combine = 'c', .packages = c("DescTools")) %dopar% {
        set.seed(1238487 + k^2)
        X <- MASS::mvrnorm(n = N, 
                           mu = rep(0, R), 
                           Sigma = Covs[[j]])
        Cov.Estimate <- t(X) %*% X / N
        Prec.Estimate <- PrecMu(Cov.Estimate, p, Param)
        Omega.True <- OmegaTrue[[j]][[(p/10)]]
        AreaUnderCurveTrapezoid(Prec.Estimate, Omega.True)
    }
    
    # Stopping cluster
    stopCluster(cl)
    
    message(c("AUC-Done: (j, p, N) = ( ", as.character(j), 
                                ", ", as.character(p),
                                ", ", as.character(N), ")"))
    return(AUCs)
}
######################################################################
######################################################################
CovID       <- data.frame(j = c(1, 2, 3, 4, 5))
Resolution  <- data.frame(p = c(20, 30, 40))
NoofSamples <- data.frame(N = c(50, 100, 200))
Specs       <- tidyr::crossing(CovID, Resolution, NoofSamples)
######################################################################
######################################################################
AUCData <- Specs %>% 
    dplyr::mutate(AUCs = purrr::pmap(list(x = j, y = p, z = N), function(x, y, z) {AUCs(x, y, z, Covs, OmegaTrue)}))
save(AUCData, file="AUCDataComplete.RData")

AUCSummary <- AUCData %>% 
    dplyr::mutate(Median = sapply(AUCs, median)) %>%
    dplyr::mutate(Median = round(Median, 2)) %>%
    dplyr::mutate(MAD = sapply(AUCs, mad)) %>%
    dplyr::mutate(MAD = round(MAD, 2)) %>%
    dplyr::mutate(Summary = paste(Median, MAD, sep = "\u00B1")) %>%
    subset(select = -c(AUCs, Median, MAD)) %>%
    tidyr::pivot_wider(names_from = N, values_from = Summary)
save(AUCSummary, file="AUCSummaryComplete.RData")
######################################################################
######################################################################
