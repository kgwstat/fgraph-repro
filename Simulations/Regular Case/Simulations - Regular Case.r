######################################################################
library(dplyr)
library(tidyr)
library(Matrix)
library(purrr)
library(MASS)
library(DescTools)
library(doParallel)
######################################################################
#source("./covs.r")
source("./covsReg.r")
source("./metrics.r")
source("./precn.r")
######################################################################
# R = 601
# p = 20, 30, 40
# N = 50, 100, 200
# noise = 0, 0.01, 0.1
######################################################################
load("./CovDetail.RData")
AUCs <- function(j, p, N, noise, CovsReg, OmegaTrue) {
    R <- 601
    
    # Starting cluster
    cl <- makeCluster(25)
    registerDoParallel(cl)
    clusterEvalQ(cl, {library(Matrix); source("./metrics.r"); source("./precn.r"); source("./CovEstmtr.r")})

    # Choosing Tuning Parameter
    set.seed(1238487)
    X     <- MASS::mvrnorm(n = N, mu = rep(0, R), Sigma = CovsReg[[j]])
    vals  <- 10^{-seq(-1,2, length.out = 15)} * norm(Diag(CovsReg[[j]]/R, p), type = "2")
    CVErr <- parSapply(cl, vals, function(y) {CVReg(y, X, p, 5)})
    Param <- vals[which.min(CVErr)]
    
    message(c("Param-Done: (j, p, N, noise) = ( ", as.character(j), 
              ", ", as.character(p),
              ", ", as.character(N),", ", as.character(noise), ")"))
    
    # Calculating AUCs
    AUCs <- foreach(k = 1:100, .combine = 'c', .packages = c("DescTools")) %dopar% {
        set.seed(1238487 + k^2)
        X <- MASS::mvrnorm(n  = N,
                           mu = rep(0, R),
                           Sigma = CovsReg[[j]] + noise*diag(R))
        Cov.Estimate  <- CovRegEstmtr(X)
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
Resolution  <- data.frame(p = c(20, 30, 40))
NoofSamples <- data.frame(N = c(50, 100, 200))
NoiseLevel  <- data.frame(noise = c(0, 0.01, 0.1))
Specs       <- tidyr::crossing(CovID, Resolution, NoofSamples, NoiseLevel)
######################################################################
######################################################################
AUCData <- Specs %>% 
    dplyr::mutate(AUCs = purrr::pmap(list(x = j, y = p, z = N, w = noise), function(x, y, z, w) {AUCs(x, y, z, w, CovsReg, OmegaTrue)}))
save(AUCData, file="AUCDataRegular.RData")
######################################################################
AUCSummary <- AUCData %>% 
    dplyr::mutate(Median = sapply(AUCs, median)) %>%
    dplyr::mutate(Median = round(Median, 2)) %>%
    dplyr::mutate(MAD = sapply(AUCs, mad)) %>%
    dplyr::mutate(MAD = round(MAD, 2)) %>%
    dplyr::mutate(Summary = paste(Median, MAD, sep = "\u00B1")) %>%
    subset(select = -c(AUCs, Median, MAD)) %>%
    tidyr::pivot_wider(names_from = N, values_from = Summary)
save(AUCSummary, file="AUCSummaryRegular.RData")
######################################################################
######################################################################
