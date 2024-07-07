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
######################################################################
# R = 600
# p = 10, 20, 30, 40
# N = 500
# noise = 0.01, 0.1
######################################################################
cl <- makeCluster(27)
registerDoParallel(cl)

load("./CovDetail.RData")
AUCs <- function(j, p, N, noise, Covs, OmegaTrue) {
    R <- 600
    AUCs <- foreach(k = 1:100, .combine = 'c', .packages = c("DescTools")) %dopar% {
        source("./metrics.r")
        source("./precn.r")
        source("./CovEstmtr.r")
        set.seed(1238487 + k^2)
        X <- MASS::mvrnorm(n = N, mu = rep(0, R), Sigma = Covs[[j]] + noise*(sum(diag(Covs[[j]]))/R)*diag(R))
        r <- 5
        for(i in 1:N) {
            s <- sample(1:R, r)
            X[i, -s] <- 0
        } # sparsifying data
        Cov.Estimate <- CovSparseEstmtr(X, 10, r) # estimating covariance
        Prec.Estimate <- PrecReg(Cov.Estimate, p)
        Omega.True <- OmegaTrue[[j]][[p/10]]
        AreaUnderCurveTrapezoid(Prec.Estimate, Omega.True)
    }
    message(c("Done: (j, p, N, noise) = ( ", as.character(j), 
              ", ", as.character(p),
              ", ", as.character(N),", ", as.character(noise), ")"))
    return(AUCs)
}
######################################################################
######################################################################
CovID       <- data.frame(j = c(1, 2, 3, 4, 5))
Resolution  <- data.frame(p = c(10, 20, 30 , 40))
NoofSamples <- data.frame(N = c(500))
NoiseLevel  <- data.frame(noise = c(0, 0.01, 0.1))
Specs <- tidyr::crossing(CovID, Resolution, NoofSamples, NoiseLevel)
######################################################################
AUCData <- Specs %>% 
    dplyr::mutate(AUCs = purrr::pmap(list(x = j, y = p, z = N, w = noise), function(x, y, z, w) {AUCs(x, y, z, w, Covs, OmegaTrue)}))
save(AUCData, file="AUCData.RData")

stopCluster(cl)

AUCSummary <- AUCData %>% 
    dplyr::mutate(Median = sapply(AUCs, median)) %>%
    dplyr::mutate(Median = round(Median, 2)) %>%
    dplyr::mutate(MAD = sapply(AUCs, mad)) %>%
    dplyr::mutate(MAD = round(MAD, 2)) %>%
    dplyr::mutate(Summary = paste(Median, MAD, sep = "\u00B1")) %>%
    subset(select = -c(AUCs, Median, MAD)) %>%
    tidyr::pivot_wider(names_from = N, values_from = Summary)
save(AUCSummary, file="AUCSummary.RData")
######################################################################
######################################################################