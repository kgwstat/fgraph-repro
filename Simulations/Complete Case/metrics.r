######################################################################
TPR <- function(A, B) {
    if (sum(B) != 0) {sum(A & B) / sum(B)} else {0}
} # True Positive Rate
FPR <- function(A, B) {
    if (sum(!B) != 0) {sum(A & !B) / sum(!B)} else {0}
} # False Positive Rate
######################################################################
AreaUnderCurve <- function(Prec.Estimate, Omega.True) {
    p <- nrow(Omega.True)
    TPR2 <- function(A) { TPR(A, Omega.True) }
    FPR2 <- function(A) { FPR(A, Omega.True) }
    # Curried versions of TPR and FPR
    MatOfOps <- MatofOpNorms(Prec.Estimate, p)
    OmegaEstimate <- function(x) { MatOfOps > x } 
    # Estimate of Omega for given threshold x
    max <- MatOfOps %>% as.vector() %>% max()
    # Maximum of the norms of entries of Prec.Estimate
    Thresholds <- seq(0, 1, length = 1000) * max
    FPRs <- Thresholds %>% lapply(OmegaEstimate) %>% lapply(FPR2) %>% unlist()
    TPRs <- Thresholds %>% lapply(OmegaEstimate) %>% lapply(TPR2) %>% unlist()
    AUC  <- DescTools::AUC(FPRs, TPRs, method = "spline")
    # Calculates TPR and FPR for different thresholds
    return(AUC)
} # Calculates area under the ROC curve obtained
######################################################################
AreaUnderCurveTrapezoid <- function(Prec.Estimate, Omega.True) {
    p <- nrow(Omega.True)
    TPR2 <- function(A) { TPR(A, Omega.True) }
    FPR2 <- function(A) { FPR(A, Omega.True) }
    # Curried versions of TPR and FPR
    MatOfOps <- MatofOpNorms(Prec.Estimate, p)
    OmegaEstimate <- function(x) { MatOfOps > x } 
    # Estimate of Omega for given threshold x
    max <- MatOfOps %>% as.vector() %>% max()
    # Maximum of the norms of entries of Prec.Estimate
    Thresholds <- seq(0, 1, length = 1000) * max
    FPRs <- Thresholds %>% lapply(OmegaEstimate) %>% lapply(FPR2) %>% unlist()
    TPRs <- Thresholds %>% lapply(OmegaEstimate) %>% lapply(TPR2) %>% unlist()
    AUC  <- DescTools::AUC(FPRs, TPRs, method = "trapezoid")
    # Calculates TPR and FPR for different thresholds
    return(AUC)
} # Calculates area under the ROC curve obtained using trapezoidal rule.
######################################################################