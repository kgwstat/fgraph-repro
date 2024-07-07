######################################################################
library(tidyverse)
library(viridisLite)
library(RColorBrewer)
library(gridExtra)
source("./precn.r")
######################################################################
# ANALYSIS
######################################################################
# Pfizer 1 minute closing price tick data
data <- read.csv("./PFIZER__EQ__NSE__NSE__MINUTE.csv")
data$open   <- NULL
data$high   <- NULL
data$low    <- NULL
data$volume <- NULL
glimpse(data)

data <- data |> mutate(date = substr(timestamp, 1, 10),
                time = 100 * as.numeric(substr(timestamp, 12, 13)) + as.numeric(substr(timestamp, 15, 16))) # nolint

data$timestamp <- NULL
data <- data |> pivot_wider(names_from = date, values_from = close)
data$time <- NULL
data <- data |> as.matrix()
colnames(data) <- NULL

# Starting and closing times
data <- data[1:375, ]

# Applying log transformation
data <- data |> log()

# Subtracting log of initial price for every day.
data <- data |> apply(2, \(x) x - x[1])

# Computing Covariance and normalizing with maximum.
fun <- \(i, j) cov(data[i, ], data[j, ], use = "pairwise.complete.obs")
Cov <- outer(X = 1:375, Y = 1:375, FUN = Vectorize(fun))
Cov <-  Cov / max(Cov)

# Projecting to the cone of PSD matrices.
CovPos <- {
  egn <- eigen(Cov)
  D <- diag(egn$values |> (\(y) ifelse(y > 0, y, 0))())
  V <- egn$vectors
  V %*% D %*% t(V)
}

# Calculating norms of entries of the precision operator matrix P:
numbers::divisors(375)
p <- 25
P <- CovPos |> Prec(p) |> MatofOpNorms(p)

######################################################################
# PLOT
######################################################################
# Plotting the norms of entries of P
P |> log10() |> hist(breaks = 35, prob=TRUE, xlim = c(6.5, 10), col = "#cb4154", xlab = "Logarithm of Norm (base 10)", ylab = "Frequency", main = NULL)
P |> log10() |> density() |> lines()
abline(v = 8.7, col= "#00b300", lwd=3)
######################################################################

time <- c(
  ((900 + 15:59) |> (\(x) paste("0", substr(x, 1, 1), ":",substr(x, 2, 3), sep = ""))()),
  c(1000 + 0:59, 1100 + 0:59, 1200 + 0:59, 1300 + 0:59, 1400 + 0:59, 1500 + 0:29) |> 
  (\(x) paste(substr(x, 1, 2), ":", substr(x, 3, 4), sep = ""))()) |> as.matrix()

######################################################################
# PLOTS: Figure 5 (b)
######################################################################

# Cumulative Log-returns for Pfizer Limited.
set.seed(1234)
matplot(1:375, data[, sample(1:938, 75)], 
                              type = "l", 
                              xaxt = "n",
                              xlab = "Time (IST)", 
                              ylab = "Cumulative Log-return",
                              main = "Cumulative Log-returns for Pfizer Limited",
                               cex = 6)
axis(1, at = seq(1, 375, by = 15),
    labels = time[seq(1, 375, by = 15), 1],
       las = 2)

######################################################################
# PLOTS: Figure 5 (d)
######################################################################

x <- 1:25
y <- 1:25
grid <- expand.grid(x = x, y = y)
grid$z <- P %>% as.vector()
lattice::levelplot(z ~ x * y, data = grid,
                              main = "",
                              xlab = "Time (IST)", 
                              ylab = "Time (IST)",
                              scales = list(x = list(at = c(1, 7, 13, 19, 25),
                                                 labels = c("09:15", "10:45", "12:15", "13:45", "15:15"), 
                                                    rot = 90),
                                            y = list(at = c(1, 7, 13, 19, 25),
                                                 labels = c("09:15", "10:45", "12:15", "13:45", "15:15"))),
                                     colorkey = TRUE,
                                  col.regions = viridis(100), 
                                       aspect = "full")

######################################################################
# PLOTS: Figure 5 (f)
######################################################################

# Plotting the CI graph of the data.
grid$z <- (P > 10^{8.7}) %>% as.vector()
lattice::levelplot(z ~ x * y, data = grid,
                              main = "",
                              xlab = "Time (IST)", 
                              ylab = "Time (IST)",
                              scales = list(x = list(at = c(1, 7, 13, 19, 25),
                                                 labels = c("09:15", "10:45", "12:15", "13:45", "15:15"), 
                                                    rot = 90),
                                            y = list(at = c(1, 7, 13, 19, 25),
                                                 labels = c("09:15", "10:45", "12:15", "13:45", "15:15"))),
                                     colorkey = list(at = c(0, 0.5, 1),
                                                 labels = c("0", "", "1"),
                                                    col = viridis(3)),
                                  col.regions = viridis(2), 
                                       aspect = "full")

######################################################################


