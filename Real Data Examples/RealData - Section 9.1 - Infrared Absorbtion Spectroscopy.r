######################################################################
library(tidyverse)
library(viridisLite)
library(RColorBrewer)
library(gridExtra)
source("./precn.r")
######################################################################
# ANALYSIS
######################################################################
# importing data
data <- read.csv("./MIR_Fruit_purees.csv")
glimpse(data)

# wavelength
wavelength <- data[, 1]
r <- length(wavelength)

# spectra
spectra <- {
              y <- colnames(data)
              y <- y[startsWith(y, "Strawberry")]
              data[, y]
           }
colnames(spectra) <- NULL
N <- ncol(spectra)

# normalizing spectra
spectra <- apply(spectra, 2, \(x) (x * 235 / sum(x)))

######################################################################
# PLOTS: Figure 5 (a)
######################################################################
matplot(wavelength, spectra, type = "l",
                             xlab = "Wavelength (nm)", 
                             ylab = "Normalized Intensity",
                             main = "Normalized Absorption Spectra from Strawberry Purees", 
                              cex = 6)
######################################################################
# ANALYSIS
######################################################################

# calculating mean and covariance
mean    <- rowMeans(spectra)
spectra <- apply(spectra, 2, \(x) (x - mean))
Cov     <- spectra %*% t(spectra) / N

# we are going to bring the grid resolution from 235 to 234 so as to ensure easy division by p = 39.
Cov <- Cov[-235, -235]
numbers::divisors(234)

# calculating the precision
p <- 39
P <- Cov |> Prec(p) |> MatofOpNorms(p)
P |> log10() |> hist(breaks = 200)
(P > 10^{7.6}) |> lattice::levelplot()
P |> lattice::levelplot()

######################################################################
# PLOTS: Figure 5 (c)
######################################################################

x <- seq(900, 1800, length = p)
y <- seq(900, 1800, length = p)
grid <- expand.grid(x = x, y = y)
grid$z <- P %>% as.vector()
lattice::levelplot(z ~ x * y, 
                data = grid,
                main = "",
                xlab = "Wavelength (nm)", 
                ylab = "Wavelength (nm)",
            colorkey = TRUE,
         col.regions = viridis(100), 
              aspect = "full")

######################################################################
# PLOTS: Figure 5 (e)
######################################################################

grid$z <- (P > 10^{7.4}) %>% as.vector()
lattice::levelplot(z ~ x * y, 
                data = grid,
                main = "",
                xlab = "Wavelength (nm)", 
                ylab = "Wavelength (nm)",
            colorkey = list(at = c(0, 0.5, 1),
                        labels = c("0", "", "1"),
                           col = viridis(3)),
         col.regions = viridis(2), 
              aspect = "full")
              
######################################################################

