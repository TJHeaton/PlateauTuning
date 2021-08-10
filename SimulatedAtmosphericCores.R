# We create simulated, hypothetical, atmospheric cores
# Select cores with similar propoerties to Lake Suigetsu 14C record (in terms of density and precision) 
# Use IntCal20 curve as baseline and add observational noise ro create simulated record
# Try and identify plateaus from these simulated cores
# This code will recreate plots for Figures 5 and 6
# Choose bandsd = 50 for Figure 5 and bandsd = 100 for Figure 6


# Edits:
# 10th August 2021
# Improve comments to upload to Github

# 19th April 2021
# Post review add gradient lines at 0, 0.5 and 1 

# Choose seed so that reproducible 
set.seed(14)
par(las = 1)

# Select parameters for simulation study 
range <- c(12000, 13900)
kneigh <- 60 # Number of neighbouring observations to use in local gradient estimation

# Select bandwidth for use in estimating the local gradient - larger means greater weights to further points
bandsd <- 50 # 50 for Figure 5 and 100 for Figure 6


source("PlateauMatchingFunctions.R")

# Install colorspace for plotting 
library(colorspace)
pal <- terrain.colors(4)[-4]
pal[3] <- "black" # Choose black as one of the palette

# Load IntCal20 curve as atmospheric baseline 
calcurve <- read.csv("IntCal20TreeRing.csv", header = TRUE)

# Read in Suigetsu 14C data and extract number of observations in range and c14sd (so simulated cores can mimic it)
Suigetsu <- read.csv("Suigetsu2013.csv", header = TRUE, sep = " ")
Suiginrange <- which(Suigetsu$calage < range[2] & Suigetsu$calage > range[1])
Suign <- length(Suiginrange)
Suigobssd <- Suigetsu$c14sig[Suiginrange]
rm(Suiginrange)

# Create simulated pseudo-Suigetsu atmospheric cores based upon even sampled calendar ages
TestCoreA <- CreateObsDataEvenSpacing(Suign, calcurve = calcurve, range = range, obssd = Suigobssd)
# Estimate the local gradient 
TestGradA <- sapply(range[1]:range[2], function(x) 
  FindGrad(x, CoreInfo = TestCoreA, k = kneigh, bandsd = bandsd))


# Create another two pseudo-Suigetsu atmospheric cores
TestCoreB <- CreateObsData(Suign, calcurve = calcurve, range = range, obssd = Suigobssd)
TestGradB <- sapply(range[1]:range[2], function(x) 
  FindGrad(x, CoreInfo = TestCoreB, k = kneigh, bandsd = bandsd))

# Create another pseudo core 
TestCoreC <- CreateObsData(Suign, calcurve = calcurve, range = range, obssd = Suigobssd)
TestGradC <- sapply(range[1]:range[2], function(x) 
  FindGrad(x, CoreInfo = TestCoreC, k = kneigh, bandsd = bandsd))


####################################################################################################
# Create plot showing each core with gradient estimate along the bottom
# For figures in paper we selected plotting window with width = 12, height = 8) 

# Change plotting layout to allow 4 graphs on single plot 
par(mfrow = c(2,2), mar = c(5, 6, 4, 4.5) + 0.1, las = 1, mgp = c(3.5, 1, 0))
# CHoose suitable range for y-axis for plots
yplotrange <- approx(x = calcurve[,1], y = calcurve[,2], xout = range)$y + (3 * max(Suigobssd) * c(-1, 1)) + c(-700, 0)


# Plot simulated core A
plot(calcurve[,1], calcurve[,2], xlim = range, ylim = yplotrange, type = "l",
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")))
points(TestCoreA$theta, TestCoreA$c14age, pch = 19, cex = 0.6, col = "blue")
par(new = T, las = 0)
plot(TestGradA[1,], TestGradA[2,], col = "blue", type = "l", axes = FALSE, xlab = NA, ylab = NA, ylim = c(-1,3))
polygon(c(rev(TestGradA[1,]), TestGradA[1,]), 
        c(rev(TestGradA[2,] + 2 * TestGradA[3,]),  TestGradA[2,] - 2 * TestGradA[3,]), 
        col=rgb(0,0,1,.5), border=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Estimated gradient')
mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
      line = -1.3)
abline(h = c(0, 0.5, 1), col = pal, lty = 2, lwd = 2)
text(x = 13500, y = 0.21, "Plateau Gradient Thresholds", col = "black")

# Now plot the plateaus according to the threshold
par(new = TRUE, las = 0)
plot(TestGradA[1,], TestGradA[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThresh(TestGradA, thresh = 0, col = pal[1], lwd = 3)
PlotBelowThresh(TestGradA, thresh = 0.5, col = pal[2], lwd = 3)
PlotBelowThresh(TestGradA, thresh = 1, col = pal[3], lwd = 3)
text(x = 13250, y = 3.5, expression("Identified Periods Below Gradient Thresholds"), col = "black")


# Plot simulated core B
par(las = 1)
plot(calcurve[,1], calcurve[,2], xlim = range, ylim = yplotrange, type = "l",
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")))
points(TestCoreB$theta, TestCoreB$c14age, pch = 19, cex = 0.6, col = "red")
par(new = T, las = 0)
plot(TestGradB[1,], TestGradB[2,], col = "red", type = "l", axes = FALSE, xlab = NA, ylab = NA, ylim = c(-1,3))
polygon(c(rev(TestGradB[1,]), TestGradB[1,]), 
        c(rev(TestGradB[2,] + 2 * TestGradB[3,]),  TestGradB[2,] - 2 * TestGradB[3,]), 
        col=rgb(1,0,0,.5), border=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Estimated gradient')
mtext(paste0("(", letters[2], ")"), side = 3, adj = 0.05, 
      line = -1.3)
abline(h = c(0, 0.5, 1), col = pal, lty = 2, lwd = 2)
text(x = 13500, y = 0.21, "Plateau Gradient Thresholds", col = "black")
# Now plot the plateaus according to the threshold
par(new = TRUE, las = 0)
plot(TestGradB[1,], TestGradB[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThresh(TestGradB, thresh = 0, col = pal[1], lwd = 3)
PlotBelowThresh(TestGradB, thresh = 0.5, col = pal[2], lwd = 3)
PlotBelowThresh(TestGradB, thresh = 1, col = pal[3], lwd = 3)
text(x = 13250, y = 3.5, expression("Identified Periods Below Gradient Thresholds"), col = "black")

# Plot simulated core C
par(las = 1)
plot(calcurve[,1], calcurve[,2], xlim = range, ylim = yplotrange, type = "l",
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")))
points(TestCoreC$theta, TestCoreC$c14age, pch = 19, cex = 0.6, col = "green")
par(new = T, las = 0)
plot(TestGradC[1,], TestGradC[2,], col = "green", type = "l", axes = FALSE, xlab = NA, ylab = NA, ylim = c(-1,3))
polygon(c(rev(TestGradC[1,]), TestGradC[1,]), 
        c(rev(TestGradC[2,] + 2 * TestGradC[3,]),  TestGradC[2,] - 2 * TestGradC[3,]), 
        col=rgb(0,1,0,.5), border=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Estimated gradient')
mtext(paste0("(", letters[3], ")"), side = 3, adj = 0.05, 
      line = -1.3)
abline(h = c(0, 0.5, 1), col = pal, lty = 2, lwd = 2)
text(x = 13500, y = 0.21, "Plateau Gradient Thresholds", col = "black")
# Now plot the plateaus according to the threshold
par(new = TRUE, las = 0)
plot(TestGradC[1,], TestGradC[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThresh(TestGradC, thresh = 0, col = pal[1], lwd = 3)
PlotBelowThresh(TestGradC, thresh = 0.5, col = pal[2], lwd = 3)
PlotBelowThresh(TestGradC, thresh = 1, col = pal[3], lwd = 3)
text(x = 13250, y = 3.5, expression("Identified Periods Below Gradient Thresholds"), col = "black")

# Now plot the overlay of the estimated gradients 
par(las = 0)
plot(TestGradA[1,], TestGradA[2,], type = "l", col = "blue",
     xlab = "Calendar Age (cal yr BP)", ylab = "Linear Model Gradient", ylim = c(-1,3))
polygon(c(rev(TestGradA[1,]), TestGradA[1,]), 
        c(rev(TestGradA[2,] + 2 * TestGradA[3,]),  TestGradA[2,] - 2 * TestGradA[3,]), 
        col=rgb(0,0,1,.5), border=NA)
lines(TestGradB[1,], TestGradB[2,], col = "red")
polygon(c(rev(TestGradB[1,]), TestGradB[1,]), 
        c(rev(TestGradB[2,] + 2 * TestGradB[3,]),  TestGradB[2,] - 2 * TestGradB[3,]), 
        col=rgb(1,0,0,.5), border=NA)
lines(TestGradC[1,], TestGradC[2,], col = "green")
polygon(c(rev(TestGradC[1,]), TestGradC[1,]), 
        c(rev(TestGradC[2,] + 2 * TestGradC[3,]),  TestGradC[2,] - 2 * TestGradC[3,]), 
        col=rgb(0,1,0,.5), border=NA)
mtext(paste0("(", letters[4], ")"), side = 3, adj = 0.05, 
      line = -1.3)
abline(h = 1, col = "black", lty = 2, lwd = 2)
text(x = 13380, y = 0.71, "Sarnthein Plateau Gradient Threshold", col = "black")

# Now plot as a rug along the bottom the plateaus thresholds
par(new = TRUE, las = 0)
plot(TestGradC[1,], TestGradC[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThreshJitter(TestGradA, thresh = 1, col = "blue", lwd = 3, yval = 1)
PlotBelowThreshJitter(TestGradB, thresh = 1, col = "red", lwd = 3, yval = 0.5)
PlotBelowThreshJitter(TestGradC, thresh = 1, col = "green", lwd = 3, yval = 0)
text(x = 13400, y = 3.5, expression("Identified Periods With Gradient"<" 1"), col = "black")
