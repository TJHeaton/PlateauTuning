### NOTE: Must run SimulatedMarineCores.R first (and then not remove anything)  
# This uses bandsd as defined in SimulatedMarineCores.R

# Create simulated, hypothetical, pseudo-marine cores which will be compared with simulated atmospheric cores
# Select cores with similar propoerties to Cariaco 14C record (in terms of density and precision) 
# Use IntCal20 curve as baseline and add observational noise to create simulated record
# Try and identify plateaus from these simulated cores
# This code will recreate plots for Figure 7

# When estimating the local gradient , we assume our simulated marine cores don't have any depth information
# Instead we just use a set of neighbours (assume observations equally spaced in depth along the core)

# Edits:
# 10th August 2021
# Improve comments to upload to Github

# 19th April 2021
# Post review add gradient lines at 0, 0.5 and 1 

# 29th January 2020
# This version uses just the uncertainty in the 14C determinations (does not use MRA uncertainty) 

set.seed(63)

# Read in Cariaco data and use this to create marine cores with similar density and precision 
Cariaco <- read.csv("Cariaco2013Raw.csv", header = TRUE, sep = "\t")
# Extract Cariaco informatio form similar range (slightly older as no Cariaco in chosen range) 
Cariacoinrange <- which(Cariaco$calage < (2000 + range[2]) & Cariaco$calage > (2000 + range[1]))
Cariacon <- length(Cariacoinrange) # Number of observations we will sample in our simulated cores
Cariacoobssd <- Cariaco$c14sig[Cariacoinrange]
rm(Cariacoinrange)

# Create two simulated pseudo-Cariaco marine cores based upon even sampled calendar ages
CariacoCoreA <- CreateObsDataEvenSpacing(Cariacon, calcurve = calcurve, range = range, obssd = Cariacoobssd)
CariacoCoreB <- CreateObsDataEvenSpacing(Cariacon, calcurve = calcurve, range = range, obssd = Cariacoobssd)

# Adjust the c14age by a constant MRA of 400 C14yrs
CariacoCoreA$c14age <- CariacoCoreA$c14age + 400 
CariacoCoreB$c14age <- CariacoCoreB$c14age + 400

# Shift pseduo-Cariaco onto a depth only timescale where only know ordering
# truetheta is the genuine underlying calendar age from which we have sampled
# theta is pseudo-age under assumption of even depth spacing of 14C observations in core 
CariacoCoreA$truetheta <- CariacoCoreA$theta
CariacoCoreA$theta <- seq(range[1], range[2], length = Cariacon)
CariacoCoreB$truetheta <- CariacoCoreB$theta
CariacoCoreB$theta <- seq(range[1], range[2], length = Cariacon)

# Find the implied sedimentation rate that would lead to even depth spacing
CariacoCoreA$sedrate <- c(1/(CariacoCoreA$truetheta[-1] - CariacoCoreA$truetheta[-Cariacon]), NA)   
CariacoCoreA$sedrate <- CariacoCoreA$sedrate / mean(CariacoCoreA$sedrate, na.rm = TRUE)

CariacoCoreB$sedrate <- c(1/(CariacoCoreB$truetheta[-1] - CariacoCoreB$truetheta[-Cariacon]), NA)   
CariacoCoreB$sedrate <- CariacoCoreB$sedrate / mean(CariacoCoreB$sedrate, na.rm = TRUE)

# Estimate the local pseudo-gradient on the even-depth spacing scale
CariacoGradA <- sapply(range[1]:range[2], function(x) 
  FindGrad(x, CoreInfo = CariacoCoreA, k = Cariacon, bandsd = bandsd))
CariacoGradB <- sapply(range[1]:range[2], function(x) 
  FindGrad(x, CoreInfo = CariacoCoreB, k = Cariacon, bandsd = bandsd))

# Convert pseudo-gradient back to the original true Cariaco timescale
temptA <- seq(min(CariacoCoreA$truetheta), max(CariacoCoreA$truetheta), length = 600) # Underlying true cal age grid 
pseudotA <- approx(x = CariacoCoreA$truetheta, y = CariacoCoreA$theta, xout = temptA)$y # Corresponding pseudo-ages depths
gradestA <- approx(x = CariacoGradA[1,], y = CariacoGradA[2,], xout = pseudotA)$y # Local gradient at these pseudo-depths

temptB <- seq(min(CariacoCoreB$truetheta), max(CariacoCoreB$truetheta), length = 600)
pseudotB <- approx(x = CariacoCoreB$truetheta, y = CariacoCoreB$theta, xout = temptB)$y
gradestB <- approx(x = CariacoGradB[1,], y = CariacoGradB[2,], xout = pseudotB)$y


####################################################################################################
# Create plot showing each core with gradient estimate along the bottom
# For figures in paper we selected plotting window with width = 12, height = 8) 

# Change plotting layout to allow 4 graphs on single plot 
par(mfrow = c(2,2), mar = c(5, 6, 4, 4.5) + 0.1, las = 1, mgp = c(3.5, 1, 0))
# Choose suitable range for y-axis for plots
yplotrange <- approx(x = calcurve[,1], y = calcurve[,2], xout = range)$y + (3 * max(Cariacoobssd) * c(-1, 1)) + c(0, 400)
ylimsedrate <- c(0, 10*max(c(CariacoCoreA$sedrate, CariacoCoreB$sedrate), na.rm = TRUE))


# Plot A: The true data (on true calendar age scale)
plot(calcurve[,1], calcurve[,2], xlim = range, ylim = yplotrange, type = "l",
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")))
points(CariacoCoreA$truetheta, CariacoCoreA$c14age, pch = 19, cex = 0.6, col = "purple")
points(CariacoCoreB$truetheta, CariacoCoreB$c14age, pch = 19, cex = 0.6, col = "orange")
legend("bottomright", legend = c("Pseudo marine core A", "Pseudo marine core B"), pch = 19, col = c("purple", "orange"))
mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
      line = -1.3)

# Plot B: Converted to equi-depth (i.e. where we do not know the calendar age scale)
plot(CariacoCoreA$theta, CariacoCoreA$c14age, pch = 19, cex = 0.6, col = "purple", axes = FALSE,
     xlab = "Ordering Only", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")), ylim = yplotrange)
points(CariacoCoreB$theta, CariacoCoreB$c14age, pch = 19, cex = 0.6, col = "orange")
axis(side = 2)
par(new = T, las = 0)
plot(CariacoGradA[1,], CariacoGradA[2,], type = "l", col = "purple", axes = FALSE, 
     xlab = NA, ylab = NA, ylim = c(-2,4))
lines(CariacoGradB[1,], CariacoGradB[2,], col = "orange")
mtext(side = 4, line = 3, 'Estimated pseudo-gradient')
box()
polygon(c(rev(CariacoGradA[1,]), CariacoGradA[1,]), 
        c(rev(CariacoGradA[2,] + 2 * CariacoGradA[3,]),  CariacoGradA[2,] - 2 * CariacoGradA[3,]), 
        col=rgb(1,0,1,.5), border=NA)
polygon(c(rev(CariacoGradB[1,]), CariacoGradB[1,]), 
        c(rev(CariacoGradB[2,] + 2 * CariacoGradB[3,]),  CariacoGradB[2,] - 2 * CariacoGradB[3,]), 
        col=rgb(1,0.647,0,.5), border=NA)
mtext(paste0("(", letters[2], ")"), side = 3, adj = 0.05, 
      line = -1.3)
axis(side = 4, at = -1:4)
abline(h = c(0, 0.5, 1), col = pal, lty = 2, lwd = 2)
text(x = 13500, y = 0.21 - 0.5, "Plateau Gradient Thresholds", col = "black")

# Plot the plateaus under the gradient thresholds
par(new = TRUE, las = 0)
plot(TestGradA[1,], CariacoGradA[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(-31,4))
PlotBelowThresh(CariacoGradA, thresh = 1, col = "purple", lwd = 3)
PlotBelowThreshJitter(CariacoGradB, thresh = 1, col = "orange", lwd = 3, yval = 0.5)
text(x = 13400, y = 3, expression("Identified Periods With Gradient"<" 1"), col = "black")

# Plot sedimentation rates
par(new = TRUE, las = 0)
plot(CariacoCoreA$theta, CariacoCoreA$sedrate, type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = ylimsedrate)
segments(x0 = CariacoCoreA$theta[-Cariacon], x1 = CariacoCoreA$theta[-1], y0 = CariacoCoreA$sedrate[-Cariacon], col = "purple", lwd = 2, lty = 1)
segments(x0 = CariacoCoreB$theta[-Cariacon], x1 = CariacoCoreB$theta[-1], y0 = CariacoCoreB$sedrate[-Cariacon], col = "orange", lwd = 2, lty = 1)
axis(4, at = seq(0, 3, length =4), labels = NA, col = "blue", padj = -1)
abline(h = 0, col = "black")
abline(h = 4, col = "black")
mtext(side = 4, line = 1, 'Rel. Sed.', col ="blue", adj = -0.05, cex = 0.7)


# Plot C: Plot transformed gradient on the original scale
par(las = 1)
plot(calcurve[,1], calcurve[,2], xlim = range, ylim = yplotrange, type = "l",
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")))
points(CariacoCoreA$truetheta, CariacoCoreA$c14age, pch = 19, cex = 0.6, col = "purple")
points(CariacoCoreB$truetheta, CariacoCoreB$c14age, pch = 19, cex = 0.6, col = "orange")
par(new = T, las = 0)
plot(temptA, gradestA, col = "purple", type = "l", axes = FALSE, xlab = NA, ylab = NA, ylim = c(-2,4))
lines(temptB, gradestB, col = "orange")
mtext(side = 4, line = 3, 'Estimated pseudo-gradient')
mtext(paste0("(", letters[3], ")"), side = 3, adj = 0.05, 
      line = -1.3)
axis(side = 4)
abline(h = c(0, 0.5, 1), col = pal, lty = 2, lwd = 2)
# Plot the plateaus under the gradient thresholds
par(new = TRUE, las = 0)
plot(TestGradC[1,], TestGradC[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThresh(rbind(temptA, gradestA), thresh = 1, col = "purple", lwd = 3)
PlotBelowThreshJitter(rbind(temptB, gradestB), thresh = 1, col = "orange", lwd = 3, yval = 0.5)
text(x = 13400, y = 3.5, expression("Identified Periods With Gradient"<" 1"), col = "black")

# Plot D: Overlay the chosen plateau vs plateaus identified in simulated pseudo-Suigetsu atmospheric record A
plot(TestGradA[1,], TestGradA[2,], type = "l", col = "blue",
     xlab = "Calendar Age (cal yr BP)", ylab = "Linear Model Gradient", ylim = c(-2,4))
polygon(c(rev(TestGradA[1,]), TestGradA[1,]), 
        c(rev(TestGradA[2,] + 2 * TestGradA[3,]),  TestGradA[2,] - 2 * TestGradA[3,]), 
        col=rgb(0,0,1,.5), border=NA)
lines(temptA, gradestA, col = "purple")
lines(temptB, gradestB, col = "orange")
abline(h = 1, col = "black", lty = 2, lwd = 2)
text(x = 13400, y = 0.71, "Sarnthein Plateau Gradient Threshold", col = "black")
legend("topright", legend = c("Pseudo atmosphere", "Pseudo marine core A", "Pseudo marine core B"), 
        lty = 1, col = c("blue", "purple", "orange"))
mtext(paste0("(", letters[4], ")"), side = 3, adj = 0.05, 
        line = -1.3)
  
# Now plot as a rug along the bottom the plateaus thresholds
par(new = TRUE, las = 0)
plot(TestGradC[1,], TestGradC[2,], type = "n", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0,35))
PlotBelowThreshJitter(TestGradA, thresh = 1, col = "blue", lwd = 3, yval = 1)
PlotBelowThreshJitter(rbind(temptA, gradestA), thresh = 1, col = "purple", lwd = 3, yval = 0.5)
PlotBelowThreshJitter(rbind(temptB, gradestB), thresh = 1, col = "orange", lwd = 3, yval = 0)  
text(x = 13400, y = 3.5, expression("Identified Periods With Gradient"<" 1"), col = "black")



