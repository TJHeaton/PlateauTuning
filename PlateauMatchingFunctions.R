# 10th August 2021
# Functions to create simulated cores and estimate local gradients using local linear model
###############################################


# Function to create hypothetical 14C observations based on calibration curve
# Will sample n calendar ages Unif[range] (and sorts them into ascending order)
# Then create 14C observations using calibration curve (and uncertainty) 
# Arguments:
# n - the number of observations to create
# calcurve - the calibration curve to use
# range - the calendar range over which to sample the hypothetical data 
# obssd - the chosen observational sd
CreateObsData <- function(n, calcurve, range, obssd) {
  # Replace obssd with vector if not of the correct length 
  if(length(obssd) != n) obssd <- rep(obssd[1], n)
  
  # Sample some calendar ages and order them     
  samptheta <- sort(sample(range[1]:range[2], n, replace = TRUE))
  sampmu <- approx(x = calcurve[,1], y = calcurve[,2], xout = samptheta)$y
  obsy <- rnorm(n, mean = sampmu, sd = obssd) 
  retlist <- list(theta = samptheta, c14age = obsy, c14sig = obssd)
  return(retlist)
}

# Function to create hypothetical 14C observations but with more evenly spaced calendar ages 
# As above but create more even calendar age spacing:
# Sampling 2n ~ Unif[range] then choosing every 2nd even value 
CreateObsDataEvenSpacing <- function(n, calcurve, range, obssd) {
  # Replace obssd with vector if not of the correct length 
  if(length(obssd) != n) obssd <- rep(obssd[1], n)
  
  # Sample some calendar ages and order them     
  samptheta <- sort(sample(range[1]:range[2], 2 * n, replace = TRUE))[2*(1:n)]
  sampmu <- approx(x = calcurve[,1], y = calcurve[,2], xout = samptheta)$y
  obsy <- rnorm(n, mean = sampmu, sd = obssd) 
  retlist <- list(theta = samptheta, c14age = obsy, c14sig = obssd)
  return(retlist)
}


# Function to estimate the local gradient of the observed values
# Fits a linear model to the nearest k observations
# Arguments:
# t - the calendar age at which we want to estimate the gradient
# CoreInfo - the information on the core, needs to contain: 
#                   theta - the calendar ages 
#                   c14sig - the measurement uncertainty
# k - the number of 14C observations to use in estimating the gradient
# bandsd - the bandwith for the kernel to create weights N(0, bandsd^2)
# Returns:
# grad - local gradient estimate based upon linear model
FindGrad <- function(t, CoreInfo, k = 15, bandsd = 100) {
  # Centre the calendar ages around the time of interest
  CoreInfo$shiftcal <- CoreInfo$theta - t
  useid <- order(abs(CoreInfo$shiftcal))[1:k]
  # Reduce the core info to the nearest k observations 
  subcore <- as.data.frame(CoreInfo)[useid,]
  weights <-  dnorm(subcore$shiftcal, mean = 0, sd = bandsd) / subcore$c14sig^2
  
  # Now fit a linear model to this data and extract the gradient (and se)
  lm1 <- lm(c14age ~ shiftcal, data = subcore, weights = weights)  
  retlist <- as.vector(c(t = t, grad = lm1$coefficients[2], se = sqrt(vcov(lm1)[2,2])))
  return(retlist)
}