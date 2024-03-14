#The simplest Binomial N-Mixture Models

# Choose sample sizes and prepare observed data array C
set.seed(24) #So we all get same data set
M <- 100 # Number of sites
J <- 3 # Number of abu. measurements per site (rep.counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
lambda <- 4 # Expected abundance
p <- 0.6 # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}
# Look at data
# The truth ....
table(N) # True abundance distribution
N

sum(N) # True total population size at M sites
sum(N>0) # True number of occupied sites
mean(N) # True mean abundance (estimate of lambda)

table(apply(C, 1, max)) # Observed abundance distribution (max count)

sum(apply(C, 1, max)) [1] # Observed total population size at M sites # Observed number of occupied sites
sum(apply(C, 1, max)>0) [1] # Observed mean "relative abundance"
mean(apply(C, 1, max)) [1]  # Observed mean "relative abundance"

head(cbind(N=N, count1=C[,1], count2=C[,2])) # First 6 sites
cor(C)[1,2]
#[1] 0.4023582 Note excellent agreement with value of p

library(unmarked) # Load package
umf <- unmarkedFramePCount(y = C) # Create um data frame
summary(umf) # Summarize
(fm1 <- pcount(~1 ~1, data = umf)) # Fit model: get estimates on link scale 
# Get estimates on natural scale
backTransform(fm1, "state")
backTransform(fm1, "det")

#Binomial N-Mixture Models with covariate

set.seed(24) #So we all get same data set
M <- 150 # Number of sites
J <- 3 # Number of abu. measurements per site (rep.counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for abundance model and compute lambda
beta0 <- 0 # Log-scale intercept
beta1 <- 2 # Log-scale slope for vegHt
lambda <- exp(beta0 + beta1 * vegHt) # Expected abundance
plot(vegHt, lambda, type = "l", lwd = 3) # Expected abundance

# Draw local abundance and look at data so far
N <- rpois(M, lambda)
points(vegHt, N) # Add realized abundance to plot
table(N)

# Plot the true system state
#par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, N, xlab="Vegetation height", ylab="True abundance (N)", frame = F, cex = 1.5)
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")

# Create a co variate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2
alpha1 <- -3
p <- plogis(alpha0 + alpha1 * wind) # plot(p ~ wind, ylim = c(0,1))
# Logit-scale intercept
# Logit-scale slope for wind # Detection probability
# Look at relationship
# Take J [ 3 abundance measurements at each site
for(j in 1:J) {
C[,j] <- rbinom(M, N, p[,j]) }

# Plot observed data and effect of wind on det. probability
plot(wind, C/max(C), xlab="Wind", ylab="Scaled counts: C/max(C)", frame = F, cex = 1.5)
lines(seq(-1,1,,100), plogis(alpha0 + alpha1*seq(-1,1,,100)), lwd=3, col="red")

# Expected (lambda) and realized abundance (N) and measurements (C)
cbind(lambda=round(lambda,2), N=N, C1=C[,1], C2=C[,2], C3=C[,3])

time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 50), rep("B", 50), rep("C", 50)) # assumes M = 150

library(unmarked)
umf <- unmarkedFramePCount(
  y = C,# Counts matrix
  siteCovs = data.frame(vegHt = vegHt, hab = hab),  # Site covariates
  obsCovs = list(time = time, wind = wind)) # Observation covs

summary(umf)

# Fit model and extract estimates
# linear model for p follows first tilde, then comes linear model for lambda
summary(fm.Nmix1 <- pcount(~wind ~vegHt, data=umf, control=list(trace=T, REPORT=1)))

fm.Nmix2 <- pcount(~wind ~vegHt, data=umf, mixture="NB", control=list(trace=TRUE, REPORT=5))
fm.Nmix3 <- pcount(~wind ~vegHt, data=umf, mixture="ZIP", control=list(trace=TRUE, REPORT=5))
cbind(AIC.P=fm.Nmix1@AIC, AIC.NB=fm.Nmix2@AIC, AIC.ZIP=fm.Nmix3@AIC)

# Predictions of lambda for specified values of vegHt, say 0.2 and 2.1
#newdat <- data.frame(vegHt=c(0.2, 1.1))
#predict(fm.Nmix1, type="state", newdata=newdat, append = T)

# ... or of p for values of wind of -1 to 1
#newdat <- data.frame(wind=seq(-1, 1, , 5))
#predict(fm.Nmix1, type="det", newdata=newdat, append = T)

# Predict lambda and detection for actual data set
(lambda.hat <- predict(fm.Nmix1, type="state")) # lambda at every site
(p.hat <- predict(fm.Nmix1, type="det")) # p during every survey

# Fit detection-naive GLM to counts and plot comparison
summary(fm.glm <- glm(c(C) ~ rep(vegHt, 3), family=poisson)) # p-naive model
matplot(vegHt, C, xlab="Vegetation height", ylab="Counts", frame = F, cex = 1.5, pch = 1, col = "black")
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")
curve(exp(coef(fm.glm)[1]+coef(fm.glm)[2]*x), -1, 1, type ="l", lwd=3, add=TRUE)
lines(vegHt, predict(fm.Nmix1, type="state")[,1], col = "blue", lwd = 3)
legend(-1, 7, c("Truth", "'Poisson GLM' with p", "Poisson GLM without p"), col=c("red", "blue", "black"), lty = 1, lwd=3, cex = 1.2)

################################################################################
# The SIMPLEST POSSIBLE SITE-OCCUPANCY MODEL
# Choose sample sizes and prepare observed data array y

set.seed(24) #So we all get same data set
M <- 100 # Number of sites
J <- 2 # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8 # Probability of occupancy or presence 
p <- 0.5 # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}
# Look at data
sum(z) # True number of occupied sites
sum(apply(y, 1, max)) # Observed number of occupied sites

head(cbind(z = z, y1 = y[,1], y2 = y[,2])) # Truth and measurements for first 6 sites

library(unmarked)
umf <- unmarkedFrameOccu(y = y) # Create unmarked data frame 
summary(umf) # Summarize data frame
(fm1 <- occu(~1 ~1, data = umf)) # Fit model

backTransform(fm1, "state") # Get estimates on probability scale
backTransform(fm1, "det")

#SITE-OCCUPANCY MODEL WHIT COVARIATE

# Choose sample sizes and prepare obs. data array y
set.seed(1) # So we all get same data set
M <- 100 # Number of sites
J <- 3 # Number of presence/absence measurements 
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # Sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0 # Logit-scale intercept
beta1 <- 3 # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability

# plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship

# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi) # True presence/absence

# Look at data so far
table(z) 
z

# Plot the true system state
#par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2 # Logit-scale intercept
alpha1 <- -3 # Logit-scale slope for wind 
p <- plogis(alpha0 + alpha1 * wind) # Detection probability 
# plot(p ~ wind, ylim = c(0,1)) # Look at relationship

# Take J [ 3 presence/absence measurements at each site 
for(j in 1:J) {
y[,j] <- rbinom(M, z, p[,j]) }
sum(apply(y, 1, max)) # Number of sites with observed presences 

# Plot observed data and true effect of wind on detection probability
plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")

# Look at the data: occupancy, true presence/absence (z), and measurements (y)
cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34)) # Must have M = 100

# Load unmarked, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(
  y = y,                                            # Pres/Abs measurements
  siteCovs = data.frame(vegHt = vegHt, hab = hab),  # site-specific covs.
  obsCovs = list(wind = wind, time = time))         # obs-specific covs.
summary(umf) 

# Fit model and extract estimates
# Detection covariates follow first tilde, then occupancy covariates 
summary(fm1.occ <- occu(~wind ~vegHt, data=umf))

# Predict occupancy and detection as function of covs (with 95% CIs)
# Add truth from data simulation  
newdat <- data.frame(vegHt=seq(-1, 1, 0.01))
pred.occ <- predict(fm1.occ, type="state", newdata=newdat)
newdat <- data.frame(wind=seq(-1, 1, 0.1))
pred.det <- predict(fm1.occ, type="det", newdata=newdat)

# ... for values of wind of -1 to 1
newdat <- data.frame(wind=seq(-1, 1, , 5)) 
predict(fm1.occ, type="det", newdata=newdat, append = T)

# Fit detection-naive GLM to observed occurrence and plot comparison
summary(fm.glm <- glm(apply(y, 1, max) ~ vegHt, family=binomial))
plot(vegHt, apply(y, 1, max), xlab="Vegetation height", ylab="Observed occurrence ('ever observed ?')", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")
lines(vegHt, predict(fm.glm,,"response"), type = "l", lwd = 3)
lines(vegHt, predict(fm1.occ, type="state")[,1], col = "blue", lwd = 3)
legend(-1, 0.9, c("Truth", "'LR' with p", "LR without p"), col=c("red", "blue", "black"), lty = 1, lwd=3, cex = 1.2)
