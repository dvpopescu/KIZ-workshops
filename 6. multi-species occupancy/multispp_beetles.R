# NOTE: run with the working directory set to this script's folder
# (RStudio: Session > Set Working Directory > To Source File Location);
# the two CSV files are read from this folder.

library(spOccupancy)
library(coda)
library(tidyverse)
library(stars)
library(ggplot2)


# READ AND PREPARE DATA FOR spOccupancy -----------------------------------

# read covariates
vars <- read.csv("covariates_multi spp.csv", header=T)
str(vars)

#observation/detection covariate
Julian_Day <- vars %>%
  select(starts_with("Julian"))

names(Julian_Day) <- NULL
Julian_Day <- as.matrix(Julian_Day) # we need the data as a matrix

Precipitation <- vars %>%
  select(starts_with("Precip"))

names(Precipitation) <- NULL
Precipitation <- as.matrix(Precipitation) # we need the data as a matrix

# create a LIST of detection covariates
det.covs <- list(Julian = Julian_Day, Precip = Precipitation)
str(det.covs)

# create a LIST of occupancy (site-level) covariates
# (covariate rows are in the same trap order E1-E25, M1-M25 as the species data)
occ.covs <- list(Canopy = vars[,2], Groundcover = vars[,4], VolDeadWood = vars[,9])

# read and organize species  data
spp_data <- read.csv("all_species_occupancy.csv", header=T)
str(spp_data)

# identify species names
spp.names <- unique(spp_data$Species)
spp.names

# identify site names
site.code <- unique(spp_data$Trap)
site.code

# calculate number of species
N <- length(spp.names)
# Maximum number of replicates at a site
K <- 12
# Number of sites
J <- length(site.code)


# manipulate species observations into an 3D ARRAY with dimensions:
# 1) number of species (28); 2) number of sites (50) and 3) number of visits (12)
# this operation may differ based on the current format of your data, but in the 
# end it has to have the same structure as the 'hbef2015' data

# keep only columns with species name and detections and covert all abundance 
# numbers to 1's
df <- spp_data %>%
  select(!(Site:SpeciesID)) %>%
  mutate_if(is.numeric, ~1 * (. > 0))
df

# The steps below convert the data to an ARRAY of species x sites x visits 
#trial <- simplify2array(by(df,df$Species,as.matrix))
#str(trial)
#trial1 <- aperm(trial, c(3,1,2))
#str(trial1)
#trial2 <- trial1[,,-1]
#str(trial2)

# a more direct way is to combine all the step above into a single line of code 
# to convert the data to an ARRAY of species x sites x visits 
df1 <- aperm(simplify2array(by(df,df$Species,as.matrix)), c(3,1,2))[,,-1]
str(df1)

# the 1's and 0's are now characters (as.matrix() on a data frame that still
# contains the character 'Species' column gives a character matrix), so we have
# to covert them back to numbers
df2 <- as.numeric(df1)
str(df2)

# when we convert the data to numbers, it creates a very long vector of 
# 1's and 0's, so we have to place the data back into an ARRAY with the
# corect dimensions (species x sites x visits)
y <- array(df2, dim = c(N, J, K))
str(y)

# give species and site names to their respective ARRAY structures 
dimnames(y)[[1]] <- spp.names
dimnames(y)[[2]] <- site.code

str(y) 

### PHEW, we did it!!!!!! 


# Create spOccupancy object -----------------------------------------------


# create spOccupancy object
beetle_data <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
str(beetle_data)



# Multispecies occupancy analysis -----------------------------------------

# formula for detection and occupancy 
occ.ms.formula <- ~ scale(VolDeadWood) + scale(Canopy) + scale(Groundcover)
det.ms.formula <- ~ scale(Julian) + scale(Precip) + I(scale(Julian)^2)

# initial values for our coefficients (from the spOccupancy vignette)
N <- dim(beetle_data$y)[1]
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(beetle_data$y, c(1, 2), max, na.rm = TRUE))

# formulate prior distributions for coefficients; there se are done at a community-level,
# where all individual coefficients share the same distributions, thus species 
# with little data "borrow" information from species with good data 
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# RUN THE MODEL!
out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = beetle_data, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 4, 
                  verbose = TRUE, 
                  n.report = 6000, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 4)

# inspect summary
summary(out.ms)
str(out.ms)

# Examine trace plots
plot(out.ms, 'beta', density = TRUE) # Occupancy parameters.
plot(out.ms, 'alpha', density = TRUE) # Detection parameters.

# posterior predictive check
ppc.ms.out <- ppcOcc(out.ms, 'chi-squared', group = 1)
summary(ppc.ms.out)

# For multi-species models, fit.y and fit.y.rep are matrices
# (posterior samples x species), so we plot the check for one species at a time
# (here species 10; see spp.names for the order)
sp.ppc <- 10
ppc.df <- data.frame(fit = ppc.ms.out$fit.y[, sp.ppc],
                     fit.rep = ppc.ms.out$fit.y.rep[, sp.ppc],
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')


# Examine trace plots
plot(out.ms, 'beta', density = TRUE) # Occupancy parameters.
plot(out.ms, 'alpha', density = TRUE) # Detection parameters.

waicOcc(out.ms, by.sp = TRUE)



# Predictions -------------------------------------------------------------


deadwood.pred <- seq(min(scale(beetle_data$occ.covs$VolDeadWood)), 
                  max(scale(beetle_data$occ.covs$VolDeadWood)), 0.1)
# other covariates held at their (scaled) mean; length matches deadwood.pred
X.1 <- cbind(1, deadwood.pred,
             rep(mean(scale(beetle_data$occ.covs$Canopy)), length(deadwood.pred)),
             rep(mean(scale(beetle_data$occ.covs$Groundcover)), length(deadwood.pred)))
out.pred1 <- predict(out.ms, X.1)

str(out.pred1)

plot.dat1 <- data.frame(deadwood.pred, 
                        mean.psi = apply(out.pred1$psi.0.samples[,10,], 2, mean), 
                        upper.psi = apply(out.pred1$psi.0.samples[,10,], 2, 
                                          function(x) quantile(x, probs = 0.975)),
                        lower.psi = apply(out.pred1$psi.0.samples[,10,], 2, 
                                          function(x) quantile(x, probs = 0.025)),
                        stringsAsFactors = FALSE)

ggplot(data=plot.dat1, aes(x=deadwood.pred, y=mean.psi)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=plot.dat1$lower.psi, ymax=plot.dat1$upper.psi), linetype=2, alpha=0.1) +
  labs(title=paste(spp.names[10], "occupancy probability by VOlume of Dead Wood"),
       x ="Dead Wood (scaled)", y = "Occupancy probability")




