
### This code analyzes data on longicorn beetle Plagionotus detritus collected
### in Summer 2015 the Iron Gates Natural Park, Romania, and published in 
### Brodie, B. S., V. D. Popescu, R. Iosif, et al. 2019. Non-lethal monitoring of 
### longicorn beetle communities using generic pheromone lures and occupancy models. 
### Ecological Indicators 101:330-340.

setwd('...')

beetles <- read.csv("cerambycidae.csv")
names(beetles)
attach(beetles)
str(beetles)

# for Ruben to understand, run abundance models for species: Xarv, Pdet, Ptest, Mneb, Aclav, Cmut

library(unmarked)
#library(emdbook)
library(MuMIn)
options(scipen=999)


### setup 'unmarked' dataframe using site-specific and survey-specific covariates ###
# survey-specific covariates are: JulianDay (day of year) and Precipitation
# site-specific covariates are Site, Canopy cover, Ground cover, Tree density,
# Number of tree species, Dominant tree species, Average tree diameter,
# Number of stumps, and Volume of downed wood

# bind the survey-specific covariates and standardize them
JulianDay = matrix(cbind(Julian1, Julian2, Julian3, Julian4,
                         Julian5, Julian6, Julian7, Julian8,
                         Julian9, Julian10, Julian11, Julian12))
JulianDay = matrix(scale(JulianDay), nrow=50, byrow=F)

Precipitation = matrix(cbind(Precip1, Precip2, Precip3, Precip4,
                             Precip5, Precip6, Precip7, Precip8,
                             Precip9, Precip10, Precip11, Precip12))
Precipitation = matrix(scale(Precipitation), nrow=50, byrow=F)

# create an unmarked objects for theone of the focal species, Plagionotus detritus

PdetUMF <- with(beetles, {
  unmarkedFramePCount(
    y = cbind(Pdet1, Pdet2, Pdet3, Pdet4, Pdet5, Pdet6, Pdet7, Pdet8, Pdet9, Pdet10, Pdet11, Pdet12),
    siteCovs = data.frame(Site, Lures, scale(Canopy), scale(GroundCover), scale(TreeDens), 
                          scale(NoTreeSp), Dominant.1, scale(AvgDiameter), scale(NoStumps), scale(VolDeadWood1)),
    obsCovs = list(Julian = JulianDay,
                   Precip = Precipitation))
  
})

plot(PdetUMF)

# test whether Poisson or Negative Binomial distribution fits the data better
mP <- pcount(~Julian+Precip ~scale.AvgDiameter. + scale.GroundCover., mixture="P", K=100, data=PdetUMF)
mNB <- pcount(~Julian+Precip ~scale.AvgDiameter. + scale.GroundCover., mixture="NB", K=100, data=PdetUMF)

modelList <- fitList("modP" = mP, "modNB" = mNB)
modSel_Pdet <- modSel(modelList)
modSel_Pdet


### BUILD AND RUN DETECTION MODELS (using a full model for occupancy)
# run several models for detection

m0 <- pcount(~1 ~scale.AvgDiameter.+ scale.GroundCover., mixture="NB", K=100, data=PdetUMF)
m1 <- pcount(~Julian ~scale.AvgDiameter.+ scale.GroundCover., mixture="NB", K=100, data=PdetUMF)
m2 <- pcount(~Precip ~scale.AvgDiameter.+ scale.GroundCover., mixture="NB", K=100, data=PdetUMF)
m3 <- pcount(~Julian+Precip ~scale.AvgDiameter.+ scale.GroundCover., mixture="NB", K=100, data=PdetUMF)

# detection model selection 
modelList <- fitList(Null = m0, "mod1" = m1, "mod2" = m2, "mod3" = m3)
modSel_Pdet <- modSel(modelList, nullmod = 'Null')
modSel_Pdet



### BUILD AND RUN MODELS FOR ABUNDANCE ###

Pdet_Null <- pcount(~Julian+Precip  ~1, mixture="NB", K=100, data=PdetUMF)
Pdet_m0 <- pcount(~Julian+Precip  ~scale.GroundCover., mixture="NB", K=100, data=PdetUMF) 
Pdet_m1 <- pcount(~Julian+Precip  ~scale.AvgDiameter., mixture="NB", K=100, data=PdetUMF)
Pdet_m2 <- pcount(~Julian+Precip  ~scale.AvgDiameter. + scale.GroundCover., 
                  mixture="NB", K=100, data=PdetUMF) 

# rank models by AIC for PDET
modelList_Pdet <- fitList(Null = Pdet_Null, m0 = Pdet_m0, m1 = Pdet_m1, m2 = Pdet_m2)
modSel_Pdet <- modSel(modelList_Pdet, nullmod = 'Null')
modSel_Pdet


# extract coeficients for best model
summary(Pdet_m0)
coef(Pdet_m0)

# estimate model averaged individuals per site
Pdet_abund <- predict(Pdet_m0, type="state")
Pdet_abund
mean(Pdet_abund$Predicted)

# compute overall abundance estimates
lcab <- linearComb(Pdet_m0, c(1, 2), type="state")
lcab <- backTransform(lcab)
lcab

# estimate model averaged probability of detection #
Pdet_det <- predict(Pdet_m0, type="det")
Pdet_det
mean(Pdet_det$Predicted)

# compute overall detection estimates
lcdet <- linearComb(Pdet_m0, c(1, 1, 1), type="det")
lcdet <- backTransform(lcdet)
lcdet

# END MODELS

### PLOTS

Date_grid <- seq(min(PdetUMF@obsCovs$Julian), 
                 max(PdetUMF@obsCovs$Julian), 
                 length = 50)
grid <- expand.grid(Julian = Date_grid, Precip = mean(PdetUMF@obsCovs$Precip))
Date_pred <- predict(Pdet_m0, type = "det", as.data.frame(grid))
str(Date_pred)

plot(Date_grid,
     Date_pred$Predicted, 
     type = "l", 
     lwd = 3,
     ylim = c(0,0.3),
     xlab = "scaled Julian Date", 
     ylab = "Detection probability")
lines(Date_grid, Date_pred$lower)
lines(Date_grid, Date_pred$upper)

# predictions for Precip
Precip_grid <- seq(min(PdetUMF@obsCovs$Precip), 
                max(PdetUMF@obsCovs$Precip), 
                length = 50)
grid <- expand.grid(Julian = mean(PdetUMF@obsCovs$Julian), Precip = Precip_grid)
Precip_pred <- predict(Pdet_m0, type = "det", as.data.frame(grid))
str(Precip_pred)

plot(Precip_grid,
     Precip_pred$Predicted, 
     type = "l", 
     lwd = 3, 
     ylim = c(0,0.3),
     xlab = "scaled Precipitation", 
     ylab = "Detection probability")
lines(Precip_grid, Precip_pred$lower)
lines(Precip_grid, Precip_pred$upper)

### Plot predictions for abundance as a function of covariate Ground Cover
GC_grid <- seq(min(PdetUMF@siteCovs$scale.GroundCover.), 
               max(PdetUMF@siteCovs$scale.GroundCover.), 
               length = 50)
scale.GroundCover. <- GC_grid
GC_pred <- predict(Pdet_m0, type = "state", as.data.frame(GC_grid))
str(GC_pred)

plot(GC_grid,
     GC_pred$Predicted,
     type = "l", 
     lwd = 3,
     xlab = "scaled Ground COver", 
     ylab = "Abundance")

lines(GC_grid, GC_pred$lower)
lines(GC_grid, GC_pred$upper)
