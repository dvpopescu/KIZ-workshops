
### This code analyzes occupancy data for mink frogs (Lithobates septentrionalis) collected 
### at 66 ponds in the Adirondack Mountains, New York State, USA in Summer 2012.

setwd('....')

mink <- read.csv("minkfrog_dataset.csv", header=T)

names(mink)
str(mink)

# we need to 'attach' the file for unmarked to work; usually it is not 
# recommended to use this function
attach(mink) 

library(unmarked)
library(MuMIn)

# calculate naive occupancy (based on observed data, not accounting for imperfect detection)
minkfrog_naive <- sum(apply(mink[,3:6], 1, sum) > 0) / nrow(mink)
minkfrog_naive

# create an unmarked object for the focal species, the mink frog (Lithobates septentrionalis)
# survey-specific covariates are: Date (day of year), Wind and Sky
# site-specific covariates are DO (pond Dissolved Oxygen in mg/l) and 
# JulyTemp (mean historical July temperature in degree Celsius)

minkUMF <- with(mink, {
  unmarkedFrameOccu(
    
    y = cbind(LISE1, LISE2, LISE3, LISE4),
    
    siteCovs = data.frame(scale(JulyTemp), scale(DO), scale(DistMFbreed), Beaver),
    
    obsCovs = list(Date = scale(cbind(Julian1, Julian2, Julian3, Julian4)),
                   Wind = cbind(Wind1, Wind2, Wind3, Wind4),
                   Sky = cbind(Sky1, Sky2, Sky3, Sky4)))
})

# run models to evaluate the best detection variables
mink.d1 <- occu(~Date + Wind + Sky ~scale.JulyTemp. + scale.DO., data=minkUMF)
mink.d2 <- occu(~Date + Sky ~scale.JulyTemp. + scale.DO., data=minkUMF)
mink.d3 <- occu(~Date  ~scale.JulyTemp. + scale.DO., data=minkUMF)
mink.d4 <- occu(~Sky  ~scale.JulyTemp. + scale.DO., data=minkUMF)
mink.d5 <- occu(~Date + Wind  ~scale.JulyTemp. + scale.DO., data=minkUMF)

# model selection for detection models
modlist_d <- list (DateWindSky = mink.d1, DateSky = mink.d2, 
                   Date = mink.d3, Sky = mink.d4, DateWind = mink.d5)

selection_d = model.sel(modlist_d)
selection_d

# run occupancy models using the best detection variables
mink.o1 <- occu(~Date + Sky ~scale.JulyTemp., data=minkUMF)
mink.o2 <- occu(~Date + Sky ~scale.DO., data=minkUMF)
mink.o3 <- occu(~Date + Sky ~scale.JulyTemp. + scale.DO., data=minkUMF)
mink.Null <- occu(~Date + Sky ~1, data=minkUMF)

# model selection for occupancy variables
modlist_o <- list(JTemp = mink.o1, DO = mink.o2, JTempDO = mink.o3,
                  Null = mink.Null)

selection_o = model.sel(modlist_o)
selection_o

# predict occupancy and detection probabilities base don the top model (DO)
detpred = predict(mink.o2, type="det")
detpred
hist(detpred$Predicted)

psipred = predict(mink.o2, type="state")
psipred
hist(psipred$Predicted)

# compute overall occupancy and detection estimates

lcpsi <- linearComb(mink.o2, c(1, 0), type="state")
lcpsi <- backTransform(lcpsi)

lcdet <- linearComb(mink.o2, c(1, 0, 0), type="det")
lcdet <- backTransform(lcdet)


### Plot predictions for probability of detection as a function of covariates Date and Sky
# predictions for Date
Date_grid <- seq(min(minkUMF@obsCovs$Date), 
                 max(minkUMF@obsCovs$Date), 
                 length = 50)
grid <- expand.grid(Date = Date_grid, Sky = mean(minkUMF@obsCovs$Sky))
Date_pred <- predict(mink.o2, type = "det", as.data.frame(grid))
str(Date_pred)

plot(Date_grid,
     Date_pred$Predicted, 
     type = "l", 
     lwd = 3,
     ylim = c(0,1),
     xlab = "scaled Date", 
     ylab = "Detection probability")
lines(Date_grid, Date_pred$lower)
lines(Date_grid, Date_pred$upper)

# predictions for SKY
Sky_grid <- seq(min(minkUMF@obsCovs$Sky), 
                max(minkUMF@obsCovs$Sky), 
                length = 50)
grid <- expand.grid(Date = mean(minkUMF@obsCovs$Date), Sky = Sky_grid)
Sky_pred <- predict(mink.o2, type = "det", as.data.frame(grid))
str(Sky_pred)

plot(Sky_grid,
     Sky_pred$Predicted, 
     type = "l", 
     lwd = 3, 
     ylim = c(0,1),
     xlab = "scaled Sky", 
     ylab = "Detection probability")
lines(Sky_grid, Sky_pred$lower)
lines(Sky_grid, Sky_pred$upper)

### Plot predictions for probability of occupancy as a function of covariate DO
DO_grid <- seq(min(minkUMF@siteCovs$scale.DO.), 
               max(minkUMF@siteCovs$scale.DO.), 
               length = 50)
scale.DO. <- DO_grid
DO_pred <- predict(mink.o2, type = "state", as.data.frame(DO_grid))
str(DO_pred)

plot(DO_grid,
     DO_pred$Predicted,
     type = "l", 
     lwd = 3, 
     ylim = c(0,1),
     xlab = "scaled Dissolved Oxygen", 
     ylab = "Occupancy probability")

lines(DO_grid, DO_pred$lower)
lines(DO_grid, DO_pred$upper)
