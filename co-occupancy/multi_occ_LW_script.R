# R code for co-occupancy analysis of lynx and wolf in the Romanian Carpathians. This analysis uses data from 
# a multi-species (wolf, lynx and wildcat) occupancy analysis is published in Dyck et al., 2022

# Dyck, M. A., Iosif, R., Promberger–Fürpass, B., & Popescu, V. D. (2022). Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians. *Ecology and evolution*, 12(5), e8921.


# Summary -----------------------------------------------------------------

# This code is a simplified version of the code the authors used for Dyck et al., 2022. Current script contains analyses for the winter season only. Autumn season was conducted following the same workflow and can be duplicated using the code below and the autumn data provided. Some exploratory analyses and variables determined to be not useful are not included. We do not include full model selection for detection function but if you would like to re-create the full model selection for detection please contact the authors for a list of models tested. Please reference the README.md file for information about the data files. This file can be opened and knit directly in R Studio using R Markdown or viewed on GitHub.

# Libraries ---------------------------------------------------------------

# install libraries first using install.packages

library(unmarked)
library(dplyr)
library(ggplot2)
library(PerformanceAnalytics)


# Source ------------------------------------------------------------------

# this section uses the package 'unmarked'

# source code provided from Ken Kellner to fix package bugs related to the predict function
source('om_predict_fix.R')


# Data (winter) -------------------------------------------------------------

# read in data

# detection history
winter_spp <- 
  read.csv('lynx_wolf_matrix_winter.csv') %>% 
  
  # set Trapcode as factor
  mutate(TrapCode = as.factor(TrapCode))


# trap effort/observation covaraites
winter_traps <- 
  read.csv('trap_effort_winter.csv') %>% 
  
  # set Trapcode as factor
  mutate(TrapCode = as.factor(TrapCode))


# site covariates
winter_sites <- 
  read.csv('cams_data_winter.csv') %>% 
  
  # alter variable structure
  mutate(TrapCode = as.factor(TrapCode),
         Z = as.numeric(Z),
         Impact = as.factor(Impact),
         CLC2018 = as.factor(CLC2018))


# Format data (winter) ----------------------------------------------------

# this section uses package 'dplyr'

#creating and unmarkedFrameOccuMulti
#unmarkedFrameOccuMulti(y, siteCovs=NULL, obsCovs=NULL, mapInfo)

#creating y: A list (optionally a named list) of length S where each element is an MxJ matrix of the detection, non-detection data for one species, where M is the number of sites, J is the maximum number of sampling periods per site, and S is the number of species in the analysis.
y_multi_winter <- 
  list(
    matrix(unlist(winter_spp[2:9]),   ncol = 8,  byrow = F), # Lynx
    matrix(unlist(winter_spp[10:17]), ncol = 8,  byrow = F) # Wolf
    ) 

# add species names
names(y_multi_winter) <- c("Lynx", "Wolf")

# check that the columns are correct using head or print functions
print(y_multi_winter$Lynx)

# observation covariates/trap effort
obs_covs_winter <- 
  list(winter_traps[, 2:9])

# site covariates
# combining impact at site for detection covariate and calling it 'some'. Not enough samples per each category to keep them separate, combining also simplifies variable reducing potential bias given that some rangers may report "low" while others more specifically report "logging" etc.
some <- c("Isolated buildings", "Logging", "Village", "Low")

# Combine specific CORINE Land Cover types into broader categories (CLC_ag, CLC_open, CLC_forest)
winter_sites <-
  winter_sites %>% 
  mutate(CLC_ag = (CLC211 + CLC222 + CLC242 + CLC243),# proportion agriculture
         CLC_open = (CLC231 + CLC321 + CLC322 + CLC324), # proportion open habitat
         CLC_forest = (CLC311 + CLC312 + CLC313),# proportion CLC
         Impact.2 = case_when(Impact == "None" ~ 0,
                              Impact %in% some ~ 1),
         
         # change Impact.2 to factor
         Impact.2 = as.factor(Impact.2))


# create new object winter_sites.scaled for scaled variables 
winter_sites.scaled <- winter_sites %>% 
  
  # use mutate with across and where to scale all numeric varaibles
  mutate(across(where(is.numeric), scale))

# check data
head(winter_sites.scaled)

# create site covariates data for the unmarkedFrameOccuMulti 
winter_site_covs <- 
  data.frame(winter_sites.scaled)

# adding lynx presence/absence as potential covariate for detection of other species
Lynx <- data.frame(winter_spp[2:9]) # pull lynx detection history out of full detection history

# adding wolf presence/absence as potential covariate for detection of other species
Wolves <- data.frame(winter_spp[10:17])


# combine observation covariates plus species detection histories
winter_obs_covs_species <- 
  c(list(effort = data.frame(winter_traps[,2:9])),
    list(Lynx = Lynx),
    list(Wolf = Wolves))

# create unmarkedFrameOccuMulti object for analysis
winter_occ_data <- 
  unmarkedFrameOccuMulti(
    y_multi_winter, 
    siteCovs = winter_site_covs, 
    obsCovs = winter_obs_covs_species)


# Explore data (winter) ---------------------------------------------------

summary(winter_occ_data)
# naive occupancy for each species (sites with at least 1 detection / total sites)
# lynx 0.67
# wildcat 0.34
# wolf 0.51

plot(winter_occ_data)

# Look at f parameter design matrix
winter_occ_data@fDesign


# Correlations (winter) ---------------------------------------------------

# create a subset of numeric covariates we are interested in using for analysis to test for correlations between variables
winter_sites.corr <- 
  winter_sites.scaled %>%
  select(denslocalr, 
         distlocalr, 
         distnatlro, 
         distsettle, 
         diststream, 
         TRI5, 
         CLC_ag, 
         CLC_open, 
         CLC_forest, 
         CLC311, 
         CLC312, 
         CLC313, 
         Z)

# correlation matrix - Pearson
chart.Correlation(winter_sites.corr, 
                  histogram = TRUE, 
                  method = "pearson")


# use function cor.test() to get Pearson's correlation for highly correlated variables >= 0.70
cor.test(winter_sites.scaled$distnatlro, winter_sites.scaled$distsettle) # 0.8242827 
cor.test(winter_sites.scaled$CLC_ag, winter_sites.scaled$CLC_forest) # -0.7000685 
cor.test(winter_sites.scaled$distsettle, winter_sites.scaled$CLC312) #0.7533168 


# Detection function (winter) ---------------------------------------------

winter_LW_det <- c(
  '~Wolf + diststream', 
  '~Lynx + diststream') 


# Species marginal occupancy (winter) ---------------------------------------------
# I ran these as species specific analyses first, holding marginal occupancy covariates of other 2 species constant as well as co-occupancy covariates constant to see what covariates were important for each species. Variables were selected from full variable set (not included) based on a-priori hypotheses. 

# Lynx (winter) --------------------------------------------------------------------

# lynx TRI5
winter_L_OF_TRI5 <- c('~TRI5', '~1', # marginal occupancy lynx, wolf
                      '~1')    # co-occupancy lynx:wolf

winter_L_TRI5 <- occuMulti(winter_LW_det, winter_L_OF_TRI5, winter_occ_data)
summary(winter_L_TRI5)

# lynx elevation
winter_L_OF_z <- c('~Z', '~1',
                   '~1')

winter_L_z <- occuMulti(winter_LW_det, winter_L_OF_z, winter_occ_data)
summary(winter_L_z)

# lynx CLC
winter_L_OF_CLC <- c('~CLC_forest', '~1',
                     '~1')
                     
winter_L_CLC <- occuMulti(winter_LW_det, winter_L_OF_CLC, winter_occ_data)
summary(winter_L_CLC)

# lynx human disturbance/local rd density
winter_L_OF_rd <- c('~denslocalr', '~1',
                    '~1') 
                    
winter_L_rd <- occuMulti(winter_LW_det, winter_L_OF_rd, winter_occ_data)
summary(winter_L_rd)

# create fitlist for model selection table
winter_L <- fitList(winter_L_CLC, winter_L_TRI5, winter_L_z, winter_L_rd)

# run model selection
modSel(winter_L)


# Wolf (winter) --------------------------------------------------------------------

# wolf road density
winter_W_OF_rd <- c('~1', '~denslocalr', 
                    '~1')
                    
winter_W_rd <- occuMulti(winter_LW_det, winter_W_OF_rd, winter_occ_data)
summary(winter_W_rd)

# wolf CLC
winter_W_OF_CLC <- c('~1', '~CLC_forest',
                     '~1')

winter_W_CLC <- occuMulti(winter_LW_det, winter_W_OF_rd, winter_occ_data)
summary(winter_W_CLC)

# wolf elevation
winter_W_OF_z <- c('~1', '~Z', 
                   '~1') 
                   
winter_W_z <- occuMulti(winter_LW_det, winter_W_OF_z, winter_occ_data)
summary(winter_W_z)

# wolf TRI5
winter_W_OF_TRI5 <- c('~1', '~TRI5',
                      '~1')  
                      
winter_W_TRI5 <- occuMulti(winter_LW_det, winter_W_OF_TRI5, winter_occ_data)
summary(winter_W_TRI5)

# create fitlist for model selection table
winter_W <- 
  fitList(winter_W_CLC, winter_W_TRI5, winter_W_z, winter_W_rd)

# run model selection
modSel(winter_W)


# Marginal occupancy models (winter) --------------------------------------

winter_OF_1 <- c('~CLC_forest', '~Z',
                 '~1') 
                 
winter_LW_1 <- occuMulti(winter_LW_det, winter_OF_1, winter_occ_data)
summary(winter_LW_1)

winter_OF_2 <- c('~denslocalr', '~Z',
                 '~1')
                 
winter_LW_2 <- occuMulti(winter_LW_det, winter_OF_2, winter_occ_data)
summary(winter_LW_2)

winter_OF_3 <- c('~CLC_forest+denslocalr', '~Z',
                 '~1')
                 
winter_LW_3 <- occuMulti(winter_LW_det, winter_OF_3, winter_occ_data)
summary(winter_LW_3)

winter_OF_4 <- c('~CLC_forest+denslocalr', '~Z+CLC_forest',
                 '~1')

winter_LW_4 <- occuMulti(winter_LW_det, winter_OF_4, winter_occ_data)
summary(winter_LW_4)

# create fitlist
winter_occ_mods <- 
  fitList(winter_LW_1, winter_LW_2, winter_LW_3, winter_LW_4)

# model selection
modSel(winter_occ_mods)


# Co-occupancy models (winter) --------------------------------------------

# felid elevation + habitat
winter_LW_OF_200 <- c('~CLC_forest', '~Z',
                       '~Z')
        
winter_LW_200 <- occuMulti(winter_LW_det, winter_LW_OF_200, winter_occ_data)
summary(winter_LWW_200) 

# habitat
winter_LW_OF_201 <- c('~CLC_forest', '~Z',
                       '~CLC_forest')
                    
winter_LW_201 <- occuMulti(winter_LW_det, winter_LW_OF_201, winter_occ_data)
summary(winter_LW_201) 

# wildcat habitat, lynx:wolf movement
winter_LW_OF_202 <- c('~CLC_forest', '~Z',
                       '~TRI5')
                       
winter_LW_202 <- occuMulti(winter_LW_det, winter_LW_OF_202, winter_occ_data)
summary(winter_LW_202) 

# habitat wildct:wolf roads
winter_LW_OF_203<- c('~CLC_forest', '~Z',
                      '~denslocalr') 
                      
winter_LW_203 <- occuMulti(winter_LW_det, winter_LW_OF_203, winter_occ_data)
summary(winter_LW_203)

#lynx:wildcat roads + habitat
winter_LW_OF_204 <- c('~CLC_forest', '~Z',
                       '~denslocalr+CLC_forest')
                       
winter_LW_204 <- occuMulti(winter_LW_det, winter_LW_OF_204, winter_occ_data)
summary(winter_LW_204)

#lynx:wildcat movement, wolf habitat
winter_LW_OF_205 <-  c('~CLC_forest', '~Z',
                        '~denslocalr+CLC_forest+TRI5')

winter_LW_205 <- occuMulti(winter_LW_det, winter_LW_OF_205, winter_occ_data)
summary(winter_LW_205)

#lynx:wildcat movement, wolf:lynx ~ CLC, wildcat:wolf ~elevation
winter_LW_OF_206 <-  c('~CLC_forest', '~Z',
                       '~denslocalr+CLC_forest+Z') 
                        
winter_LW_206 <- occuMulti(winter_LW_det, winter_LW_OF_206, winter_occ_data)
summary(winter_LW_206)

#movement 
winter_LW_OF_207 <-  c('~CLC_forest', '~Z',
                       '~CLC_forest+Z+TRI5')
winter_LW_207 <- occuMulti(winter_LW_det, winter_LW_OF_207, winter_occ_data)
summary(winter_LW_207)

# create fitlist
winter_co.occ_mods <- 
  fitList(winter_LW_200, winter_LW_201, winter_LW_202, winter_LW_203, 
          winter_LW_204, winter_LW_205, winter_LW_206, winter_LW_207)

# model selection
modSel(winter_co.occ_mods)


# Final model (winter) ----------------------------------------------------

#lynx:wildcat movement, wolf habitat
winter_LW_OF_203<- c('~CLC_forest', '~Z',
                     '~denslocalr') 

winter_LW_203 <- occuMulti(winter_LW_det, winter_LW_OF_203, winter_occ_data)
summary(winter_LW_203)


# Marginal occupancy results (winter) -------------------------------------

# this section uses package ggplot2

# calculate marginal occupancy using the predict function

# lynx
winter_L_mo<- 
  (predict(winter_LW_203, #the model we want to use to make predictions
           'state', #specifies we want to predict the occupancy state
           species ='Lynx')) #specifies which species

# saving predicted values for each species to plot
winter_L_pred <- 
  mean(winter_L_mo$Predicted)

winter_L_low <- 
  mean(winter_L_mo$lower)

winter_L_up <- 
  mean(winter_L_mo$upper)


# wolf
winter_W_mo <- 
  (predict(winter_LW_203,
           'state',
           species='Wolf'))

winter_W_pred <- 
  mean(winter_W_mo$Predicted)

winter_W_low <- 
  mean(winter_W_mo$lower)

winter_W_up <- 
  mean(winter_W_mo$upper)

# save predicted values for all three species and assign variable names to plot together
winter_pred <- 
  c(winter_L_pred ,winter_W_pred)

winter_upper <- 
  c(winter_L_up, winter_W_up)

winter_lower <- 
  c(winter_L_low, winter_W_low)

# create variables names for graph 
winter_species <- 
  c("Lynx", "Wolf")

winter_season <- 
  c("Winter", "Winter")

# combine into data frame
winter_pred.df<-
  data.frame(cbind(winter_species, winter_season, winter_pred, winter_upper, winter_lower))

# change variable structure
winter_pred.df$winter_pred <- 
  as.numeric(winter_pred.df$winter_pred)

winter_pred.df$winter_upper <- 
  as.numeric(winter_pred.df$winter_upper)

winter_pred.df$winter_lower <- 
  as.numeric(winter_pred.df$winter_lower)

ggplot(
  data = winter_pred.df, 
  aes(x = winter_species, 
      y = winter_pred))+ 
  geom_point()+
  geom_errorbar(
    data = winter_pred.df, 
    aes(x = winter_species, 
        ymin = winter_lower, 
        ymax = winter_upper), 
    width = 0.25,
    color = '#003E83')+
  scale_y_continuous(
    breaks = seq(0,1,0.1),
    limits = c(0,1),
    expand = c(0,0))+
  labs(
    title = 'Predicted marginal occupancy for lynx and wolf in Winter',
    x = 'Species',
    y = 'Occupancy')+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5,
                              size = 18),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15))


# Marginal detection results (winter) -------------------------------------

# calculate marginal detection using the predict function

# lynx
winter_L_md <- 
  (predict(winter_LW_203, # the model we want to use to make predictions
           'det', # specifies we want the predicted detection
           species='Lynx')) # specifies the species

mean(winter_L_md$Predicted)

# wolf
winter_W_md <- 
  (predict(winter_LW_203,
           'det',
           species ='Wolf'))

mean(winter_W_md$Predicted)



# Conditional occupancy results (winter) ----------------------------------

#calculate conditional occupancy results (predicted occupancy of one species conditional on the presence/absence of another species) using the predict function

# Lynx

# lynx | wolf present
winter_L_co_W <- 
  (predict(winter_LW_203,
           'state',
           species='Lynx',
           cond='Wolf')) 

mean(winter_L_co_W$Predicted)

# lynx | wolf absent
winter_L_noW <- 
  (predict(winter_LW_203,
           'state',
           species ='Lynx',
           cond ='-Wolf'))

mean(winter_L_noW$Predicted)


# Wolf

# wolf | lynx present
winter_W_co_L<- 
  (predict(winter_LW_203,
           'state',
           species ='Wolf',
           cond ='Lynx'))

mean(winter_W_co_L$Predicted)

# wolf | lynx absent
winter_W_noL <- 
  (predict(winter_LW_203,
           'state',
           species ='Wolf',
           cond ='-Lynx')) 

mean(winter_W_noL$Predicted)


# Marginal occupancy graphs (winter) --------------------------------------

# this section uses package ggplot2 and package rphylopic

# I have provided code to make a marginal occupancy graph for one species based on one variable. You can use the template below to graph other species predictions and other variables


# data

# create new data frame holding all variables constant except the one we want to plot using expand.grid function. This new data frame must include all variables in the model including both detection on occupancy variables for all species and species combinations
winter_203_Z.df <- 
  data.frame(
    expand.grid(
      Z = seq(
        min(winter_sites.scaled$Z), 
        max(winter_sites.scaled$Z), 
        0.01),
      CLC_forest = mean(winter_sites.scaled$CLC_forest),
      denslocalr = mean(winter_sites.scaled$denslocalr),
      diststream = mean(winter_sites.scaled$diststream)))

# calculate predicted occupancy using the chosen model for a species to add to the new data frame created above
winter_Epsi_203_W <- 
  predict(winter_LW_203, 
          type ="state", 
          species = 'Wolf',  
          newdata = winter_203_Z.df)

# combine predicted data from above with new data frame
winter_Epsi_203_W <- 
  data.frame(winter_Epsi_203_W, 
             denslocalr = winter_203_Z.df$denslocalr,
             Z = winter_203_Z.df$Z,
             CLC_forest = winter_203_Z.df$CLC_forest,
             diststream = winter_203_Z.df$diststream)


winter_Epsi_203_W <- 
  winter_Epsi_203_W %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE)

# add raw data for x-axis to new data frame for plotting instead of scaled variable which is hard to interpret
alt <- 
  matrix(
    seq(from = 663, 
        to = 1600, 
        length.out = 447),
    ncol = 1)

# combined
winter_Epsi_203_W <- 
  cbind(winter_Epsi_203_W, alt)


# graph

winter_W_Z.plot <- ggplot(data = winter_Epsi_203_W, aes(x = alt, 
                                                          y = Predicted)) +
  
  # add predicted line
  geom_line(linewidth = 1, 
            color = "black") +
  
  # add error ribbon around predicted line
  geom_ribbon( aes(ymin = lwr, 
                   ymax = upr), 
               alpha = 0.45, 
               fill = "skyblue3") +
  
  # alter axis labels
  xlab ("Altitude (meters)")+
  ylab ("Occupancy Probability - Wolf") +
  
  # standardize plot area and set breaks
  coord_cartesian(ylim = c(0,1),
                  xlim = c(706,1550)) +
  scale_x_continuous(breaks = seq(700,1600, 200)) +
  
  # adjust theme elements for pub
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 18))

winter_W_Z.plot

# Co-occupancy graphs -----------------------------------------------------

# this sections uses package ggplot 


# data

# create new data frame holding all variables constant except the one we want to plot using expand.grid function. This new data frame must include all variables in the model including both detection on occupancy variables for all species and species combinations
winter_203_denslocalr.df <- 
  data.frame(
    expand.grid(
      denslocalr = seq(
        min(winter_sites.scaled$denslocalr), 
        max(winter_sites.scaled$denslocalr), 
        0.01),
      Z = mean(winter_sites.scaled$Z),
      CLC_forest = mean(winter_sites.scaled$CLC_forest),
      diststream = mean(winter_sites.scaled$diststream)))

# add raw data for x-axis to new data frame for plotting instead of scaled variable which is hard to interpret
roads <- 
  matrix(
    seq(from = 0.21, 
        to = 0.35, 
        length.out = 401),
    ncol = 1)

# get predicted values for co-occurrence of lynx and wolf with new data frame only varying values for CLC_forest
winter_Epsi_coOcc_LW<- 
  predict(winter_LW_203, 
          type = "state", 
          species = c('Lynx', 'Wolf'), 
          newdata = winter_203_denslocalr.df)

# combine predicted data from above with new data frame
winter_Epsi_coOcc_LW <- 
  data.frame(winter_Epsi_coOcc_LW,
             denslocalr = winter_203_denslocalr.df$denslocalr,
             Z = winter_203_denslocalr.df$Z,
             CLC_forest = winter_203_denslocalr.df$CLC_forest)

# combine raw CLC_forest values with predicted values for lynx and wolf in new data frame
winter_Epsi_coOcc_LW <- 
  cbind(winter_Epsi_coOcc_LW, roads)

winter_Epsi_coOcc_LW <- 
  winter_Epsi_coOcc_LW %>% 
  mutate(lwr = Predicted - SE,
         upr = Predicted + SE)


# graph
winter_L_W_coOcc.plot <- ggplot(data = winter_Epsi_coOcc_LW, aes(x = roads, 
                                                                 y = Predicted)) + 
  
  # add predicted line
  geom_line(size = 1) +
  
  # add error ribbon
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr), 
              alpha = .45, 
              fill = "skyblue3") +
  
  # alter axis labels
  labs(x = expression ("Road density (km/sqkm)"),
       y = "Predicted Co-occurrence") +
  
  # standardize plot area
  coord_cartesian(ylim = c(0,0.955),
                  xlim = c(0.2,0.36)) +
  
  # title with species common names
  ggtitle("Lynx & Wolf") +
  
  # adjust theme elements for pub
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 18),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18))
  

winter_L_W_coOcc.plot


# Conditional occupancy graph ---------------------------------------------

# data

# can use the new data frame created earlier for roads
roads <- 
  matrix(
    seq(from = 0.21, 
        to = 0.35, 
        length.out = 401),
    ncol = 1)

# create variable for predicted occupancy of each species conditional on every other species presence/absence

# lynx | wolf
winter_Epsi_L_W <- 
  predict(winter_LW_203, 
          type ="state", 
          species = 'Lynx', 
          cond = 'Wolf', 
          newdata = winter_203_denslocalr.df)

# lynx | wolf absent
winter_Epsi_L_noW <- 
  predict(winter_LW_203, 
          type ="state", 
          species = 'Lynx', 
          cond = '-Wolf', 
          newdata= winter_203_denslocalr.df)


# wolf | lynx
winter_Epsi_W_L <- 
  predict(winter_LW_203, 
          type ="state", 
          species = 'Wolf', 
          cond = 'Lynx', 
          newdata = winter_203_denslocalr.df)

# wolf | lynx absent
winter_Epsi_W_noL <- 
  predict(winter_LW_203, 
          type ="state", 
          species = 'Wolf', 
          cond = '-Lynx', 
          newdata = winter_203_denslocalr.df)


# add conditional column for each of the conditional occupancy data frames that specifies present/absent and another for species, as well as column for 1 SE, and column for the unscaled variable CLC_forest

# lynx | wolf
winter_Epsi_L_W <- 
  winter_Epsi_L_W %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Present",
         c.spp = "Wolf")

winter_Epsi_L_W <- 
  cbind(winter_Epsi_L_W, roads)

# lynx | wolf absent
winter_Epsi_L_noW <- 
  winter_Epsi_L_noW %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Absent",
         c.spp = "Wolf")

winter_Epsi_L_noW <- 
  cbind(winter_Epsi_L_noW, roads)

# wolf | lynx
winter_Epsi_W_L <- 
  winter_Epsi_W_L %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Present",
         c.spp = "Lynx")

winter_Epsi_W_L <- 
  cbind(winter_Epsi_W_L, roads)

# wolf | lynx absent
winter_Epsi_W_noL <- 
  winter_Epsi_W_noL %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Absent",
         c.spp = "Lynx")

winter_Epsi_W_noL <- 
  cbind(winter_Epsi_W_noL, roads)



winter_Epsi_condOcc_all <- 
  rbind(winter_Epsi_L_W, 
        winter_Epsi_L_noW, 
        winter_Epsi_W_L, 
        winter_Epsi_W_noL) 

View(winter_Epsi_condOcc_all)


# graph

# labels and color vectors
colors <- 
  c("plum4", "lightseagreen")

label <- 
  c("Absent", "Present")

lines <- 
  c("dotted", "solid")

occ.labels <- 
  c("Probability of lynx", 
    "Probability of wolf")

names(occ.labels) <- 
  c("Lynx", 
    "Wolf")

cond.labels <- 
  c("Conditional on lynx", 
    "Conditional on wolf")

names(cond.labels) <- 
  c("Lynx",
    "Wolf")

winter_condOcc.plot <- ggplot(data = winter_Epsi_condOcc_all, aes(x = roads, 
                                                                  y = Predicted, 
                                                                  group = conditional)) +
  
  # add error ribbon
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr, 
                  fill = conditional), 
              alpha = 0.7) +
  
  # add predicted lines
  geom_line(aes(x = roads, 
                y = Predicted, 
                linetype = conditional), 
            size = 0.5, 
            color = "black") +
  
  # rename axis labels
  labs(x = expression ("Road density (km/sqkm)"),
       y = "Occupancy Probability") +
  
  # use facet grid to divide the plots and label the panels
  facet_grid(c.spp ~ spp, 
             labeller = labeller(c.spp = cond.labels, 
                                 spp = occ.labels)) +
  
  # change the line type to be different for present and absent data
  scale_linetype_manual(values = lines, 
                        labels = label) +
  
  
  # change colors of error ribbon to be different for present and absent data
  scale_fill_manual(values = colors, 
                    labels = label) +
  
  # set plotting area
  coord_cartesian(ylim = c(0,1.1)) +
  
  # specify breaks for y -axis
  scale_y_continuous(breaks = seq(from = 0, to =1, by = 0.25)) +
  
  # adjust theme elements for pub
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 15)); winter_condOcc.plot #print

# columns are spp it is predicting occupancy for, and rows are the species the predictions are conditional on

