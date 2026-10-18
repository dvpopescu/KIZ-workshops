###########################################################################
####  Species distribution modeling for saproxylic beetles in Romania  ####
###########################################################################

# WORKING DIRECTORY: all paths in this script are relative to the "ensembleSDM" folder
# (the folder that contains OCC/ and ENV/). Open this script from that folder, or use
# Session > Set Working Directory > To Source File Location in RStudio, or:
# setwd("path/to/5. SDM/ensembleSDM")
# biomod2 writes its outputs into a sub-folder named after the species ("cerambyx/"),
# and the csv/png outputs below are also written into this folder.

#to reload a previously saved ensemble model (example from the author's original project):
#myBiomodEM <- get(load('cerambyx/cerambyx.cerambyxFirstModeling.ensemble.models.out'))
#install the packages required for the script
#install.packages(c("biomod2","raster","terra","dplyr","rasterVis","beepr"))
# (rgdal was retired from CRAN in 2023 and is not needed: raster/terra read the .tif files)

#load the libraries
library(biomod2)
library(raster)
library(dplyr)
library(terra)
library(rasterVis)

#check packages
(.packages())

#remove list of global environments

#rm(list=ls())

#load our species data
#in the OCC folder only csv for the species can be found. The csv is composed of species name, latitude and longitude.
#for every entry 1 is attributed to the species column. 1 - true presence, 0 - true absence. No true absence are in our data.


DataSpecies <- read.csv("./OCC/cerambyx.csv")         
head(DataSpecies)

#name of the column with the presence data (1 - presence, 0 - true absence)

myRespName <- 'cerambyx'                                                                                                                 

#presence/absence is loaded as binary where 1 is pres and 0 is abs (the absence needs to be generated)
myResp <- as.numeric(DataSpecies[,myRespName])

#the XY coordinates of species data
myRespXY <- DataSpecies[,c("longitude","latitude")]
myRespXY

###########################################################################
#load the current environmental raster files
#env data for Current and 6 projection for HadGEM3-GC31-LL and 6 for Miroc6
#ssp126 (2021-2040,2041-2060,2061-2080), 
#ssp585 (2021-2040,2041-2060,2061-2080))


######## CURRENT ########
ras_current <- list.files(path = "./ENV/ENV_C", pattern = ".tif", full.names = T);
env_current <- stack(ras_current)
proj4string(env_current)

######## Had ########

#load the future environmental raster files H1
ras_h126_2140 <- list.files(path = "./ENV/ENV_H1", pattern = ".tif", full.names = T);
env_h126_2140 <- stack(ras_h126_2140)


#load the future environmental raster files H2
ras_h126_4160 <- list.files(path = "./ENV/ENV_H2", pattern = ".tif", full.names = T);
env_h126_4160 <- stack(ras_h126_4160)


#load the future environmental raster files H3
ras_h126_6180 <- list.files(path = "./ENV/ENV_H3", pattern = ".tif", full.names = T);
env_h126_6180 <- stack(ras_h126_6180)


#load the future environmental raster files H4
ras_h585_2140 <- list.files(path = "./ENV/ENV_H4", pattern = ".tif", full.names = T);
env_h585_2140 <- stack(ras_h585_2140)


#load the future environmental raster files H5
ras_h585_4160 <- list.files(path = "./ENV/ENV_H5", pattern = ".tif", full.names = T);
env_h585_4160 <- stack(ras_h585_4160)


#load the future environmental raster files H6
ras_h585_6180 <- list.files(path = "./ENV/ENV_H6", pattern = ".tif", full.names = T);
env_h585_6180 <- stack(ras_h585_6180)



######## Miroc ########


#load the future environmental raster files M1
ras_m126_2140 <- list.files(path = "./ENV/ENV_M1", pattern = ".tif", full.names = T);
env_m126_2140 <- stack(ras_m126_2140)


#load the future environmental raster files M2
ras_m126_4160 <- list.files(path = "./ENV/ENV_M2", pattern = ".tif", full.names = T);
env_m126_4160 <- stack(ras_m126_4160)


#load the future environmental raster files M3
ras_m126_6180 <- list.files(path = "./ENV/ENV_M3", pattern = ".tif", full.names = T);
env_m126_6180 <- stack(ras_m126_6180)


#load the future environmental raster files M4
ras_m585_2140 <- list.files(path = "./ENV/ENV_M4", pattern = ".tif", full.names = T);
env_m585_2140 <- stack(ras_m585_2140)


#load the future environmental raster files M5
ras_m585_4160 <- list.files(path = "./ENV/ENV_M5", pattern = ".tif", full.names = T);
env_m585_4160 <- stack(ras_m585_4160)


#load the future environmental raster files M6
ras_m585_6180 <- list.files(path = "./ENV/ENV_M6", pattern = ".tif", full.names = T);
env_m585_6180 <- stack(ras_m585_6180)


###########################################################################


#CALCULATE CORRELATIONS BETWEEN VARIABLES
#pearson's correlation is parametric
#spearman and Kendall are not parametric
#reject variables with > |0.75| If enough variables, reject variables > |0.70|
#around 5-8 parameters is a good number
#split your data in training records and testing records - 2 predictors → 50:50
#2 predictors → 50:50 ... 5 predictors → 67:33 ... >10 predictors → 75:25

corr <- raster::layerStats(env_current, 'pearson', na.rm = T)
corr

#transform correlation results in a data frame
corr.df <- as.data.frame(corr)
head(corr.df)
#transform correlation matrix to distances
var.dist <- abs(as.dist(corr.df))
var.dist

var.dist.data <- as.data.frame(as.matrix(var.dist))
write.csv(var.dist.data, file = "Pearson_correlation_cerambyx.csv")
#calculate dendrogram based on distance (less distance  =  more correlation)
var.cluster  <-  hclust(1 - var.dist)

#plot dendrogram
png("correlation_dendrogram_cerambyx.png", width = 1200, height = 900)
plot(var.cluster)
abline(h = 0.25, lty = 2, lwd = 2) # variables that have a correlation < 0.75
dev.off()

###########################################################################

#number of pseudo absences to be used
p_abs <- nrow(myRespXY)
p_abs <- p_abs * 3
p_abs <- round(p_abs, digits = 0)
p_abs 

#?BIOMOD_FormatingData

#formatting data to make sure they can be used within biomod2
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = env_current,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.strategy = "random",
                                     PA.nb.rep = 10,             #How many sets of pseudo absences (10)
                                     PA.nb.absences = p_abs)     #How many pseudo-absence records




#print formatting data
myBiomodData
#plot(myBiomodData)

#defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()


#computing the models 
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    models = c('GLM','GAM','GBM','RF','MAXNET'),
                                    bm.options = myBiomodOption,
                                    CV.nb.rep = 10,                                               #Number of replicates  (100) 
                                    data.split.perc = 80,
                                    prevalence= 0.5,
                                    metric.eval = c('TSS', 'ROC', 'KAPPA'),
                                    var.import = 10,
                                    do.full.models = FALSE,
                                    # save.output = TRUE,
                                    nb.cpu = 1,
                                    seed.val = 13,                                             #Seed value
                                    modeling.id =  paste(myRespName,"FirstModeling",sep=""))



#modeling summary
##myBiomodModelOut


#get_evaluations(myBiomodModelOut, as.data.frame = TRUE)

#get all models evaluation (GLM, GAM, GBM, RF and MAXNET)
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodModelEval
table_evaluation <- myBiomodModelEval
write.csv(table_evaluation, file = "evaluation_1_cerambyx.csv")


#data frame containing variable importance 
get_variables_importance(myBiomodModelOut)
table_importance <- get_variables_importance(myBiomodModelOut)
write.csv(table_importance, file = "table_importance_cerambyx.csv")

###########################################################################
#ensemble modeling containing information relative to an biomod2 ensemble modeling procedure
#TSS scores range from −1 to 1 (Allouche et al. 2006) . TSS values below 0.5 indicate an inaccurate model, 
#values from 0.6 to 0.8 indicate a fairly accurate model, with values larger than 0.8 indicating good model performance (Coetzee et al.)
#eval.metric.quality.threshold = c(value), selects only the models that exceed the threshold

myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by='all',
  metric.select = c('TSS'),
  metric.select.thresh = c(0.4), # THIS IS USSUALY SET ON 0.7
  metric.eval = c('TSS', 'ROC','KAPPA'),
  var.import = 10,
  nb.cpu = 1,
  seed.val = 13,
  do.progress = TRUE,
  em.algo = c('EMmean','EMmedian','EMwmean')
)


#print summary of modeling ensemble
##myBiomodEM



#data frame containing ensemble models evaluation scores
get_evaluations(myBiomodEM)
evaluation_table <- get_evaluations(myBiomodEM)
write.csv(evaluation_table, file = "evaluation_2_cerambyx.csv")

get_variables_importance(myBiomodEM)
variable_importance_ensamble <- get_variables_importance(myBiomodEM)
write.csv(variable_importance_ensamble, file = "table_importance_2_cerambyx.csv")


###########################################################################
#projection on environmental conditions

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_current,
  proj.name = 'current',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)


myCurrentProj <- get_predictions(myBiomodProj)
##plot(myCurrentProj)


#ensemble forecasting
myBiomodEF <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all',
  binary.meth = 'TSS',
  projection.output = myBiomodProj)


#ensemble forecasting plotting results
##plot(myBiomodEF)


##Future 

########## HAD ###########

###### H1 ######
myBiomodProjH1 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h126_2140,
  proj.name = 'H1',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)


ProjH1 <- get_predictions(myBiomodProjH1)

EnsH1 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH1,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH1)


###### #H2 #######
myBiomodProjH2 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h126_4160,
  proj.name = 'H2',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjH2 <- get_predictions(myBiomodProjH2)

EnsH2 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH2,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH2)


####### H3 #######
myBiomodProjH3 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h126_6180,
  proj.name = 'H3',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjH3 <- get_predictions(myBiomodProjH3)

EnsH3 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH3,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH3)


####### H4 #######
myBiomodProjH4 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h585_2140,
  proj.name = 'H4',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjH4 <- get_predictions(myBiomodProjH4)

EnsH4 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH4,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH4)


####### H5 #######
myBiomodProjH5 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h585_4160,
  proj.name = 'H5',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjH5 <- get_predictions(myBiomodProjH5)

EnsH5 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH5,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH5)


####### H6 #######
myBiomodProjH6 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_h585_6180,
  proj.name = 'H6',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjH6 <- get_predictions(myBiomodProjH6)

EnsH6 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjH6,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjH6)



########## MIROC ###########

####### M1 #######
myBiomodProjM1<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m126_2140,
  proj.name = 'M1',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM1 <- get_predictions(myBiomodProjM1)

EnsM1 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM1,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM1)


####### M2 #######
myBiomodProjM2<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m126_4160,
  proj.name = 'M2',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM2 <- get_predictions(myBiomodProjM2)

EnsM2 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM2,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM2)


####### M3 #######
myBiomodProjM3<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m126_6180,
  proj.name = 'M3',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM3 <- get_predictions(myBiomodProjM3)

EnsM3 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM3,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM3)


####### M4 #######
myBiomodProjM4<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m585_2140,
  proj.name = 'M4',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM4 <- get_predictions(myBiomodProjM4)

EnsM4 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM4,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM4)

####### M5 #######
myBiomodProjM5<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m585_4160,
  proj.name = 'M5',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM5 <- get_predictions(myBiomodProjM5)

EnsM5 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM5,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM5)


####### M6 #######
myBiomodProjM6<- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  new.env = env_m585_6180,
  proj.name = 'M6',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = TRUE,
  build.clamping.mask = TRUE,
  output.format = '.tif',
  do.stack = FALSE,
  seed.val = 13)

ProjM6 <- get_predictions(myBiomodProjM6)

EnsM6 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProjM6,
  models.chosen = 'all',
  metric.binary = 'all',
  projection.output = myBiomodProjM6)


###########################################################################


#optional: play a sound when the script is finished (only if the 'beepr' package is installed)
if (requireNamespace("beepr", quietly = TRUE)) beepr::beep("ready")

