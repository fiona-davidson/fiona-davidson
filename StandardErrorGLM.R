# Setting Data up for GLM Testing (Confidence Interval Testing) 
# June 2019, Fiona 
# Using BC Hecate extents, three env. Vars. = alkalinity, oxygen, silicate. Extracted raster values in ArcMap

setwd("D:/Fiona/BootstrapTesting/")

# IMPORT LIBRARIES
library('readr')
library(rgdal)
library(raster)
library(Rcpp)
library(dismo)
library(tibble)
library(dplyr)
library(vip)
library(ggplot2)
#install.packages("vip")
library(vip)
if (!requireNamespace("devtools")) {
  install.packages("devtools")
}
library(parallel)
library(snow)

# IMPORT DATA

data <- read.csv("D:/Fiona/BootstrapTesting/BC_Hecate/BC_hecate_species_thinned.csv")
pairs(data[,4:6], cex=0.1, fig=TRUE)
#plot(data$Oxygen, data$Silicate)
#plot(data$Alkalinity, data$Oxygen)

data <- mutate(data, presenceab = as.logical(PrsncAb))

# MODEL FITTING

m1 <- glm(formula = presenceab ~ Alkalinity + Oxygen + Silicate, family= binomial(), data = data)
summary(m1)

## Some data to predict at: 100 values over the range of alkalinity 
#ndata <- with(data, data_frame(Alkalinity = seq(min(Alkalinity), max(Alkalinity), length = 100)))
set.seed(364574)
ndata = data[sample(nrow(data), 100), ]
ndata = data

## Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata, fit = predict(m1, newdata = ndata, type = 'response'))

## Plot
plt <- ggplot(ndata, aes(x = Alkalinity, y = fit)) +
  geom_line() +
  theme_set(theme_minimal()) +
  geom_rug(aes(y = PrsncAb, colour = presenceab), data = data) +
  scale_colour_discrete(name = 'Presence Absence') +
  labs(x = 'Alkalinity (umol l-1)', y = 'Probability of species')
plt

plt <- ggplot(ndata, aes(x = Oxygen, y = fit)) +
  geom_line() +
  theme_set(theme_minimal()) +
  geom_rug(aes(y = PrsncAb, colour = presenceab), data = data) +
  scale_colour_discrete(name = 'Presence Absence') +
  labs(x = 'Oxygen (ml l-1)', y = 'Probability of species')
plt

plt <- ggplot(ndata, aes(x = Silicate, y = fit)) +
  geom_line() +
  theme_set(theme_minimal()) +
  geom_rug(aes(y = PrsncAb, colour = presenceab), data = data) +
  scale_colour_discrete(name = 'Presence Absence') +
  labs(x = 'Silicate (umol l -1)', y = 'Probability of species')
plt

#___________________calculating confidence intervals for models_________________#


# Adding confidence interval to model 
fam <- family(m1)
fam
str(fam)
ilink <- fam$linkinv
ilink
poisson()$linkinv

## grad the inverse link function 
ilink <- family(m1)$linkinv 

## add fit and se.fit on the link scale 
ndata <- bind_cols(ndata, setNames(as_tibble(predict(m1, ndata, se.fit = TRUE)[1:2]), 
                                   c('fit_link','se_link'))) 

## create the interval and backtransform
ndata <- mutate(ndata, 
                fit_resp = ilink(fit_link), 
                right_upr = ilink(fit_link + (2 * se_link)), 
                right_lwr = ilink(fit_link - (2 * se_link))) 

## show 
ndata

## Plot
plt + geom_ribbon(data = ndata,
                  aes(ymin = right_lwr, ymax = right_upr),
                  alpha = 0.1)

# ___________________________Predict as raster_________________________#

require(graphics)

# Read rasters into stack 
bcEnv.dir <- ("D:/Fiona/BootstrapTesting/env/")
results.dir <- ("D:/Fiona/BootstrapTesting/results_thinned/")
bootstrap.dir = ("D:/Fiona/BootstrapTesting/boot/")

files_bc <- list.files(path = bcEnv.dir, pattern = "\\.tif$")
s_bc <- stack(paste(bcEnv.dir, files_bc[1], sep=""))

for (i in 2:length(files_bc))
{
  s_bc = addLayer(s_bc, raster(paste(bcEnv.dir, files_bc[i], sep="")))
}
names(s_bc) = c("Alkalinity", "Oxygen", "Silicate")


# Predict to data points
pred_data <- predict(m1, newdata = ndata, type = 'response', se.fit = TRUE)

summary(pred_data$fit)
summary(pred_data$se.fit)
#summary(ndata$fit)

# Predict to raster
predfun <- function(model, data) {
  v <- predict(model, data, se.fit=TRUE)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}

predfun <- function(model, data) {
  v <- predict(model, data, type='response', se.fit=TRUE)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}


plot(pred_map)
plot(pred_map_se)
writeRaster(pred_map, paste(results.dir, "Dec2pred.tif", sep=""), format="GTiff", overwrite=T)
writeRaster(pred_map_se, paste(results.dir, "Dec2pred_se.tif", sep=""), format="GTiff", overwrite=T)

pred_map_nose = predict(s_bc, m1)


# ___________________________Bootstrapping_________________________#

createBootstrapSamples = function(original_samples, nboot)
{
  # Thanks to Jessica Nephin for providing awesome code for this function
  s = sample(row.names(original_samples), nrow(original_samples)*nboot, replace=TRUE)
  boot_s = matrix(s, nrow=nrow(original_samples), ncol=nboot)
  return (boot_s)
}
calibrateModels = function(boot)
{
  dataSample = data[bootsamples[,boot],]
  m1 <- glm(formula = PrsncAb ~ Alkalinity + Oxygen + Silicate, data = dataSample)
  return(m1)
}
applyModels = function(boot)
{
  modelToUse = modelList[[boot]]
  prediction_raster = raster::predict(s_bc, modelToUse)
  return(prediction_raster)
}

se <- function(x) sqrt(var(x)/length(x))

nboot = 200
set.seed(67563)
cat(paste("1-3, Time:", Sys.time(), "\n"))
# Create bootstrap samples of origin_sponges
bootsamples = createBootstrapSamples(data, nboot)

# Calibrate models based on the bootstrap samples. Note that each entry in modelList contains all models in "models"
ncores = detectCores() - 4
cl = parallel::makeCluster(ncores)  
facVars = ""
parallel::clusterExport(cl, varlist=c("bootsamples", "data", "modelList", "s_bc", "facVars"))
modelList = parallel::parLapply(cl, 1:nboot, calibrateModels) # This is a list of the models calibrated with the bootstrap samples
parallel::stopCluster(cl)

# Use calibrated models to make predictions for target.
cl = parallel::makeCluster(ncores)  
facVars = ""
parallel::clusterExport(cl, varlist=c("s_bc", "facVars", "modelList"))
extrapolationList = parallel::parLapply(cl, 1:nboot, applyModels) # This is a list of the maps produced with the bootstrap samples
parallel::stopCluster(cl)


predictionStack = stack()
for (i in 1:nboot)
{
  predictionStack = stack(predictionStack, extrapolationList[[i]])
}

# And then, also for each model, calculate mean, sd, and CoV of extrapolations
meanPred = mean(predictionStack)
beginCluster(ncores)
sdPred = raster::calc(predictionStack, fun=sd)
sePred = raster::calc((predictionStack*sqrt(200)), fun=se)
endCluster()

# Write to file
writeRaster(meanPred, paste(results.dir, "glm_meanboot200.tif", sep=""), format="GTiff", overwrite=T)
writeRaster(sdPred, paste(results.dir, "Dec2glm_sdboot200.tif", sep=""), format="GTiff", overwrite=T)
writeRaster(sePred, paste(results.dir, "Dec2glm_seboot200.tif", sep=""), format="GTiff", overwrite=T)


# ___________________________Running all Biomod Models & Calculating SE_________________________#

library(biomod2);
library(raster);
library(RColorBrewer);
library(dismo);
library(rgdal);
library(ggplot2);
library(dplyr)

setwd("D:/Fiona/DFO_2/BC_Hecate")
sponge_data <- read.csv(file="BC_hecate_species_thinned.csv",header=T)
head(sponge_data)

############################################
#Preparations for Formating 

# Name of studied species
spc_name <- 'GlassSponge'

# Presence/Absences data for the species
myResp<- as.numeric(sponge_data[,'PrsncAb'])
lat_lon<- (sponge_data[,2:3])


intersect_mask = function(x)
{
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

target_predictors = stack(mask(s_bc, intersect_mask(s_bc)))

target_predictors

plot(target_predictors$Alkalinity)

############################################
# BIOMOD STEPS

# 1. Formating data file for Biomod2
setwd("D:/Fiona/DFO_2/BC_Hecate/")
bmData <- BIOMOD_FormatingData(resp.name = spc_name,
                               resp.var = myResp,
                               resp.xy = lat_lon, 
                               expl.var = target_predictors,
                               PA.nb.rep=0,
                               PA.strategy = 'random'
                               
);
bmData

plot(bmData)

# 2. Defining Models Options using default options
myBiomodOption <- BIOMOD_ModelingOptions()
myBiomodOption

# 3. Running models
setwd("D:/Fiona/DFO_2/BC_Hecate/")

# Get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
# Print the dimnames of this object
dimnames(myBiomodModelEval)

## Print metric scores of all selected models
# KAPPA 
myBiomodModelEval["KAPPA","Testing.data",,,]

# TSS
myBiomodModelEval["TSS","Testing.data",,,]

# ROC
myBiomodModelEval["ROC","Testing.data",,,]

# BIAS
myBiomodModelEval["BIAS","Testing.data",,,]

# Relative importance of the explanatory variables
# Print variable importances
var.imp <- get_variables_importance(myBiomodModelOut, nb_rand=1)
var.imp

### Projections 

rasterOptions(tmpdir='D:/Fiona/temp')

intersect_mask = function(x)
{
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

target_predictors = stack(mask(s_bc, intersect_mask(s_bc)))

target_predictors

myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = target_predictors,
                                  proj.name = 'suitability',
                                  selected.models = 'all',
                                  binary.meth = 'ROC',
                                  filtered.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = T, 
                                  do.stack = T,
                                  output.format = '.grd')

#plot(myBiomodProj, str.grep = 'GlassSponge_AllData_RUN1_GAM')
#plot(myBiomodProj, str.grep = 'GlassSponge_AllData_RUN1_GBM')
#plot(myBiomodProj, str.grep = 'GlassSponge_AllData_RUN1_MAXENT.Phillips')


## Load the projections of the models
proj.stk <- get_predictions(myBiomodProj)
## Identify which layer corresponds to which model
names(proj.stk) 
## Choose which layer you want to export, usually 1 = GAM, 2 = GBM and 3 = MaxEnt.
writeRaster(subset(proj.stk, 1), filename="GLM.tif")
writeRaster(subset(proj.stk, 2), filename="GAM.tif")
writeRaster(subset(proj.stk, 3), filename="GBM.tif")
writeRaster(subset(proj.stk, 4), filename="FDA.tif")
writeRaster(subset(proj.stk, 5), filename="MARS.tif")
writeRaster(subset(proj.stk, 6), filename="RF.tif")
writeRaster(subset(proj.stk, 7), filename="MaxEntP.tif")
writeRaster(subset(proj.stk, 8), filename="MaxEntT.tif")

# Calculating Standard Error across above predictions

# A) Read rasters into a stack

model_output.dir <- ("D:/Fiona/DFO_2/BC_Hecate/ModelOutputsThinned/")
model_outputs <- list.files(path = model_output.dir, pattern = "\\.tif$")
S_models <- stack(paste(model_output.dir, model_outputs[1], sep=""))

for (i in 2:length(model_outputs))
{
  S_models = addLayer(S_models, raster(paste(model_output.dir, model_outputs[i], sep="")))
}


# B) Specify RMSE and SE functions

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm=TRUE))

se <- function(x) { sd(x, na.rm=TRUE) / sqrt(length(x[!is.na(x)])) }
#se <- function(x) sqrt(var(x)/length(x))

meanPred = mean(S_models)
beginCluster(ncores)
#rmsePred = raster::calc(S_models, fun=rmse)
#sdPred = raster::calc(S_models, fun=sd)
se_modelsPred = raster::calc(S_models, fun=se)
endCluster()

# Write to file
#writeRaster(meanPred, paste(results.dir, "glm_meanboot200.tif", sep=""), format="GTiff", overwrite=T)
#writeRaster(sdPred, paste(model_output.dir, "stage3_SD.tif", sep=""), format="GTiff", overwrite=T)
writeRaster(se_modelsPred, paste(model_output.dir, "stage3b_SE.tif", sep=""), format="GTiff", overwrite=T)
