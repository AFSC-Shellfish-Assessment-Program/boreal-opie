## summarize events across each model
## to estimate preindustrial probability of borealization

library(tidyverse)
library(brms)

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1854-1949 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1850-1949 for each model preindustrial and hist_ssp585
# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Record the positive and negative outcomes for each preindustrial run
# for each model-anomaly combination (>= this anomaly or not)

# load sst-borealization brms model
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")

# load CMIP6 anomalies
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea",
         experiment == "piControl") %>% 
  select(year, model, annual.unsmoothed) %>%
  rename(sst.anomaly = annual.unsmoothed)

# estimate borealization for each CMIP6 SST anomaly from model posteriors
cmip.preind.borealization <- data.frame()

for(i in 1:nrow(cmip.anom)){
  
  temp <- data.frame(model = cmip.anom$model[i],
                     borealization_index = posterior_predict(sst_boreal_brm,
                                                             newdata = cmip.anom[i,],
                                                             ndraws = 1)) # reducing to 1 draw each, as we have n = 250 preindustrial observations for each model
  cmip.preind.borealization <- rbind(cmip.preind.borealization,
                                  temp)
  
}

nrow(cmip.preind.borealization) # 5750! 23 models * 250 years 
# # load CMIP6 model weights
# model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# # load estimated warming level timing for each model
# timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch outcomes for preindustrial runs
preindustrial.boreal.outcomes <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1


    # separate model of interest
    pre.temp <- cmip.preind.borealization %>% 
      filter(model == models[i])
    
    # loop through each year of observation
    for(k in 1:nrow(observed.borealization)){ # start k loop (observed draws from brms model)
    # k <- 1
    
    annual.event <- ifelse(pre.temp$borealization_index >= observed.borealization$borealization_index[k], 1, 0)
    # annual.event = 1 if cmip value (pre.temp) > observed, = 0 if cmip < observed
    
    # add to df
    preindustrial.boreal.outcomes <- rbind(preindustrial.boreal.outcomes,
                              data.frame(model = models[i],
                                         period = "preindustrial",
                                         observed.year = observed.borealization$year[k],
                                         
                                         borealization_index = observed.borealization$borealization_index[k],
                                         annual.event = annual.event))
    print(paste("i = ", i, "; k = ", k, sep = "")) 

    
      } # close k loop (ersst years)
    

  
  } # close i loop (models)

nrow(preindustrial.boreal.outcomes) # 4.2M rows!


#  and save

for(i in 1:length(models)){
  
  temp <- preindustrial.boreal.outcomes %>%
    filter(model == models[i]) 
  
  write.csv(temp, file = paste("./output/", models[i], "_preindustrial_boreal_outcomes.csv", sep = ""), row.names = F)
  
  
}

