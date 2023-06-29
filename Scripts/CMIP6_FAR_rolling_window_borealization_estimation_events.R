## summarize events across each model-region combination
## to estimate historical probability

library(tidyverse)
library(brms)

theme_set(theme_bw())

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1854-1949 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1850-1949 for each model preindustrial and hist_ssp585
# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Record the positive and negative outcomes for each preindustrial run
# for each model-anomaly combination (>= this anomaly or not)
# this was done in "CMIP6_FAR_estimation_events.R"

## STEP 4 ---------------------------------------------
# Record the positive and negative outcomes for each hist_ssp585 run
# for each model-anomaly combination (>= this anomaly or not)
# using 15-year rolling window of estimated warming as "present"
# THIS SCRIPT

# load sst-borealization brms model
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")

# load observed borealization values (drawn from brms posteriors in "CMIP6_FAR_borealization_estimation_events.R")
observed.borealization <- read.csv("./output/observed_borealization_from_posteriors.csv")

# ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")
# 
# # limit to 1950-2022 (period of interest for attribution)
# ersst.anom <- ersst.anom %>%
#   filter(year %in% 1950:2022)

# # load CMIP6 anomalies
# cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv")

# load CMIP6 anomalies
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea",
         experiment == "hist_ssp585",
         year %in% 1940:2040) %>% # (to capture warming that matches period of interest for attribution)
  select(year, model, annual.unsmoothed) %>%
  rename(sst.anomaly = annual.unsmoothed)

# estimate borealization for each CMIP6 SST anomaly from model posteriors
cmip.hist.borealization <- data.frame()

for(i in 1:nrow(cmip.anom)){
  # i <- 1
  temp <- data.frame(model = cmip.anom$model[i],
                     year = cmip.anom$year[i],
                     borealization_index = posterior_predict(sst_boreal_brm,
                                                             newdata = cmip.anom[i,],
                                                             ndraws = 4))
  cmip.hist.borealization <- rbind(cmip.hist.borealization,
                                     temp)
  
}

# save
write.csv(cmip.hist.borealization, "./output/cmip_historical_borealization_draws.csv", row.names = F)


cmip.hist.borealization <- read.csv("./output/cmip_historical_borealization_draws.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# load CMIP6 model weights
model.weights <- read.csv("./data/normalized_CMIP6_weights.csv") 

# load brms estimate of warming trend for each model
model.warming.trends <- read.csv("./data/CMIP6_brms_warming_rate.csv")

# check model names
check <- data.frame(models = models,
                    trend.names = unique(model.warming.trends$model))
check # line up fine

# and load predicted warming rate across all models 
predicted.warming <- read.csv("./data/brms_predicted_North_Pac_warming.csv")


## record historical outcomes ----------------
# approach for each observation:
# -select window of year-7 to year+7 for each observation
# -find range of predicted warming during that window in predicted.warming
# -limit each model to the relevant warming range in model.warming.trends
# calculate proportion as big as or larger than ersst anomaly


  # create df of historical outcomes
  historical.rolling.window.borealization.outcomes <- data.frame()  

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1


  # separate model of interest
  hist.temp <- cmip.hist.borealization %>% 
      filter(model == models[i])

  
  # loop through each borealization draws from observed SST
  for(k in 1:nrow(observed.borealization)){ # start k loop (draws)
    # k <- 1
    
    # define 15-year window
    window <- (observed.borealization$year[k] - 7) : (observed.borealization$year[k] + 7)
    
    # define range of predicted warming values for this window
    warming.range <- range(predicted.warming$pred_mean[predicted.warming$year %in% window])
    
    # find years for the model of interest that fall into this warming range
    use <- model.warming.trends %>%
      filter(model == models[i],
             warming >= warming.range[1] & warming <= warming.range[2])
    
    hist.temp.use <- hist.temp %>%
      filter(year %in% use$year)
    
    # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    
    annual.events <- ifelse(hist.temp.use$borealization_index >= observed.borealization$borealization_index[k], 1, 0)
    

    # add to df
    if(length(annual.events) > 0) { # in case there are no cmip years within the required warming rate
    
    historical.rolling.window.borealization.outcomes <- rbind(historical.rolling.window.borealization.outcomes,
                                    data.frame(model = models[i],
                                               period = "historical",
                                               ersst.year = observed.borealization$year[k],
                                               
                                               borealization_index = observed.borealization$borealization_index[k],
                                               annual.events = annual.events))
    
    print(paste("i = ", i, "; k = ", k, sep = ""))  
                                             
    }
    
   } # close k loop (ersst years)
  
  
  } # close i loop (models)



# break into separate objects for each model and save

for(i in 1:length(models)){

  temp <- historical.rolling.window.borealization.outcomes %>%
    filter(model == models[i])

  write.csv(temp, file = paste("./output/historical_borealization_outcomes/",
                               models[i], "_historical_borealization_outcomes_rolling_window.csv", sep = ""), row.names = F)

}

