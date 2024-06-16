## create borealization PDFs for different climate states

library(tidyverse)
library(brms)


# load sst-borealization brms model
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")

## preindstrial runs---------------------------------------
# load preindustrial CMIP6 anomalies
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
                                                             ndraws = 1)) #  1 draw each, as we have n = 250 preindustrial observations for each model
  cmip.preind.borealization <- rbind(cmip.preind.borealization,
                                  temp)
  
}

nrow(cmip.preind.borealization) # 5750! 23 models * 250 years 

### now historical / ssp558 runs------------------------------
# load CMIP6 anomalies
cmip.anom <- read.csv("./Data/CMIP6.anomaly.time.series.csv") %>%  
  filter(region == "Eastern_Bering_Sea") %>% 
  select(year, experiment, model, annual.unsmoothed) %>%
  rename(sst.anomaly = annual.unsmoothed)

# load estimated warming level timing for each model
timing <- read.csv("./Data/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch outcomes for extreme runs
borealization.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i])
  
  
  # record anomalies
  borealization.pdfs <- rbind(borealization.pdfs,
                        data.frame(model = rep(models[i], 250),
                                   period = rep("preindustrial", 250),
                                   borealization_index = as.vector(posterior_predict(sst_boreal_brm,
                                                                           newdata = pre.temp,
                                                                           ndraws = 1))))
  
  
 


## record outcomes using different warming levels from hist.585 

  
  # separate model 
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i])
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # estimate borealization values
  borealization.pdfs <- rbind(borealization.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   borealization_index = as.vector(posterior_predict(sst_boreal_brm,
                                                                                     newdata = hist.temp.use,
                                                                                     ndraws = 5)))) # using 5 draws for historical periods
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # estimate borealization values
  borealization.pdfs <- rbind(borealization.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   borealization_index = as.vector(posterior_predict(sst_boreal_brm,
                                                                                     newdata = hist.temp.use,
                                                                                     ndraws = 5))))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  borealization.pdfs <- rbind(borealization.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   borealization_index = as.vector(posterior_predict(sst_boreal_brm,
                                                                                     newdata = hist.temp.use,
                                                                                     ndraws = 5))))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  borealization.pdfs <- rbind(borealization.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   borealization_index = as.vector(posterior_predict(sst_boreal_brm,
                                                                                     newdata = hist.temp.use,
                                                                                     ndraws = 5))))
  
  
  
} # close i loop (models)

# save
write.csv(borealization.pdfs, "./output/borealization_pdfs.csv", row.names = F)

## plot --------------------

borealization.pdfs <- read.csv("./output/borealization_pdfs.csv")

# check outcomes per model
check <- borealization.pdfs %>%
  group_by(period, model) %>%
  summarise(count = n())

check # instances differ among models, need to account for this!

# check outcomes per period
check <- borealization.pdfs %>%
  group_by(period) %>%
  summarise(count = n())

check

borealization.pdfs <- left_join(borealization.pdfs, check)

# load CMIP6 model weights
model.weights <- read.csv("./Data/normalized_CMIP6_weights.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(model, normalized_weight) 


borealization.pdfs <- left_join(borealization.pdfs, model.weights) %>%
  mutate(total_weight = normalized_weight/count)

# resample to weight models
resample.pdf <- resample_density <- data.frame()

periods <- unique(borealization.pdfs$period)

set.seed(99)
for(i in 1:length(periods)){
  # i <- 1
  
  temp <- borealization.pdfs[borealization.pdfs$period == periods[i],]
  
  # for(j in 1:1000){ # 1000 density estimates from each period (each defined by 1000 draws from the relevant borealization_index estimates)

  
  temp_borealization_index <- sample(temp$borealization_index, 100000, replace = T, prob = temp$total_weight)
    
  resample.pdf <- rbind(resample.pdf,
  data.frame(period = periods[i],
  borealization_index = temp_borealization_index))

  
  # } # close j loop (density estimates)

} # close i loop (periods)



# save
write.csv(resample.pdf, "./output/resampled_borealization_pdfs.csv", row.names = F)

# load 
resample.pdf <- read.csv("./output/resampled_borealization_pdfs.csv")

# get summary statistics (proportion of years < -1 or > 2 by era)
sum <- resample.pdf %>%
  group_by(period) %>%
  summarize(proportion_high = sum(borealization_index > 2) / n(),
            proportion_low = sum(borealization_index < -1) / n(),
            ratio_high_to_low = proportion_high / proportion_low) 

sum

