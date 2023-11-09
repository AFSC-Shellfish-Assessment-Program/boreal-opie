## estimate preindustrial and historical 
## probabilities for given borealization value

library(tidyverse)


theme_set(theme_bw())

# load observed borealization events 
# load DFA model
mod <- readRDS("./output/DFA_model.rds")
summary(mod)

# calculate DFA trend (borealization index)
observed.borealization <- data.frame(year=1972:2022,
                    borealization_index = as.vector(mod$states))


# load cmip6 historical draws from the same model
historical.borealization <- read.csv("./output/cmip_historical_borealization_draws.csv")

nrow(historical.borealization)

# load preindustrial CMIP6 anomalies and make draws from posteriors
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
    filter(region == "Eastern_Bering_Sea",
           experiment == "piControl") %>% 
    select(year, model, annual.unsmoothed) %>%
    dplyr::rename(sst.anomaly = annual.unsmoothed)

# load brms model
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")

# estimate borealization for each CMIP6 SST anomaly from model posteriors
preindustrial.borealization <- data.frame()

for(i in 1:nrow(cmip.anom)){
    print(i)
    temp <- data.frame(model = cmip.anom$model[i],
                       borealization_index = posterior_predict(sst_boreal_brm,
                                                               newdata = cmip.anom[i,],
                                                               nsamples = 1)) #  1 draw each, as we have n = 250 preindustrial observations for each model
    preindustrial.borealization <- rbind(preindustrial.borealization,
                                       temp)
    
}

nrow(preindustrial.borealization) 

# save
write.csv(preindustrial.borealization,
          "./output/cmip_preindustrial_historical_borealization_draws.csv", row.names = F)

preindustrial.borealization <- read.csv("./output/cmip_preindustrial_historical_borealization_draws.csv")

# load model weights 
weights <- read.csv("./data/normalized_CMIP6_weights.csv")

weights <- weights %>%
    filter(region == "Eastern_Bering_Sea") %>%
   select(model, normalized_weight) %>%
   dplyr::rename(model_weight = normalized_weight)

# save for paper
write.csv(weights, "./output/model_weights.csv", row.names = F)


# load brms estimate of warming trend for each model
model.warming.trends <- read.csv("./data/CMIP6_brms_warming_rate.csv")

# check model names
# get vector of model names
models <- unique(cmip.anom$model)

check <- data.frame(models = models,
                    trend.names = unique(model.warming.trends$model))
check # line up fine

# and load predicted warming rate across all models 
predicted.warming <- read.csv("./data/brms_predicted_North_Pac_warming.csv")

###################
# load estimated warming level timing for each model
timing <- read.csv("./Data/model.north.pacific.warming.timing.csv")


# create df to catch preindustrial and historical outcomes (greater than or not) for each borealization observation
outcomes <- data.frame()

# loop through each observation
for(h in 1:nrow(observed.borealization)){
    # h <- 1

    obs.temp <- observed.borealization[h,]
    
    # define 15-year window
    window <- (obs.temp$year - 7) : (obs.temp$year + 7)
    
    # define range of predicted warming values for this window
    warming.range <- range(predicted.warming$pred_mean[predicted.warming$year %in% window])
    
    
    # loop through each model
for(i in 1:length(models)){ # start i loop (models)
    # i <- 1
    
    # separate model of interest from preindustrial runs
    pre.temp <- preindustrial.borealization %>% 
        filter(model == models[i])
    
    
    # record preindustrial outcomes
    outcomes <- rbind(outcomes,
                                data.frame(year = obs.temp$year,
                                           observed_borealization = obs.temp$borealization_index,
                                            model = models[i],
                                           period = "preindustrial",
                                           n = 250,
                                           outcome = if_else(pre.temp$borealization_index >= obs.temp$borealization_index, 1, 0)))
    
    
    ## record historical outcomes 
    # approach for each observation:
    # -select window of year-7 to year+7 for each observation
    # -find range of predicted warming during that window in predicted.warming
    # -limit each model to the relevant warming range in model.warming.trends
    # calculate proportion as big as or larger than observed borealization
    
    # find years for the model of interest that fall into this warming range
    use <- model.warming.trends %>%
        filter(model == models[i],
               warming >= warming.range[1] & warming <= warming.range[2])
    
    hist.temp.use <- historical.borealization %>%
        filter(model == models[i],
               year %in% use$year)
    
    if(nrow(hist.temp.use) > 0){ # record outcomes only if we have observations in 1940:2040 matching this warming window
    
    outcomes <- rbind(outcomes,
                   data.frame(year = obs.temp$year,
                              observed_borealization = obs.temp$borealization_index,
                              model = models[i],
                              n = nrow(hist.temp.use),
                              period = "historical",
                              outcome = if_else(hist.temp.use$borealization_index >= obs.temp$borealization_index, 1, 0)))
    }
    
} # close i loop (models) 
    
} # close h loop (observations)
    

# save
write.csv(outcomes, "./output/borealization_outcomes_for_attribution.csv", row.names = F)

outcomes <- read.csv("./output/borealization_outcomes_for_attribution.csv")

# join with model weights
outcomes <- left_join(outcomes, weights) %>%
    mutate(total_weight = model_weight/n) 


# examine sample size for each (preindustrial/historical)
outcomes.summary <- outcomes %>%
    group_by(period) %>%
    dplyr::summarise(count = n())

outcomes.summary


historical.summary <- outcomes %>%
    filter(period == "historical") %>%
    group_by(year) %>%
    dplyr::summarise(count = n())

historical.summary

preindustrial.summary <- outcomes %>%
    filter(period == "preindustrial") %>%
    group_by(year) %>%
    dplyr::summarise(count = n())

preindustrial.summary


# dividing weight by n - all models have same n for preindustrial
# outcomes, so no effect; but different n for historical windows, so
# need to account for this in weights when resampling data to avoid
# over-weighting models with more historical observations in a given window

# resample and calculate FAR!

# create df for output
resampled_FAR <- data.frame()


# now resample
periods <- unique(outcomes$period)

set.seed(99)
for(y in 1972:2022){
    # y <- 1972
    for(x in 1:1000){ # go through 1000 iterations
        # x <- 1
    # start with preindustrial probability
     p <- 1
    
     temp <- outcomes %>%
        filter(period == periods[p],
               year == y) 
     
     # resample outcomes 1000 times
     
     preind_resample <- sample(temp$outcome, 1000, replace = T, prob = temp$total_weight)
     preind_prob <- sum(preind_resample)/1000
     
     # now historical probability
     p <- 2
     
     temp <- outcomes %>%
         filter(period == periods[p],
                year == y) 
     
     # resample outcomes 1000 times
     
     hist_resample <- sample(temp$outcome, 1000, replace = T, prob = temp$total_weight)
     
     hist_prob <- sum(hist_resample)/1000
     
     # and record results for this iteration
     resampled_FAR <- rbind(resampled_FAR,
                            data.frame(year = y,
                                       draw = x,
                                       preind_prob = preind_prob,
                                       hist_prob = hist_prob,
                                       FAR = 1 - (preind_prob/hist_prob)))
     
     
    }
    
}



# and summarize results
attribution_stats <- data.frame()

for(y in 1972:2022){
    
temp <- resampled_FAR %>%
    filter(year == y)

attribution_stats <- rbind(attribution_stats,
                           data.frame(year = y,
                                      preind_prob = mean(temp$preind_prob),
                                      hist_prob = mean(temp$hist_prob),
                                      FAR = mean(temp$FAR),
                                      LCI_FAR = quantile(temp$FAR, 0.025), 
                                      UCI_FAR = quantile(temp$FAR, 0.975)))

}

ggplot(attribution_stats, aes(year, FAR)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin = LCI_FAR,
                      ymax = UCI_FAR), fill = "dark grey", alpha=0.3)

# RR
attribution_stats <- attribution_stats %>%
    mutate(RR = 1/(1-FAR),
           LCI_RR = 1/(1-LCI_FAR),
           UCI_RR = 1/(1-UCI_FAR))


ggplot(attribution_stats, aes(year, RR)) + # no error ribbon b/c UCI is undefined if FAR = 1
    geom_point() +
    geom_line() 


# hist prob
ggplot(attribution_stats, aes(year, hist_prob)) + # no error ribbon b/c UCI is undefined if FAR = 1
    geom_point() +
    geom_line() 


# preind prob
ggplot(attribution_stats, aes(year, preind_prob)) + # no error ribbon b/c UCI is undefined if FAR = 1
    geom_point() +
    geom_line() 


ggplot(attribution_stats, aes(year, preind_prob/hist_prob)) + # no error ribbon b/c UCI is undefined if FAR = 1
    geom_point() +
    geom_line() 

# save 

write.csv(attribution_stats, "./output/probabilistic_attribution_stats.csv", row.names = F)
 
