## estimate preindustrial and historical 
## probabilities for given borealization value

library(tidyverse)


theme_set(theme_bw())

# load observed borealization events 
# drawn from SST-borealization posteriors applied to observed SST values
observed.borealization <- read.csv("./output/observed_borealization_from_posteriors.csv")

# load cmip6 historical draws from the same model
historical.borealization <- read.csv("./output/cmip_historical_borealization_draws.csv")

# load preindustrial CMIP6 anomalies and make draws from posteriors
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
    filter(region == "Eastern_Bering_Sea",
           experiment == "piControl") %>% 
    select(year, model, annual.unsmoothed) %>%
    rename(sst.anomaly = annual.unsmoothed)

# estimate borealization for each CMIP6 SST anomaly from model posteriors
preindustrial.borealization <- data.frame()

for(i in 1:nrow(cmip.anom)){
    
    temp <- data.frame(model = cmip.anom$model[i],
                       borealization_index = posterior_predict(sst_boreal_brm,
                                                               newdata = cmip.anom[i,],
                                                               ndraws = 1)) #  1 draw each, as we have n = 250 preindustrial observations for each model
    preindustrial.borealization <- rbind(preindustrial.borealization,
                                       temp)
    
}

nrow(preindustrial.borealization) 


# load model weights 
weights <- read.csv("./data/normalized_CMIP6_weights.csv")

weights <- weights %>%
    filter(region == "Eastern_Bering_Sea") %>%
   select(model, normalized_weight) %>%
   rename(model_weight = normalized_weight)


# load brms estimate of warming trend for each model
model.warming.trends <- read.csv("./data/CMIP6_brms_warming_rate.csv")

# check model names
check <- data.frame(models = models,
                    trend.names = unique(model.warming.trends$model))
check # line up fine

# and load predicted warming rate across all models 
predicted.warming <- read.csv("./data/brms_predicted_North_Pac_warming.csv")

###################
# load estimated warming level timing for each model
timing <- read.csv("./Data/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch preindustrial and historical probabilities for each borealization observation
probs <- data.frame()

# loop through each observation
for(h in 1:nrow(observed.borealization)){
    # h <- 1

    obs.temp <- observed.borealization[h,]
    
    # define 15-year window
    window <- (obs.temp$year[k] - 7) : (obs.temp$year[k] + 7)
    
    # define range of predicted warming values for this window
    warming.range <- range(predicted.warming$pred_mean[predicted.warming$year %in% window])
    
    
    # loop through each model
for(i in 1:length(models)){ # start i loop (models)
    # i <- 1
    
    # separate model of interest from preindustrial runs
    pre.temp <- preindustrial.borealization %>% 
        filter(model == models[i])
    
    
    # record preindustrial probability
    probs <- rbind(probs,
                                data.frame(year = obs.temp$year,
                                           observed_borealization = obs.temp$borealization_index,
                                            model = models[i],
                                           period = "preindustrial",
                                           n = 250,
                                           probabililty = sum(pre.temp$borealization_index >= obs.temp$borealization_index)/250))
    
    
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
    
    probs <- rbind(probs,
                   data.frame(year = obs.temp$year,
                              observed_borealization = obs.temp$borealization_index,
                              model = models[i],
                              period = "historical",
                              n = nrow(hist.temp.use),
                              probabililty = sum(hist.temp.use$borealization_index >= obs.temp$borealization_index)/nrow(hist.temp.use)))
    
} # close i loop (models) 
    
} # close h loop (observations)
    


# save
write.csv(probs, "./output/borealization_probabilities.csv", row.names = F)

probs <- read.csv("./output/borealization_probabilities.csv")


# summarize!

# create df for output
summarized_probabilities <- data.frame()

# remove n = 0 cases and calcuate # of cases of greater borealization
probs <- probs %>%
    filter(n > 0) %>%
    mutate(positive_cases = n * probabililty)

summarized_probabilities <- probs %>%
    group_by(year, period, model) %>%
    summarise(n = sum(n), 
              positive_cases = sum(positive_cases)) %>%
    mutate(probabililty = positive_cases / n)


summarized_probabilities <- left_join(summarized_probabilities, weights) %>%
    mutate(total_weight = n*model_weight)

# now resample
periods <- unique(summarized_probabilities$period)
resampled_probs <- data.frame()

for(p in 1:2){
    # p <- 1
    for(y in 1950:2022){
     # y <- 1950
    
     temp <- summarized_probabilities %>%
        filter(period == periods[p],
               year == y) 
        
     resampled_probs <- rbind(resampled_probs,
                              data.frame(year = y,
                                         period = periods[p],
                                         probability = weighted.mean(temp$probabililty, w = temp$total_weight)))   
        
    }
    
}


wide_probs <- resampled_probs %>%
    pivot_wider(names_from = period, 
                values_from = probability) %>%
    mutate(FAR = 1 - (preindustrial/historical),
           RR = 1/(1-FAR))

# save 

write.csv(wide_probs, "./output/probabilistic_attribution_stats.csv", row.names = F)

ggplot(wide_probs, aes(year, FAR)) +
    geom_point() +
    geom_line()

ggplot(wide_probs, aes(year, RR)) +
    geom_point() +
    geom_line()



wide_probs$FAR <- 1 - (preindustrial)

resample.pdf <- data.frame()

periods <- unique(borealization.pdfs$period)

for(i in 1:length(periods)){
    # i <- 1
    
    temp <- borealization.pdfs[borealization.pdfs$period == periods[i],]
    
    resample.pdf <- rbind(resample.pdf,
                          data.frame(period = periods[i],
                                     borealization_index = sample(temp$borealization_index, 1000, replace = T, prob = temp$normalized_weight)))
    
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)


resample.pdf <- left_join(resample.pdf, plot.order) %>%
    mutate(period =  reorder(period, order))


sum <- resample.pdf %>%
    group_by(period) %>%
    summarize(proportion = sum(borealization_index > 2) / n()) # using 1.5 as a critical value

sum

# and plot
theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

boreal.pdf.plot <- ggplot(resample.pdf, aes(period, borealization_index)) +
    geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
    coord_flip(ylim = c(-4.5, 4.5)) +
    xlab("North Pacific warming") +
    ylab("Borealization index") +
    geom_hline(yintercept = 2, lty = 2) +
    scale_x_discrete(labels = c("Preindustrial",
                                "1950 to 0.5°",
                                "0.5° to 1.0°  
                              (~2003 - 2019)",
                                "1.0° to 1.5°  ,
                              (~2019 - 2038)",
                                "1.5° to 2.0°  
                              (~2038 - 2060)")) +
    scale_y_continuous(breaks = seq(-4,4, by = 2))


boreal.pdf.plot


###################
## Calc FAR
far <- 1 - (pre_prob / his_prob)
range(far, na.rm = TRUE)


far_pred <- data.frame(annual.anomaly.3yr = nd_pre$annual.anomaly.3yr,
                       prob = apply(far, 2, mean),
                       lower = apply(far, 2, quantile, probs = 0.025),
                       upper = apply(far, 2, quantile, probs = 0.975))

g <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = annual.anomaly.3yr, y = prob), size = 0.8) +
    geom_ribbon(aes(x = annual.anomaly.3yr, ymin = lower, ymax = upper), alpha = 0.15)
print(g)
