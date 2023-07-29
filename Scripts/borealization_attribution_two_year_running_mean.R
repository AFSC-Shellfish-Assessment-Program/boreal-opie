# Attribution of Bering Sea borealization - two-year running mean SST and borealization index

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(zoo)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

mc.cores = parallel::detectCores()

# load borealization DFA trend and sst anomaly time series

trend <- read.csv("./output/dfa_trend.csv") %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate) %>%
  mutate(two.yr.index = zoo::rollmean(trend, 2, align = "right", fill = NA))

sst <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(year, annual.anomaly.unsmoothed, annual.anomaly.two.yr.running.mean) %>%
  rename(annual.sst = annual.anomaly.unsmoothed,
         two.yr.sst = annual.anomaly.two.yr.running.mean) %>%
  mutate(test.2yr = zoo::rollmean(annual.sst, 2, align = "right", fill = NA)) 
# so we want right-aligned running means, attribution study uses left-aligned!

dat <- left_join(trend, sst) %>%
  select(year, test.2yr, two.yr.index) %>% # using right-alinged sst!
  rename(two.yr.sst = test.2yr)

ggplot(dat, aes(two.yr.sst, two.yr.index)) +
  geom_point() +
  geom_smooth(method = "gam", color = "blue", se = F) +
  geom_smooth(method = "lm", color = "red", se = F)

# fit brms model ------------------------------

sst_boreal_2yr_formula <-  bf(two.yr.index ~ s(two.yr.sst))

## Show default priors
get_prior(sst_boreal_2yr_formula, dat)

## fit 
sst_boreal_2yr_brm <- brm(sst_boreal_2yr_formula,
                      data = dat,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.9999, max_treedepth = 10))

saveRDS(sst_boreal_2yr_brm, file = "output/sst_boreal_2yr_running_mean_brm.rds")

# run diagnostics
sst_boreal_2yr_brm <- readRDS("./output/sst_boreal_2yr_running_mean_brm.rds")
check_hmc_diagnostics(sst_boreal_2yr_brm$fit)
neff_lowest(sst_boreal_2yr_brm$fit)
rhat_highest(sst_boreal_2yr_brm$fit)
summary(sst_boreal_2yr_brm)
bayes_R2(sst_boreal_2yr_brm)

# plot
## 95% CI
ce1s_1 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sst_boreal_2yr_brm, effect = "two.yr.sst", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$two.yr.sst
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$two.yr.sst[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$two.yr.sst[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$two.yr.sst[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$two.yr.sst[["lower__"]]
# dat_ce[["rug.anom"]] <- c(jitter(unique(dat$two.yr.sst), amount = 0.01),
# rep(NA, 100-length(unique(dat$two.yr.sst))))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  labs(x = "Two-year running mean SST anomaly", y = "Two-year running mean borealization index") +
  geom_text(data = dat, aes(two.yr.sst, two.yr.index, label = year), size = 3)
# geom_rug](aes(x=rug.anom, y=NULL)) 

print(g1)


### use this brms model as the basis for attribution -----------

# same approach as for annual anomalies!

# load observed borealization events 
# reload DFA model
mod <- readRDS("./output/DFA_model.rds")

# calculate DFA trend (borealization index) and get two-yr running mean
observed.borealization <- data.frame(year=1972:2022,
                                     borealization_index = as.vector(mod$states)) %>%
  mutate(two.yr.index = rollmean(borealization_index, 2, align = "right", fill = NA))


# make CMIP6 historical draws from the brms model

# load CMIP6 anomalies
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea",
         experiment == "hist_ssp585",
         year %in% 1940:2040) %>% # (to capture warming that matches period of interest for attribution)
  select(year, model, annual.unsmoothed) 

# get running mean for each model
cmip.running.mean <- data.frame()

models <- unique(cmip.anom$model)

for(m in 1:length(models)){
  # m <- 1
  
  temp <- cmip.anom %>%
    filter(model == models[m])

  cmip.running.mean <- rbind(cmip.running.mean,
                             data.frame(year = temp$year,
                                        model = temp$model,
                                        two.yr.sst = rollmean(temp$annual.unsmoothed, 2, align = "right", fill = NA)))
  
}

# and drop NAs
cmip.running.mean <- na.omit(cmip.running.mean)

## estimate borealization for each CMIP6 SST anomaly from model posteriors

# reload brms model
sst_boreal_2yr_brm <- readRDS("./output/sst_boreal_2yr_running_mean_brm.rds")

cmip.hist.2yr.borealization <- data.frame()

for(i in 1:nrow(cmip.running.mean)){
  # i <- 1
  temp <- data.frame(model = cmip.running.mean$model[i],
                     year = cmip.running.mean$year[i],
                     borealization_index = posterior_predict(sst_boreal_2yr_brm,
                                                             newdata = cmip.running.mean[i,],
                                                             nsamples = 4))
  cmip.hist.2yr.borealization <- rbind(cmip.hist.2yr.borealization,
                                   temp)
  
}

nrow(cmip.hist.2yr.borealization)

# save
write.csv(cmip.hist.2yr.borealization, "./output/cmip_historical_2yr_borealization_draws.csv", row.names = F)

# load preindustrial CMIP6 anomalies and make draws from posteriors
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea",
         experiment == "piControl") %>% 
  select(year, model, annual.two.yr.running.mean) %>% 
  # we can use these left-aligned running means as year identity does not matter for preindustrial runs
  rename(two.yr.sst = annual.two.yr.running.mean) %>%
  na.omit()


# estimate borealization for each CMIP6 SST anomaly from model posteriors
preindustrial.2yr.borealization <- data.frame()

for(i in 1:nrow(cmip.anom)){
  print(i)
  temp <- data.frame(model = cmip.anom$model[i],
                     borealization_index = posterior_predict(sst_boreal_2yr_brm,
                                                             newdata = cmip.anom[i,],
                                                             nsamples = 1)) #  1 draw each, as we have n = 250 preindustrial observations for each model
  preindustrial.2yr.borealization <- rbind(preindustrial.2yr.borealization,
                                       temp)
  
}

nrow(preindustrial.2yr.borealization) 

# save
write.csv(preindustrial.2yr.borealization,
          "./output/cmip_preindustrial_historical_2yr_borealization_draws.csv", row.names = F)

preindustrial.2yr.borealization <- read.csv("./output/cmip_preindustrial_historical_2yr_borealization_draws.csv")

# load model weights 
weights <- read.csv("./data/normalized_CMIP6_weights.csv")

weights <- weights %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(model, normalized_weight) %>%
  rename(model_weight = normalized_weight)


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

# remove NAs from observed borealization
observed.borealization <- na.omit(observed.borealization)

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
    pre.temp <- preindustrial.2yr.borealization %>% 
      filter(model == models[i])
    
    
    # record preindustrial outcomes
    outcomes <- rbind(outcomes,
                      data.frame(year = obs.temp$year,
                                 observed_borealization = obs.temp$two.yr.index,
                                 model = models[i],
                                 period = "preindustrial",
                                 n = 249,
                                 outcome = if_else(pre.temp$borealization_index >= obs.temp$two.yr.index, 1, 0)))
    
    
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
    
    hist.temp.use <- cmip.hist.2yr.borealization %>%
      filter(model == models[i],
             year %in% use$year)
    
    if(nrow(hist.temp.use) > 0){ # record outcomes only if we have observations in 1940:2040 matching this warming window
      
      outcomes <- rbind(outcomes,
                        data.frame(year = obs.temp$year,
                                   observed_borealization = obs.temp$borealization_index,
                                   model = models[i],
                                   n = nrow(hist.temp.use),
                                   period = "historical",
                                   outcome = if_else(hist.temp.use$borealization_index >= obs.temp$two.yr.index, 1, 0)))
    }
    
  } # close i loop (models) 
  
} # close h loop (observations)


# save
write.csv(outcomes, "./output/borealization_2yr_outcomes_for_attribution.csv", row.names = F)

outcomes <- read.csv("./output/borealization_2yr_outcomes_for_attribution.csv")

# join with model weights
outcomes <- left_join(outcomes, weights) %>%
  mutate(total_weight = model_weight/n)

# resample and calculate FAR!

# create df for output
resampled_2yr_FAR <- data.frame()

# now resample
periods <- unique(outcomes$period)

for(y in 1973:2022){
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
    resampled_2yr_FAR <- rbind(resampled_2yr_FAR,
                           data.frame(year = y,
                                      draw = x,
                                      FAR = 1 - (preind_prob/hist_prob)))
    
    
  }
  
}



# and summarize results
attribution_stats_2yr <- data.frame()

for(y in 1973:2022){
  
  temp <- resampled_2yr_FAR %>%
    filter(year == y)
  
  attribution_stats_2yr <- rbind(attribution_stats_2yr,
                             data.frame(year = y,
                                        FAR = mean(temp$FAR),
                                        LCI_FAR = quantile(temp$FAR, 0.025), 
                                        UCI_FAR = quantile(temp$FAR, 0.975)))
  
}

ggplot(attribution_stats_2yr, aes(year, FAR)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = LCI_FAR,
                  ymax = UCI_FAR), fill = "dark grey", alpha=0.3)

# RR
attribution_stats_2yr <- attribution_stats_2yr %>%
  mutate(RR = 1/(1-FAR),
         LCI_RR = 1/(1-LCI_FAR),
         UCI_RR = 1/(1-UCI_FAR))

ggplot(attribution_stats_2yr, aes(year, RR)) + # no error ribbon b/c UCI is undefined if FAR = 1
  geom_point() +
  geom_line() 


# save 

write.csv(attribution_stats_2yr, "./output/probabilistic_attribution_stats_two_yr_rolling_means.csv", row.names = F)

