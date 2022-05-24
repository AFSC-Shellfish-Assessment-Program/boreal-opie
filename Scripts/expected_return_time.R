# Expected return time for SST anomalies associated with negative
# predicted abundance change at high values of borealization index

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

# load borealization DFA trend 

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate) %>%
  filter(year <= 2020) # remove 2021 as we have not yet observed the snow crab response

tail(trend) # 2020 is the value we want to predict for

# load ersst values
ersst <- read.csv("./Data/ersst_time_series.csv")

# load ersst anomalies
ersst.anom <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv")

unique(ersst.anom$region)

ersst.anom <- ersst.anom %>%
  dplyr::filter(region == "Eastern_Bering_Sea")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:1999), aes(annual.anomaly.unsmoothed)) +
  geom_density(fill = "grey")

# load predicted sst - borealization relationship
pred.sst.boreal <- read.csv("./output/sst_borealization_predicted_relationship.csv")

# find sst anomaly predicted to correspond with 2020 borealization trend value
ersst.max <- pred.sst.boreal$sst.anomaly[which.min(abs(pred.sst.boreal$estimate__ - trend$trend[trend$year == 2020]))]

# load CMIP6 anomalies
cmip.anom <- read.csv("./Data/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./Data/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)


# create df to catch outcomes for extreme runs
extreme.outcomes <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
    # separate model and region of interest
    pre.temp <- cmip.anom %>% 
      filter(experiment == "piControl",
             model == models[i],
             region == "Eastern_Bering_Sea")
    
 
    
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                              data.frame(model = models[i],
                                         period = "preindustrial",
                                         count = sum(pre.temp$annual.unsmoothed >= ersst.max),
                                         N = length(pre.temp$annual.unsmoothed)))
    
    
 
  
} # close i loop (models)


head(extreme.outcomes) 

## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
    # separate model and region of interest
    hist.temp <- cmip.anom %>% 
      filter(experiment == "hist_ssp585",
             model == models[i],
             region == "Eastern_Bering_Sea")
    
    # # separate this region from ersst.max
    # ersst.temp <- ersst.max %>%
    #   filter(region == regions[j])
    
    ## pull 1950 - 0.5 degrees warming
    
    use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
    
    # and limit hist.temp to these years
    hist.temp.use <- hist.temp %>%
      filter(year %in% use)
    
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                              data.frame(model = models[i],
                                         period = "1950_to_0.5",
                                         count = sum(hist.temp.use$annual.unsmoothed >= ersst.max),
                                         N = length(hist.temp.use$annual.unsmoothed)))
    
    ## pull 0.5 - 1.0 degrees warming
    
    use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
    
    # and limit hist.temp to these years
    hist.temp.use <- hist.temp %>%
      filter(year %in% use)
    
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                              data.frame(model = models[i],
                                         period = "0.5_to_1.0",
                                         count = sum(hist.temp.use$annual.unsmoothed >= ersst.max),
                                         N = length(hist.temp.use$annual.unsmoothed)))
    
    ## pull 1.0 - 1.5 degrees warming
    
    use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
    
    # and limit hist.temp to these years
    hist.temp.use <- hist.temp %>%
      filter(year %in% use)
    
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                              data.frame(model = models[i],
                                         period = "1.0_to_1.5",
                                         count = sum(hist.temp.use$annual.unsmoothed >= ersst.max),
                                         N = length(hist.temp.use$annual.unsmoothed)))
    
    
    ## pull 1.5 - 2.0 degrees warming
    
    use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
    
    # and limit hist.temp to these years
    hist.temp.use <- hist.temp %>%
      filter(year %in% use)
    
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                              data.frame(model = models[i],
                                         period = "1.5_to_2.0",
                                         count = sum(hist.temp.use$annual.unsmoothed >= ersst.max),
                                         N = length(hist.temp.use$annual.unsmoothed)))
    

  
} # close i loop (models)


# check
check <- extreme.outcomes %>%
  group_by(period) %>%
  summarise(count = sum(count),
            N = sum(N)) %>%
  mutate(prop = count/N)

View(check)

## fit brms model to estimate probabilities-----------

# model weights for extreme events in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load CMIP6 model weights
model.weights <- read.csv("./Data/CMIP6_model_weights_by_region_window.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(window == "annual",
         region == "Eastern_Bering_Sea") %>%
  select(model, scaled.total.weight) 

# calculate EBS-specific model warming weights (based on prediction of experienced warming)

ersst <- read.csv("./Data/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./Data/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
data <- models %>% 
  filter(experiment == "hist_ssp585",
         region == "Eastern_Bering_Sea",
         year %in% 1850:2021) %>% # note that for regional warming we will calculate anomalies wrt 1950-1999 (beginning of trustworthy ERSST)
  select(year, annual.unsmoothed, model)

data <- rbind(data, ersst) 

# calculate 1950:1999 climatology for each model and ersst
climatology <- data %>%
  filter(year %in% 1850:1949) %>%
  group_by(model) %>%
  summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))

# combine climatology and data, calculate anomalies
data <- left_join(data, climatology) %>%
  mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)

# and pivot longer (ersst vs models)
ersst <- data %>%
  filter(model == "ersst") %>%
  select(year, anomaly) %>%
  rename(ersst.anomaly = anomaly)

data <- data %>%
  filter(model != "ersst") %>%
  left_join(., ersst)

# loop through and fit linear ersst - model regressions to get weights
regional_warming_weights <- data.frame()

models <- unique(data$model)


  for(m in 1:length(models)){ # loop through models
    # m <- 1
    
    temp.dat <- data %>%
      filter(model == models[m],
             year %in% 1972:2021)
    
    
    mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
    
    regional_warming_weights <- rbind(regional_warming_weights,
                                      data.frame(model = models[m],
                                                 regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
  }
  



weights <- left_join(model.weights, regional_warming_weights) %>%
  mutate(total_weight = scaled.total.weight * regional_warming_weight)


# plot to examine
ggplot(weights, aes(scaled.total.weight, regional_warming_weight)) +
  geom_point() 

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) 

extremes <- left_join(extreme.outcomes, weights) %>%
  mutate(model_fac = as.factor(model))

## brms: setup ---------------------------------------------

form <-  bf(count | trials(N) + weights(total_weight, scale = TRUE) ~
              period + (1 | model_fac))

# fit model

  extremes_brms <- brm(form,
                       data = extremes,
                       family = binomial(link = "logit"),
                       seed = 1234,
                       cores = 4, chains = 4, iter = 12000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.9, max_treedepth = 14))
  
  saveRDS(extremes_brms, "./output/extremes_binomial.rds")
  



# evaluate 

model <- readRDS("./output/extremes_binomial.rds")

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model) 
trace_plot(model$fit)

# and plot
new.dat <- data.frame(period = unique(extremes$period),
                      model = NA,
                      N = 1000) 

plot.dat <- data.frame()

  probs <- posterior_epred(model, newdata = new.dat, re_formula = NA)/1000 # dive by N to get probability
  
  plot.dat <- rbind(plot.dat,
                    data.frame(period = new.dat$period,
                               prob = apply(probs, 2, median),
                               lower = apply(probs, 2, quantile, probs = 0.025),
                               upper = apply(probs, 2, quantile, probs = 0.975)))


# calculate inverse to get expected return time
plot.dat[,c(2:4)] <- 1/plot.dat[,c(2:4)]

# and change values above 10^4 to 10^4

change <- plot.dat[,c(2:4)] > 10^4

plot.dat[,c(2:4)][change] <- 10^4

# set periods in order
period.order <- data.frame(period = unique(plot.dat$period),
                           period.order = 1:5)

plot.dat <- left_join(plot.dat, period.order) %>%
  mutate(period = reorder(period, period.order))


ggplot(plot.dat, aes(period, prob)) +
  geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3) +
  geom_point(color = "red", size = 4) +
  scale_y_continuous(breaks=c( 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000),
                     minor_breaks = c(2:9, 
                                      seq(20, 90, by = 10),
                                      seq(200, 900, by = 100),
                                      seq(2000, 9000, by = 1000))) +
  coord_trans(y = "pseudo_log") +
  ylab("Expected return time (years)") + 
  xlab("North Pacific warming") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggsave("./figs/extreme_return_time.png", width = 3, height = 6, units = 'in')        
