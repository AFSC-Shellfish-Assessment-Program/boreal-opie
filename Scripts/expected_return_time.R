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
  rename(year = t, trend = estimate) 

tail(trend) # 2020 is the value we want to predict for

# load ersst values
ersst <- read.csv("./Data/regional_north_pacific_ersst_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

# load ersst anomalies
ersst.anom <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

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

# load CMIP6 model weights
model.weights <- read.csv("./Data/normalized_CMIP6_weights.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(model, normalized_weight) 

extremes <- left_join(extreme.outcomes, model.weights) %>%
  mutate(model_fac = as.factor(model))

## brms: setup ---------------------------------------------

form <-  bf(count | trials(N) + weights(normalized_weight, scale = TRUE) ~
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


return.time <- ggplot(plot.dat, aes(period, prob)) +
  geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3) +
  geom_point(color = "red", size = 4) +
  scale_y_continuous(breaks=c( 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000),
                     labels=c( "1","2","5","10","20","50","100","200","500","1000","2000","5000",">10000"),
                     minor_breaks = c(2:9, 
                                      seq(20, 90, by = 10),
                                      seq(200, 900, by = 100),
                                      seq(2000, 9000, by = 1000))) +
  coord_trans(y = "pseudo_log") +
  ylab("Expected return time (years)") + 
  xlab("North Pacific warming") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggsave("./figs/extreme_return_time_3.575092_SD.png", width = 3, height = 6, units = 'in')    

## alternate approach: plot temperature pdfs for warming eras --------------------------------------

## different approach - plot hindcast and projected pdfs for sst anomalies --------------------------

# load ersst anomalies
ersst.anom <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

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
anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Eastern_Bering_Sea")
  
  
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "preindustrial",
                                   anomaly = na.omit(pre.temp$annual.unsmoothed)))
  
  
  
} # close i loop (models)


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
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  
} # close i loop (models)


# load CMIP6 model weights
model.weights <- read.csv("./Data/normalized_CMIP6_weights.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(model, normalized_weight) 

anomaly.pdfs <- left_join(anomaly.pdfs, model.weights) 

# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1
  
  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 1000, replace = T, prob = temp$normalized_weight)))
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)


resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))


sum <- resample.pdf %>%
  group_by(period) %>%
  summarize(proportion = sum(anomaly > ersst.max) / n())

sum

# and plot

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sst.pdfs <- ggplot(resample.pdf, aes(period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("SST anomaly wrt 1850-1949 (SD)") +
  geom_hline(yintercept = ersst.max, lty = 2) +
  scale_x_discrete(labels = c("Preindustrial",
                              "1950 to 0.5°",
                              "0.5° to 1.0°",
                              "1.0° to 1.5°",
                              "1.5° to 2.0°"))


## finally, plot hindcast / observed / projected warming rates

# load model objects

warming.245 <- readRDS("./Data/inverse_warming_brm_ssp245.rds")

warming.585 <- readRDS("./Data/inverse_warming_brm.rds")

warming.ersst<- readRDS("./Data/inverse_warming_brm_ersst.rds")


# plot

new.dat <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                      model_fac = NA, weight = 1)

pred.245 <- posterior_epred(warming.245, newdata = new.dat)

pred.585 <- posterior_epred(warming.585, newdata = new.dat)

new.dat <- data.frame(ersst.warming = c(0.5, 1.0))

pred.ersst <- posterior_epred(warming.ersst, newdata = new.dat) 


## SST anomaly predictions #### 95% CI

# 245

ce1s_245 <- conditional_effects(warming.245, effect = "warming", re_formula = NA,
                                probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_245$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot.245 <- data.frame(source = "SSP 245",
                            warming = c(0.5, 1.0, 1.5, 2.0),
                            year = ce1s_245$warming$estimate__[choose],
                            UCI = ce1s_245$warming$upper__[choose],
                            LCI = ce1s_245$warming$lower__[choose])


# 585 

ce1s_585 <- conditional_effects(warming.585, effect = "warming", re_formula = NA,
                                probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_585$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot.585 <- data.frame(source = "SSP 585",
                            warming = c(0.5, 1.0, 1.5, 2.0),
                            year = ce1s_585$warming$estimate__[choose],
                            UCI = ce1s_585$warming$upper__[choose],
                            LCI = ce1s_585$warming$lower__[choose])


# ersst
ce1s_ERSST <- conditional_effects(warming.ersst, effect = "ersst.warming", re_formula = NA,
                                  probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_ERSST$ersst.warming$ersst.warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)))

pred.plot.ERSST <- data.frame(source = "ERSST",
                              warming = c(0.5, 1.0),
                              year = ce1s_ERSST$ersst.warming$estimate__[choose],
                              UCI = ce1s_ERSST$ersst.warming$upper__[choose],
                              LCI = ce1s_ERSST$ersst.warming$lower__[choose])

pred.plot <- rbind(pred.plot.245,
                   pred.plot.585,
                   pred.plot.ERSST)
pred.plot$warming <- as.factor(pred.plot$warming)

pos_dodge = position_dodge(width = 0.2)


time.plot <- ggplot(pred.plot, aes(warming, year, color = source)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2, position = pos_dodge) +
  geom_point(size = 4, position = pos_dodge) +
  labs(x = "North Pacific warming (°C)",
       y = "Year reached") +
  scale_color_manual(values = cb[c(2,6,7)]) +
  theme(legend.position = c(0.75, 0.15),
        legend.title = element_blank())


# and combine the three plots

# get blank plot
blank <- ggplot + theme_void()


png("./Figs/extremes_probability_warming_timing.png", width = 10, height = 5, units = 'in', res = 300) 

ggpubr::ggarrange(time.plot,
                  sst.pdfs,
                  return.time,
                  labels = "auto",
                  ncol = 3,
                  widths = c(0.3, 0.4, 0.3))

dev.off()

# and a 2-panel version
png("./Figs/extremes_probability_warming_timing_2_panel.png", width = 7, height = 5, units = 'in', res = 300) 

ggpubr::ggarrange(sst.pdfs,
                  time.plot,
                  labels = "auto",
                  ncol = 2,
                  widths = c(0.4, 0.3))

dev.off()
