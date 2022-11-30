# compare estimated probability distributions
# between observation time blocks used in rebuilding projections
# and CMIP6 projections

library(tidyverse)

theme_set(theme_bw())

obs <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv")

head(obs)

# process
obs <- obs %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(anomaly = annual.anomaly.unsmoothed)

# plot to check
block1 <- obs %>%
  filter(year %in% 1982:2017) %>%
  mutate(time_block = "1982-2017")

block2 <- obs %>%
  filter(year %in% 1982:2019) %>%
  mutate(time_block = "1982-2019")

block3 <- obs %>%
  filter(year %in% 1994:2019) %>%
  mutate(time_block = "1994-2019")

block4 <- obs %>%
  filter(year %in% 2005:2019) %>%
  mutate(time_block = "2005-2019")

plot.obs <- rbind(block1, block2, block3, block4)

ggplot(plot.obs, aes(anomaly)) +
  geom_histogram(bins = 7, fill = "grey", color = "black") +
  facet_wrap(~time_block, scales = "free_y") +
  xlab("Anomaly with respect to 1854-1949")

ggsave("./Figs/rebuilding_time_block_sst_histograms.png", width = 6, height = 5, units = 'in')

# plot time series for each block

# first, for full observation time series
ggplot(obs, aes(year, anomaly)) +
  geom_smooth(se = F) + 
  geom_line() +
  geom_hline(yintercept = c(0, 4), lty = c(1, 2)) +
  geom_point() +
  ylab("Anomaly with respect to 1854-1949 (Std. dev.)") +
  scale_x_continuous(breaks = seq(1850, 2025, 25)) +
  theme(axis.title.x = element_blank())

ggsave("./Figs/observation_time_series_1854-2022.png", width = 6, height = 4, units = 'in')

# summarize proportion of years > 4 SD for each block
prop_4sd <- plot.obs %>%
  group_by(time_block) %>%
  summarize(percentage_4sd = 100*round(sum(anomaly > 4)/n(),2))

plot.obs <- left_join(plot.obs, prop_4sd) %>%
  mutate(time_block_percentage = paste(time_block, " (", percentage_4sd, "% observations > 4sd)", sep = ""))

ggplot(plot.obs, aes(year, anomaly)) +
  geom_line() +
  geom_point() +
  facet_wrap(~time_block_percentage, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 4, lty = 2) +
  ylab("Anomaly with respect to 1854-1949 (Std. dev.)") +
  theme(axis.title.x = element_blank())

ggsave("./Figs/rebuilding_block_sst_time_series.png", width = 6, height = 4, units = 'in')

# # resample to try to define underlying 
# blocks <- unique(plot.obs$time_block)
# 
# resample_obs <- data.frame()
# 
# for(b in 1:length(blocks)){
#   
#   temp <- plot.obs %>%
#     filter(time_block == blocks[b]) 
#   
#   resample_obs <- rbind(resample_obs,
#                         data.frame(time_block = blocks[b],
#                                    anomaly = sample(temp$anomaly, 1000, replace = T)))
#   
# }
# 
# ggplot(resample_obs, aes(time_block, anomaly)) +
#   geom_violin(alpha = 0.2) 

### CMIP6

# load CMIP6 anomalies
cmip.anom <- read.csv("./Data/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./Data/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# # create df to catch outcomes for extreme runs
# extreme.outcomes <- data.frame()
# 
# # loop through each model
# for(i in 1:length(models)){ # start i loop (models)
#   # i <- 1
#   
#   # separate model and region of interest
#   pre.temp <- cmip.anom %>% 
#     filter(experiment == "piControl",
#            model == models[i],
#            region == "Gulf_of_Alaska")
#   
#   
#   
#   # record how many model years are more extreme
#   extreme.outcomes <- rbind(extreme.outcomes,
#                             data.frame(model = models[i],
#                                        period = "preindustrial",
#                                        count = sum(pre.temp$annual.three.yr.running.mean >= ersst.max, na.rm = T),
#                                        N = length(na.omit(pre.temp$annual.three.yr.running.mean))))
#   
#   
#   
#   
# } # close i loop (models)
# 
# 
# head(extreme.outcomes) 
# 
# ## record outcomes using different warming levels from hist.585 ----------------
# 
# 
# # loop through each model
# for(i in 1:length(models)){ # start i loop (models)
#   # i <- 1
#   
#   # separate model and region of interest
#   hist.temp <- cmip.anom %>% 
#     filter(experiment == "hist_ssp585",
#            model == models[i],
#            region == "Gulf_of_Alaska")
#   
#   # # separate this region from ersst.max
#   # ersst.temp <- ersst.max %>%
#   #   filter(region == regions[j])
#   
#   ## pull 1950 - 0.5 degrees warming
#   
#   use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
#   
#   # and limit hist.temp to these years
#   hist.temp.use <- hist.temp %>%
#     filter(year %in% use)
#   
#   # record how many model years are more extreme
#   extreme.outcomes <- rbind(extreme.outcomes,
#                             data.frame(model = models[i],
#                                        period = "1950_to_0.5",
#                                        count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
#                                        N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))) )
#   
#   ## pull 0.5 - 1.0 degrees warming
#   
#   use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
#   
#   # and limit hist.temp to these years
#   hist.temp.use <- hist.temp %>%
#     filter(year %in% use)
#   
#   # record how many model years are more extreme
#   extreme.outcomes <- rbind(extreme.outcomes,
#                             data.frame(model = models[i],
#                                        period = "0.5_to_1.0",
#                                        count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
#                                        N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
#   
#   ## pull 1.0 - 1.5 degrees warming
#   
#   use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
#   
#   # and limit hist.temp to these years
#   hist.temp.use <- hist.temp %>%
#     filter(year %in% use)
#   
#   # record how many model years are more extreme
#   extreme.outcomes <- rbind(extreme.outcomes,
#                             data.frame(model = models[i],
#                                        period = "1.0_to_1.5",
#                                        count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
#                                        N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
#   
#   
#   ## pull 1.5 - 2.0 degrees warming
#   
#   use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
#   
#   # and limit hist.temp to these years
#   hist.temp.use <- hist.temp %>%
#     filter(year %in% use)
#   
#   # record how many model years are more extreme
#   extreme.outcomes <- rbind(extreme.outcomes,
#                             data.frame(model = models[i],
#                                        period = "1.5_to_2.0",
#                                        count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
#                                        N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
#   
#   
#   
# } # close i loop (models)
# 
# 
# # check
# check <- extreme.outcomes %>%
#   group_by(period) %>%
#   summarise(count = sum(count),
#             N = sum(N)) %>%
#   mutate(prop = count/N)
# 
# View(check)
# 
# ## fit brms model to estimate probabilities-----------
# 
# # model weights for extreme events in different periods - 
# # product of regional weighting (based on ar(1), correlation, bias) and 
# # prediction of observed N. Pac. weighting
# 
# # load CMIP6 model weights
# model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 
# 
# # clean up model weights 
# model.weights <- model.weights %>%
#   filter(window == "annual",
#          region == "Gulf_of_Alaska") %>%
#   select(model, scaled.total.weight) 
# 
# # calculate GOA-specific model warming weights (based on prediction of experienced warming)
# 
# ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")
# 
# ersst <- ersst %>%
#   select(year, annual.unsmoothed) %>%
#   mutate(model = "ersst")
# 
# models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")
# 
# # combine models and ersst observations into "data"
# data <- models %>% 
#   filter(experiment == "hist_ssp585",
#          region == "Gulf_of_Alaska",
#          year %in% 1850:2021) %>% 
#   select(year, annual.unsmoothed, model)
# 
# data <- rbind(data, ersst) 
# 
# # calculate 1850:1949 climatology for each model and ersst
# climatology <- data %>%
#   filter(year %in% 1850:1949) %>%
#   group_by(model) %>%
#   summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))
# 
# # combine climatology and data, calculate anomalies
# data <- left_join(data, climatology) %>%
#   mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)
# 
# # and pivot longer (ersst vs models)
# ersst <- data %>%
#   filter(model == "ersst") %>%
#   select(year, anomaly) %>%
#   rename(ersst.anomaly = anomaly)
# 
# data <- data %>%
#   filter(model != "ersst") %>%
#   left_join(., ersst)
# 
# # loop through and fit linear ersst - model regressions to get weights
# regional_warming_weights <- data.frame()
# 
# models <- unique(data$model)
# 
# 
# for(m in 1:length(models)){ # loop through models
#   # m <- 1
#   
#   temp.dat <- data %>%
#     filter(model == models[m],
#            year %in% 1972:2021)
#   
#   
#   mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
#   
#   regional_warming_weights <- rbind(regional_warming_weights,
#                                     data.frame(model = models[m],
#                                                regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
# }
# 
# 
# 
# 
# weights <- left_join(model.weights, regional_warming_weights) %>%
#   mutate(total_weight = scaled.total.weight * regional_warming_weight)
# 
# 
# # plot to examine
# ggplot(weights, aes(scaled.total.weight, regional_warming_weight)) +
#   geom_point() 
# 
# ggplot(weights, aes(total_weight)) +
#   geom_histogram(fill = "grey", color = "black", bins = 20) 
# 
# extremes <- left_join(extreme.outcomes, weights) %>%
#   mutate(model_fac = as.factor(model))
# 
# ## brms: setup ---------------------------------------------
# 
# form <-  bf(count | trials(N) + weights(total_weight, scale = TRUE) ~
#               period + (1 | model_fac))
# 
# # fit model
# 
# extremes_brms <- brm(form,
#                      data = extremes,
#                      family = binomial(link = "logit"),
#                      seed = 1234,
#                      cores = 4, chains = 4, iter = 12000,
#                      save_pars = save_pars(all = TRUE),
#                      control = list(adapt_delta = 0.9, max_treedepth = 14))
# 
# saveRDS(extremes_brms, "./CMIP6/brms_output/extremes_binomial.rds")
# 
# 
# # evaluate 
# 
# model <- readRDS("./CMIP6/brms_output/extremes_binomial.rds")
# 
# check_hmc_diagnostics(model$fit)
# neff_lowest(model$fit) 
# rhat_highest(model$fit)
# summary(model)
# bayes_R2(model) 
# trace_plot(model$fit)
# 
# # and plot
# new.dat <- data.frame(period = unique(extremes$period),
#                       model = NA,
#                       N = 1000) 
# 
# plot.dat <- data.frame()
# 
# probs <- posterior_epred(model, newdata = new.dat, re_formula = NA)/1000 # dive by N to get probability
# 
# plot.dat <- rbind(plot.dat,
#                   data.frame(period = new.dat$period,
#                              prob = apply(probs, 2, median),
#                              lower = apply(probs, 2, quantile, probs = 0.025),
#                              upper = apply(probs, 2, quantile, probs = 0.975)))
# 
# 
# # calculate inverse to get expected return time
# plot.dat[,c(2:4)] <- 1/plot.dat[,c(2:4)]
# 
# # and change values above 10^4 to 10^4
# 
# change <- plot.dat[,c(2:4)] > 10^4
# 
# plot.dat[,c(2:4)][change] <- 10^4
# 
# # set periods in order
# period.order <- data.frame(period = unique(plot.dat$period),
#                            period.order = 1:5)
# 
# plot.dat <- left_join(plot.dat, period.order) %>%
#   mutate(period = reorder(period, period.order))
# 
# 
# ggplot(plot.dat, aes(period, prob)) +
#   geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3) +
#   geom_point(color = "red", size = 4) +
#   scale_y_continuous(breaks=c( 1,2,5,10,20,50,100,200,500,1000,2000,5000),
#                      minor_breaks = c(2:9, 
#                                       seq(20, 90, by = 10),
#                                       seq(200, 900, by = 100),
#                                       seq(2000, 9000, by = 1000))) +
#   coord_trans(y = "pseudo_log") +
#   ylab("Expected return time (years)") + 
#   xlab("North Pacific warming") +
#   theme(axis.text.x = element_text(angle = 45,
#                                    hjust = 1))
# 
# ggsave("./CMIP6/figs/sockeye_extreme_return_time.png", width = 3, height = 6, units = 'in')        

## different approach - plot hindcast and projected pdfs for sst anomalies --------------------------
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
                                   anomaly = na.omit(pre.temp$annual.three.yr.running.mean)))
  
  
  
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
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  
} # close i loop (models)


# model weights for anomalies in different periods - 
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

# calculate 1850:1949 climatology for each model and ersst
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

anomaly.pdfs <- left_join(anomaly.pdfs, weights) 


# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1
  
  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 1000, replace = T, prob = temp$total_weight)))
  
  
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)



resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))

# and plot

# get good labels for each warming period
labs <- data.frame(period = unique(resample.pdf$period),
                   plot_period = c("Preindustrial",
                                   "1950 to 0.5°",
                                   "0.5° to 1.0°",
                                   "1.0° to 1.5°",
                                   "1.5° to 2.0°"),
                   plot_order = 1:5)

resample.pdf <- left_join(resample.pdf, labs) %>%
  mutate(plot_period = reorder(plot_period, plot_order))


cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf_plot <- ggplot(resample.pdf, aes(plot_period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("SST anomaly (Std. Dev.)") +
  geom_hline(yintercept = ersst.max, lty = 2) +
  labs(tag = "D")

pdf_plot

# and get summary statistics (proportion above ersst.max threshold)

summary_pdf <- resample.pdf %>%
  group_by(period) %>%
  summarise(proportion = sum(anomaly >= ersst.max) / length(anomaly))

summary_pdf


ggsave("./CMIP6/figs/sockeye_anomaly_pdfs.png", width = 3, height = 3, units = 'in')        

# combine plots

void_plot <- ggplot() + theme_void()

png("./CMIP6/figs/combined_sockeye_FAR_plot.png", width = 9, height = 7, units = 'in', res = 300)

ggpubr::ggarrange(ggpubr::ggarrange(far_plot,  catch_plot, ggpubr::ggarrange(void_plot, g3, ncol = 2, widths = c(0.3, 0.7)), ncol = 1),
                  pdf_plot, ncol = 2, widths = c(0.45, 0.55))

dev.off()

