# compare borealization models with bottom temp models;
# evaluate effects of imputed 2020 abundance as a covariate

library(tidyverse)
library(mgcv)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# load borealization DFA trend
trend <- read.csv("./output/dfa_trend.csv")

# add rolling means and lags to figure borealization effects over differ time scales
trend <- trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, trend = estimate) %>%
  mutate(trend2 = zoo::rollmean(trend, 2, fill = NA, align = "right"),
         trend3 = zoo::rollmean(trend, 3, fill = NA, align= "right"),
         trend3_lag1 = lag(trend3, 1),
         trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1),
         trend2_lag2 = lag(trend2, 2))

# load imputed bottom temp time series
temp <- read.csv("./output/imputed_bottom_temp.csv")

# add rolling means and lags for bottom temp
temp <- temp %>%
  dplyr::select(year, Bottom_temp) %>%
  rename(temp = Bottom_temp) %>%
  mutate(temp2 = zoo::rollmean(temp, 2, fill = NA, align = "right"),
         temp3 = zoo::rollmean(temp, 3, fill = NA, align= "right"),
         temp3_lag1 = lag(temp3, 1),
         temp_lag1 = lag(temp, 1),
         temp2_lag1 = lag(temp2, 1),
         temp2_lag2 = lag(temp2, 2))

# and join
dat <- left_join(trend, temp)

## analyze male abundance -------------------

abundance_male <- read.csv("./output/male_imputed_weighted_log_cpue.csv")

male_dat <- left_join(dat, abundance_male)

# fill in log_mean_lag1 for 2021
male_dat$log_mean_lag1[male_dat$year == 2021] <- mean(c(male_dat$log_mean_lag1[male_dat$year == 2020], male_dat$log_mean_lag1[male_dat$year == 2022]))

# fit models with dfa trend

male_boreal_mod1 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod1)

male_boreal_mod2 <- gam(log_mean ~ log_mean_lag1 + s(trend2) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod2)

male_boreal_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend3) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod3)

male_boreal_mod4 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod4)
# plot(male_boreal_mod4, resid = T, pch = 19, se = T)

male_boreal_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod5)

male_boreal_mod6 <- gam(log_mean ~  log_mean_lag1 + s(trend_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod6)

male_boreal_mod7 <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = male_dat)
summary(male_boreal_mod7)
# plot(male_boreal_mod7, resid = T, pch = 19, se = T)

boreal.aicc <- MuMIn::AICc(male_boreal_mod1, male_boreal_mod2, male_boreal_mod3, male_boreal_mod4, 
                           male_boreal_mod5, male_boreal_mod6, male_boreal_mod7)

# same model structures with bottom temperature as the covariate

male_temp_mod1 <- gam(log_mean ~  log_mean_lag1 + s(temp2_lag1) + s(year, k = 4), data = male_dat)
summary(male_temp_mod1)

male_temp_mod2 <- gam(log_mean ~ log_mean_lag1 + s(temp2) + s(year, k = 4), data = male_dat)
summary(male_temp_mod2)

male_temp_mod3 <- gam(log_mean ~ log_mean_lag1 + s(temp3) + s(year, k = 4), data = male_dat)
summary(male_temp_mod3)

male_temp_mod4 <- gam(log_mean ~  log_mean_lag1 + s(temp2_lag2) + s(year, k = 4), data = male_dat)
summary(male_temp_mod4)

male_temp_mod5 <- gam(log_mean ~  log_mean_lag1 + s(temp) + s(year, k = 4), data = male_dat)
summary(male_temp_mod5)

male_temp_mod6 <- gam(log_mean ~  log_mean_lag1 + s(temp_lag1) + s(year, k = 4), data = male_dat)
summary(male_temp_mod6)

male_temp_mod7 <- gam(log_mean ~  log_mean_lag1 + s(temp3_lag1) + s(year, k = 4), data = male_dat)
summary(male_temp_mod7)

temp.aicc <- MuMIn::AICc(male_temp_mod1, male_temp_mod2, male_temp_mod3, male_temp_mod4, male_temp_mod5, 
                         male_temp_mod6, male_temp_mod7)

# get minimum AICc value for the different models
min.aicc <- min(boreal.aicc$AICc, temp.aicc$AICc)

# combine to plot
plot.male.aicc <- rbind(boreal.aicc, temp.aicc) %>%
  mutate(`Covariate` = rep(c("Borealization index", "Bottom temperature"), each = 7),
         `Delta AICc` = AICc - min.aicc,
         plot.lags = rep(c("Year -1, -2",
                        "Year 0, -1",
                        "Year 0, -1, -2",
                        "Year -2, -3",
                        "Year 0",
                        "Year -1",
                        "Year -1, -2, -3"), 2))
         

plot.male.aicc <- plot.male.aicc %>%
  arrange(`Covariate`, desc(`Delta AICc`))
         
# add plot order (order of delta-AICc for boreal models)
plot.male.aicc$order <- order(plot.male.aicc$AICc[plot.male.aicc$Covariate == "Borealization index"])

plot.labels <- reorder(plot.male.aicc$plot.lags[1:7], plot.male.aicc$order[1:7])

plot.male.aicc$covariate.order <- if_else(plot.male.aicc$Covariate == "Borealization index", 2, 1)

plot.male.aicc$Covariate <- reorder(plot.male.aicc$Covariate, plot.male.aicc$covariate.order)

male.plot <- ggplot(plot.male.aicc, aes(order, `Delta AICc`, color = `Covariate`)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:7,
                     labels = plot.male.aicc$plot.lags[7:1]) +
  scale_color_manual(values = cb[c(8,6)]) +
  labs(x = "Covariate lag", y = expression(Delta~"AICc")) +
  theme(axis.text.x  = element_text(angle=60, hjust=1)) +
  ggtitle("b) Male")

male.plot

## analyze female abundance -------------------

abundance_female <- read.csv("./output/female_drop5_df_simple.csv", row.names = 1) %>%
  rename(log_mean = mean) %>%
  arrange(year) %>%
  mutate(log_mean_lag1 = lag(log_mean, 1))

female_dat <- left_join(dat, abundance_female)

# fill in log_mean_lag1 for 2021
female_dat$log_mean_lag1[female_dat$year == 2021] <- mean(c(female_dat$log_mean_lag1[female_dat$year == 2020], female_dat$log_mean_lag1[female_dat$year == 2022]))


# fit models with dfa trend

female_boreal_mod1 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod1)

female_boreal_mod2 <- gam(log_mean ~ log_mean_lag1 + s(trend2) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod2)

female_boreal_mod3 <- gam(log_mean ~ log_mean_lag1 + s(trend3) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod3)

female_boreal_mod4 <- gam(log_mean ~  log_mean_lag1 + s(trend2_lag2) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod4)

female_boreal_mod5 <- gam(log_mean ~  log_mean_lag1 + s(trend) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod5)

female_boreal_mod6 <- gam(log_mean ~  log_mean_lag1 + s(trend_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod6)

female_boreal_mod7 <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = female_dat)
summary(female_boreal_mod7)
# plot(female_boreal_mod7, resid = T, pch = 19, se = T)

boreal.aicc <- MuMIn::AICc(female_boreal_mod1, female_boreal_mod2, female_boreal_mod3, female_boreal_mod4, 
                           female_boreal_mod5, female_boreal_mod6, female_boreal_mod7)

# same model structures with bottom temperature as the covariate

female_temp_mod1 <- gam(log_mean ~  log_mean_lag1 + s(temp2_lag1) + s(year, k = 4), data = female_dat)
summary(female_temp_mod1)

female_temp_mod2 <- gam(log_mean ~ log_mean_lag1 + s(temp2) + s(year, k = 4), data = female_dat)
summary(female_temp_mod2)

female_temp_mod3 <- gam(log_mean ~ log_mean_lag1 + s(temp3) + s(year, k = 4), data = female_dat)
summary(female_temp_mod3)

female_temp_mod4 <- gam(log_mean ~  log_mean_lag1 + s(temp2_lag2) + s(year, k = 4), data = female_dat)
summary(female_temp_mod4)

female_temp_mod5 <- gam(log_mean ~  log_mean_lag1 + s(temp) + s(year, k = 4), data = female_dat)
summary(female_temp_mod5)

female_temp_mod6 <- gam(log_mean ~  log_mean_lag1 + s(temp_lag1) + s(year, k = 4), data = female_dat)
summary(female_temp_mod6)

female_temp_mod7 <- gam(log_mean ~  log_mean_lag1 + s(temp3_lag1) + s(year, k = 4), data = female_dat)
summary(male_temp_mod7)

temp.aicc <- MuMIn::AICc(female_temp_mod1, female_temp_mod2, female_temp_mod3, female_temp_mod4, female_temp_mod5, 
                         female_temp_mod6, female_temp_mod7)

# get minimum AICc value for the different models
min.aicc <- min(boreal.aicc$AICc, temp.aicc$AICc)

# combine to plot
plot.female.aicc <- rbind(boreal.aicc, temp.aicc) %>%
  mutate(`Covariate` = rep(c("Borealization index", "Bottom temperature"), each = 7),
         `Delta AICc` = AICc - min.aicc,
         plot.lags = rep(c("Year -1, -2",
                           "Year 0, -1",
                           "Year 0, -1, -2",
                           "Year -2, -3",
                           "Year 0",
                           "Year -1",
                           "Year -1, -2, -3"), 2))

plot.female.aicc <- plot.female.aicc %>%
  arrange(`Covariate`, desc(`Delta AICc`))

# add plot order (order of delta-AICc for boreal models)
plot.female.aicc$order <- order(plot.female.aicc$AICc[plot.female.aicc$Covariate == "Borealization index"])

plot.labels <- reorder(plot.female.aicc$plot.lags[1:7], plot.female.aicc$order[1:7])

plot.female.aicc$covariate.order <- if_else(plot.female.aicc$Covariate == "Borealization index", 2, 1)

plot.female.aicc$Covariate <- reorder(plot.female.aicc$Covariate, plot.female.aicc$covariate.order)

female.plot <- ggplot(plot.female.aicc, aes(order, `Delta AICc`, color = `Covariate`)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:7,
                     labels = plot.female.aicc$plot.lags[7:1]) +
  scale_color_manual(values = cb[c(8,6)]) + 
  labs(x = "Covariate lag", y = expression(Delta~"AICc")) +
  theme(axis.text.x  = element_text(angle=60, hjust=1)) +
  ggtitle("a) Female")

female.plot

# combine
png("./figs/Extended_Data_Fig_2.png", width = 8, height = 4, units = 'in', res = 300)

ggpubr::ggarrange(female.plot, male.plot, ncol = 2, common.legend = T)

dev.off()

## how robust are the best models to the imputation of 2020 abundance as a covariate? -----------

## begin with males

# use male_dat, loaded above

# create a sequence of possible values for 2020 abundance (2021 lagged value)
male_range <- male_dat$log_mean_lag1[male_dat$year == 2020] - male_dat$log_mean_lag1[male_dat$year == 2022]

possible_lag1_male_abundance <- seq(from = male_dat$log_mean_lag1[male_dat$year == 2020] - 0.1*male_range,
                                    to = male_dat$log_mean_lag1[male_dat$year == 2020] - 0.9*male_range,
                                    length.out = 100)


# create object to catch results
male_out <- data.frame()

for(i in 1:length(possible_lag1_male_abundance)){
  # i <- 1
  # reset the imputation
  male_dat$log_mean_lag1[male_dat$year == 2021] = possible_lag1_male_abundance[i]
  
  # fit the best model (#7) to this new version of the data
  mod <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = male_dat)
  
  male_out <- rbind(male_out,
                    data.frame(Sex = "Male",
                               Rsq = summary(mod)$r.sq,
                               log_mean_lag1 = possible_lag1_male_abundance[i],
                               P_lag_cpue = summary(mod)$p.table[2,4],
                               P_trend3_lag1 = summary(mod)$s.table[1,4],
                               P_year = summary(mod)$s.table[2,4]))
  
}

ggplot(male_out, aes(log_mean_lag1, P_trend3_lag1)) +
  geom_line()

## females

# use female_dat, loaded above

# create a sequence of possible values for 2020 abundance (2021 lagged value)
female_range <- female_dat$log_mean_lag1[female_dat$year == 2020] - female_dat$log_mean_lag1[female_dat$year == 2022]

possible_lag1_female_abundance <- seq(from = female_dat$log_mean_lag1[female_dat$year == 2020] - 0.1*female_range,
                                    to = female_dat$log_mean_lag1[female_dat$year == 2020] - 0.9*female_range,
                                    length.out = 100)


# create object to catch results
female_out <- data.frame()

for(i in 1:length(possible_lag1_female_abundance)){
  # i <- 1
  # reset the imputation
  female_dat$log_mean_lag1[female_dat$year == 2021] = possible_lag1_female_abundance[i]
  
  # fit the best model (#7) to this new version of the data
  mod <- gam(log_mean ~  log_mean_lag1 + s(trend3_lag1) + s(year, k = 4), data = female_dat)
  
  female_out <- rbind(female_out,
                    data.frame(Sex = "Female",
                               Rsq = summary(mod)$r.sq,
                               log_mean_lag1 = possible_lag1_female_abundance[i],
                               P_lag_cpue = summary(mod)$p.table[2,4],
                               P_trend3_lag1 = summary(mod)$s.table[1,4],
                               P_year = summary(mod)$s.table[2,4]))
  
}

ggplot(female_out, aes(log_mean_lag1, P_trend3_lag1)) +
  geom_line()

# combine
both_out <- rbind(male_out, female_out)

# data frame for denoting actual values

actual <- data.frame(Sex = c("Female", "Male"),
                     log_mean_lag1 = c(mean(c(female_dat$log_mean_lag1[female_dat$year == 2020], female_dat$log_mean_lag1[female_dat$year == 2022])),
                                       mean(c(male_dat$log_mean_lag1[male_dat$year == 2020], male_dat$log_mean_lag1[male_dat$year == 2022]))))

ggplot(both_out, aes(log_mean_lag1, P_trend3_lag1)) +
  geom_line() +
  facet_wrap(~Sex, scales = "free") +
  geom_vline(data = actual, aes(xintercept = log_mean_lag1), lty = 2) +
  labs(x = "Lag1 log(CPUE)",
       y = "P(borealization index)")

ggsave("./figs/Extended_Data_Fig_7.png", width = 6, height = 3, units = 'in')

