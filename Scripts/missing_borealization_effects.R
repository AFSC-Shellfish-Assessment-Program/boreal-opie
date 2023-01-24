## modify model to estimate missing abundance in 2020 ----------------------------
library(tidyverse)
library(mice)
library(rstan)
library(brms)
library(bayesplot)
library(mgcv)

source("./scripts/stan_utils.R")

theme_set(theme_bw())


# reload data

# load borealization DFA trend

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

# abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

abun1 <- read.csv("./output/male3059_drop10_df.csv", row.names = 1) %>%
  mutate(size = "30-59")

abun2 <- read.csv("./output/male6095_drop5_df.csv", row.names = 1) %>%
  mutate(size = "60-95")

# clean up and combine

abundance <- rbind(abun1, abun2) %>%
  rename(log_mean = imp_log_mean) %>%
  dplyr::select(year, log_mean, size) %>%
  pivot_longer(cols = c(-year, -size)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = size, values_from = value)

dat <- left_join(trend, abundance)

# cross-correlation for age classes: strongest at lag 1
ccf(as.vector(dat[dat$year %in% 1975:2019,3]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf

# cross-correlation for borealization index v. small size class - strongest at lag 1!
ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,3]))$acf

# cross-correlation for borealization index v. large size class - strongest at lag 0!
ccf(as.vector(dat[dat$year %in% 1975:2019,2]), as.vector(dat[dat$year %in% 1975:2019,4]))$acf

dat$trend2 <- zoo::rollmean(dat$trend, 2, fill = NA, align = "right")
dat$trend3 <- zoo::rollmean(dat$trend, 3, fill = NA, align= "right")

str(dat)

# cross-correlation for smoothed borealization index v. small size class - strongest at lag 0
ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,3]))

# cross-correlation for smoothed borealization index v. small size class - strongest at lag 0
ccf(as.vector(dat[dat$year %in% 1975:2019,5]), as.vector(dat[dat$year %in% 1975:2019,4]))

# set up lag trend for plotting
dat <- dat %>%
  mutate(trend_lag1 = lag(trend, 1),
         trend2_lag1 = lag(trend2, 1))
  
ggplot(dat, aes(trend, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend_lag1, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2, `30-59`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend_lag1, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend2_lag1, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

ggplot(dat, aes(trend3, `60-95`)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)

plot.dat <- dat %>%
  pivot_longer(cols = -year)


# exploratory GAMs - 60-95 size class as the lag of 3-50 + borealization effect
dat <- dat %>%
  rename(small_male = `30-59`,
         large_male = `60-95`) %>%
  mutate(small_male_lag1 = lag(small_male, 1))
 
# model small males 
sm_mod1 <- gam(small_male ~ s(trend_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod1)
plot(sm_mod1, resid = T, se = F, pch = 19)

sm_mod2 <- gam(small_male ~ s(trend2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod2)
plot(sm_mod2, resid = T, se = F, pch = 19)

sm_mod3 <- gam(small_male ~ s(trend, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod3)
plot(sm_mod3, resid = T, se = F, pch = 19)

sm_mod4 <- gam(small_male ~ s(trend2_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod4)
plot(sm_mod4, resid = T, se = F, pch = 19)

ggplot(dat, aes(trend2_lag1, small_male)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)


sm_mod5 <- gam(small_male ~ s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod5)
plot(sm_mod5, resid = T, se = F, pch = 19)

dat <- dat %>%
  mutate(trend_lag2 = lag(trend, 2))

sm_mod6 <- gam(small_male ~ s(trend_lag2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod6)
plot(sm_mod6, resid = T, se = F, pch = 19)

sm_mod7 <- gam(small_male ~ s(trend_lag1, k = 4) + s(trend_lag2, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod7)
plot(sm_mod7, resid = T, se = F, pch = 19)


MuMIn::AICc(sm_mod1, sm_mod2, sm_mod3, sm_mod4, sm_mod5, sm_mod6, sm_mod7) # best model is sm_mod4

## brms version of best small male model--------------------

form <- bf(small_male ~ s(trend2_lag1, k = 4))

## fit
mice_brm <- brm(form,
                         data = dat,
                         cores = 4, chains = 4, iter = 3000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(mice_brm, file = "output/mice_sm_male_brm_no_imputation.rds")

# diagnostics
mice_brm <- readRDS("./output/mice_sm_male_brm_no_imputation.rds")
check_hmc_diagnostics(mice_brm$fit)
neff_lowest(mice_brm$fit)
rhat_highest(mice_brm$fit)

bayes_R2(mice_brm)

plot(conditional_effects(mice_brm), ask = FALSE)

# y <- y %>%
#   group_by(year) %>%
#   summarise(log_abundance_lead1 = mean(log_abundance_lead1))
# 
# y <- y$log_abundance_lead1
# 
# yrep_mice_brm  <- fitted(mice_brm, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_mice_brm[sample(nrow(yrep_mice_brm), 25), ]) +
#   ggtitle("mice_brm")
# 
# ggsave("./Figs/mice_brms_ppc.png", width = 6, height = 4, units = 'in')
# 
# trace_plot(mice_brm$fit)

# plot mice_brm

## 95% CI
ce1s_1 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(trend$trend[trend$year %in% 1980:2020]), amount = 0.01),
                          rep(NA, 100-length(trend$trend[trend$year %in% 1980:2020])))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Log abundance") +
  geom_rug(aes(x=rug.anom, y=NULL))

print(g1)

ggsave("./Figs/mice_borealization_abundance_regression.png", width = 6, height = 4, units = 'in')
# large male models

large_mod1 <- gam(large_male ~ s(small_male_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod1)
plot(large_mod1, resid = T, se = F, pch = 19)


large_mod2 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod2)
plot(large_mod2, resid = T, se = F, pch = 19)


large_mod3 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend2, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod3)
plot(large_mod3, resid = T, se = F, pch = 19)


large_mod4 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod4)
plot(large_mod4, resid = T, se = F, pch = 19)

large_mod5 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod5)
plot(large_mod5, resid = T, se = F, pch = 19)

large_mod6 <- gam(large_male ~ s(small_male_lag1, k = 4) + s(trend2_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod6)
plot(large_mod6, resid = T, se = F, pch = 19)

MuMIn::AICc(large_mod1, large_mod2, large_mod3, large_mod4, large_mod5, large_mod6) # mod 2 is best

# now models with only trend
large_mod7 <- gam(large_male ~ s(trend_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod7)
plot(large_mod7, resid = T, se = F, pch = 19)

large_mod8 <- gam(large_male ~ s(trend2, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod8)
plot(large_mod8, resid = T, se = F, pch = 19)

large_mod9 <- gam(large_male ~ s(trend, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod9)
plot(large_mod9, resid = T, se = F, pch = 19)

large_mod10 <- gam(large_male ~ s(trend2_lag1, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod10)
plot(large_mod10, resid = T, se = F, pch = 19)

large_mod11 <- gam(large_male ~ s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(large_mod11)
plot(large_mod11, resid = T, se = F, pch = 19)

MuMIn::AICc(large_mod7, large_mod8, large_mod9, large_mod10, large_mod11) # mod 10


ggplot(dat, aes(trend2_lag1, small_male)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "gam", se = F)



sm_mod5 <- gam(small_male ~ s(trend3, k = 4), dat = dat, na.action = "na.omit")
summary(sm_mod5)
plot(sm_mod5, resid = T, se = F, pch = 19)


ggplot(plot.dat, aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)





# ## brms imputation -----------------------------------------
# 
# # load prediction time series
# pred_dat <- read.csv("./Data/imputation_data.csv")
# 
# # clean up and remove survey data - can consider these in a lagged application if needed
# pred_dat <- pred_dat[1:42,1:8] %>%
#   select(-survey_female_mat_biomass, -survey_male_mat_biomass)
# 
# # and log transform
# pred_dat[,2:6] <- log(pred_dat[,2:6])
# 
# ## combine all data sources
# mi_dat <- rbind(abundance, data.frame(year = 2020, log_abundance = NA))
# mi_dat <- mi_dat[order(mi_dat$year), ]
# mi_dat <- left_join(mi_dat, trend, by = "year")
# mi_dat$log_abundance_lead1 = lead(mi_dat$log_abundance)
# 
# ## we don't want to include 2022 abundance (log_abundance_lead1 for 2021)
# ## the fitted value for 2020 gives us our predicted 2021 log abundance
# mi_dat <- mi_dat[mi_dat$year %in% 1980:2020, ]
# 
# # imputed.data[[1]]
# 
# 
# mi_dat_pred <- left_join(mi_dat, pred_dat, by = "year")
# 
# 
# ## set up model with mi
# mi_form <- bf(log_abundance_lead1 | mi() ~ mi(log_abundance) + s(trend)) +
#      bf(log_abundance | mi() ~ model_female_mat_biomass +
#                                model_male_mat_biomass +
#                                model_3plus_pollock_biomass +
#                                model_female_plaice_biomass +
#                                model_2plus_yellowfin_biomass) + set_rescor(FALSE)
# 
# ## fit
# mi_brm <- brm(mi_form,
#               data = mi_dat_pred,
#               cores = 4, chains = 4, iter = 2000,
#               save_pars = save_pars(all = TRUE),
#               control = list(adapt_delta = 0.999, max_treedepth = 10))
# 
# saveRDS(mi_brm, file = "output/mi_brm.rds")
# 
# mi_brm <- readRDS("output/mi_brm.rds")
# 
# summary(mi_brm)
# 
# ## plot s(z)
# cs <- conditional_effects(mi_brm, effects = "trend")
# cs[[2]] <- NULL
# plot(cs, ask = FALSE)
# 
# 
# ## plot estimated missing values
# post = as.matrix(mi_brm)
# hist(post[ , "Ymi_logabundancelead1[38]"])
# hist(post[ , "Ymi_logabundance[39]"])
# 
# # seems that there is added uncertainty because the same value is estimated twice? (once at lag0 and once at lag1)?
# 
# # plot estimate value compared with observed time series
# estimated <- data.frame(year = 2020,
#                         log_abundance = mean(post[ , "Ymi_logabundance[39]"]),
#                         LCI = quantile(post[ , "Ymi_logabundance[39]"], 0.025),
#                         UCI = quantile(post[ , "Ymi_logabundance[39]"], 0.975))
# 
# 
# abundance.plot <- abundance %>%
#   mutate(LCI = NA,
#          UCI = NA)
# 
# abundance.plot <- rbind(abundance.plot, estimated)
# 
# g <- ggplot(abundance.plot, aes(year, log_abundance)) +
#   geom_line() +
#   geom_point() +
#   geom_errorbar(aes(ymin = LCI, ymax = UCI))
# print(g)
# ggsave("./Figs/brms_imputed_survey_abundance.png", width = 5, height = 3, units = 'in')



# estimate with mice ---------------------------------------
library(mice)

# load prediction time series
# pred_dat1 <- read.csv("./Data/imputation_data.csv")

# pred_dat <- read.csv("./output/bycatch_directfish_groundfish_ts.csv", row.names = 1) %>%
#   left_join(.,dat)

# # check cross correlations
# 
# ccf(pred_dat$male3059.directfish, pred_dat$`30-59`)
# 
# ccf(pred_dat$male3059.bycatch, pred_dat$`30-59`)
# 
# ccf(pred_dat$male6095.directfish, pred_dat$`60-95`)
# 
# ccf(pred_dat$male6095.bycatch, pred_dat$`60-95`)
# 
# # clean up and remove survey data - can consider these in a lagged application if needed
# pred_dat <- pred_dat[1:42,1:8] %>%
#   select(-survey_female_mat_biomass, -survey_male_mat_biomass)

pred.dat1 <- read.csv("./data/gf_imputation_data.csv")

# and log transform 
pred.dat1[,c(2:6)] <- log(pred.dat1[,c(2:6)])

pred.dat2 <- read.csv("./output/directed_bycatch_all_years.csv", row.names = 1)


# combine data for ccf analysis
cor.dat <- dat[,c(1,3,4)] %>%
  left_join(.,pred.dat1) %>%
  left_join(.,pred.dat2)

cor.dat <- cor.dat[cor.dat$year %in% 1975:2019,] # using base b/c tidyverse acting up

# check cross correlations
# small males first
ccf(cor.dat$small_male, cor.dat$pollockR1)$acf # some information at lag 2 (opilio leads pollock)
ccf(cor.dat$small_male[cor.dat$year >= 1978], cor.dat$codR0[cor.dat$year >= 1978])$acf # six year lag peak implausible! 
# lag 2 (cod leads opilio) is best, most usable
ccf(cor.dat$small_male[cor.dat$year >= 1978], cor.dat$codR0[cor.dat$year >= 1978])$lag

ccf(cor.dat$small_male, cor.dat$pollock3plus)$acf # 2 year lag (pollock leads opilio)
ccf(cor.dat$small_male[cor.dat$year >= 1978], cor.dat$codbiomass[cor.dat$year >= 1978])$acf # cod leads opilio by 4 years
ccf(cor.dat$small_male, cor.dat$ylfinbiomass)$acf # lag 0?

ccf(cor.dat$small_male[cor.dat$year >= 1995], cor.dat$large_log_cpue_directed[cor.dat$year >= 1995])$acf # lag 0
ccf(cor.dat$small_male[cor.dat$year >= 1995], cor.dat$large_log_cpue_bycatch[cor.dat$year >= 1988])$acf # lag 0

ccf(cor.dat$small_male[cor.dat$year >= 1995], cor.dat$small_log_cpue_directed[cor.dat$year >= 1995])$acf # lag 0 - weak / don't use!
ccf(cor.dat$small_male[cor.dat$year >= 1995], cor.dat$small_log_cpue_bycatch[cor.dat$year >= 1988])$acf # lag 3 - implausible / don't use!

# create a df for predicting small male abundance (at correct lags)
imp_small <- dat %>%
  dplyr::select(year, small_male) %>%
  left_join(.,pred.dat1) %>%
  left_join(.,pred.dat2)

# set up correct lags
imp_small <- imp_small %>%
  mutate(pollock_R1_lag3 = lag(pollockR1, 3)) %>% #  small males 2020, pollockR1 2017
  mutate(cod_R0_lead2 = lead(codR0, 2)) %>% # small males 2020, cod R0 2022
  mutate(codbiomass_lag3 = lag(codbiomass, 3)) %>% # small males 2020, cod biomass 2017
  mutate(pollockbiomass_lag2 = lag(pollock3plus, 2)) %>% # small male 2020, pollock biomass 2022
  rename(ylfinbiomass_lag0 = ylfinbiomass) %>%
  rename(large_male_directed_lag0 = large_log_cpue_directed) %>%
  rename(large_male_bycatch_lag0 = large_log_cpue_bycatch) %>%
  dplyr::select(-pollockR1, -codR0, -codbiomass, -pollock3plus, -small_log_cpue_bycatch, -small_log_cpue_directed) %>%
  filter(year >= 1975) # drop early years with NA for survey
  
  
 cor(imp_small, use = "p") 
 
 range(imp_small$year) # 1975-2022
 
 # remove year  
imp_small <- imp_small %>%
  dplyr::select(-year)
 
imp_obj_small <- mice(data = imp_small, method = "norm", m = 100)
 
 # pull out 2020 estimates
 estimated_2020 <- NA
 
 for(i in 1:100){
   # i <- 1
   estimated_2020[i] <- complete(imp_obj_small, i)$small_male[46]
   
 }
 
 
 # plot estimate value compared with observed time series
 estimated <- data.frame(year = 2020,
                         log_abundance = mean(estimated_2020),
                         LCI = quantile(estimated_2020, 0.025),
                         UCI = quantile(estimated_2020, 0.975))
 
 
 abundance.plot <- abundance %>%
   dplyr::rename(log_abundance = `30-59`) %>%
   dplyr::select(year, log_abundance) %>%
   mutate(LCI = NA,
          UCI = NA)
 
 abundance.plot <- rbind(abundance.plot, estimated)
 
 ggplot(abundance.plot, aes(year, log_abundance)) +
   geom_line() +
   geom_point() +
   geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
   ylab("log(immature snow crab abundance)") +
   theme(axis.title.x = element_blank())
 
 # save plot
 ggsave("./Figs/small_male_imputed_survey_abundance.png", width = 5, height = 3, units = 'in')

 ## process imp and pass mice imputations to brms ------------------
 
 imputed.data <- list()
 
 x.dat <- dat %>%
   dplyr::select(year, trend2_lag1)
 
 # get a df to plot y values for diagnostics
 y <- data.frame()
 
 for(i in 1:100){
   
   temp <- complete(imp, i)[1]
   
   y <- rbind(y, data.frame(
     year = 1975:2022,
     small_male_log_cpue = temp$small_male)
   )
   
   # add borealization trend
   temp <- left_join(y, x.dat)
   
   imputed.data[[i]] <- temp
   
 }
 
 # set up and run brms
 
 form <- bf(small_male_log_cpue ~ s(trend2_lag1, k = 4))
 
 ## fit
 mice_brm <- brm_multiple(form,
                          data = imputed.data,
                          cores = 4, chains = 4, iter = 3000,
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.999, max_treedepth = 12))
 
 saveRDS(mice_brm, file = "output/mice_sm_male_brm.rds")
 
 summary(mice_brm)
 
 # diagnostics
 mice_brm <- readRDS("./output/mice_sm_male_brm.rds")
 check_hmc_diagnostics(mice_brm$fit)
 neff_lowest(mice_brm$fit)
 rhat_highest(mice_brm$fit)
 
 bayes_R2(mice_brm)
 
 plot(conditional_effects(mice_brm), ask = FALSE)
 
 y <- y %>%
   group_by(year) %>%
   summarise(log_abundance_lead1 = mean(log_abundance_lead1))
 
 y <- y$log_abundance_lead1
 
 yrep_mice_brm  <- fitted(mice_brm, scale = "response", summary = FALSE)
 ppc_dens_overlay(y = y, yrep = yrep_mice_brm[sample(nrow(yrep_mice_brm), 25), ]) +
   ggtitle("mice_brm")
 
 ggsave("./Figs/mice_brms_ppc.png", width = 6, height = 4, units = 'in')
 
 trace_plot(mice_brm$fit)
 
 # plot mice_brm
 
 ## 95% CI
 ce1s_1 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                               probs = c(0.025, 0.975))
 ## 90% CI
 ce1s_2 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                               probs = c(0.05, 0.95))
 ## 80% CI
 ce1s_3 <- conditional_effects(mice_brm, effect = "trend2_lag1", re_formula = NA,
                               probs = c(0.1, 0.9))
 dat_ce <- ce1s_1$trend2_lag1
 dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
 dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
 dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
 dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
 dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
 dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]
 dat_ce[["rug.anom"]] <- c(jitter(unique(trend$trend[trend$year %in% 1980:2020]), amount = 0.01),
                           rep(NA, 100-length(trend$trend[trend$year %in% 1980:2020])))
 
 
 g1 <- ggplot(dat_ce) +
   aes(x = effect1__, y = estimate__) +
   geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
   geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
   geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
   geom_line(size = 1, color = "red3") +
   labs(x = "Borealization index", y = "Log abundance") +
   geom_rug(aes(x=rug.anom, y=NULL))
 
 print(g1)
 
 ggsave("./Figs/mice_borealization_abundance_regression.png", width = 6, height = 4, units = 'in')
 
 
## large males--------------------------
ccf(cor.dat$large_male, cor.dat$pollockR1)$acf # large males leads by 2 years
ccf(cor.dat$large_male[cor.dat$year >= 1978], cor.dat$codR0[cor.dat$year >= 1978])$acf # large male leads by 1 year
ccf(cor.dat$large_male, cor.dat$pollock3plus)$acf # lag 0
ccf(cor.dat$large_male[cor.dat$year >= 1978], cor.dat$codbiomass[cor.dat$year >= 1978])$acf # large males lead by 2 years
ccf(cor.dat$large_male, cor.dat$ylfinbiomass)$acf # lag 0

ccf(cor.dat$large_male[cor.dat$year >= 1995], cor.dat$male6095.directfish[cor.dat$year >= 1995])$acf # lag 0
ccf(cor.dat$large_male[cor.dat$year >= 1995], cor.dat$male6095.bycatch[cor.dat$year >= 1995])$acf # lag 0

ccf(cor.dat$large_male[cor.dat$year >= 1995], cor.dat$male3059.directfish[cor.dat$year >= 1995])$acf # lag 0 - weak / don't use!
ccf(cor.dat$large_male[cor.dat$year >= 1995], cor.dat$male3059.bycatch[cor.dat$year >= 1995])$acf # lag 0 - weak / don't use!


  

pred_dat[,2:6] <- log(pred_dat[,2:6])

dat <- abundance %>%
  rbind(.,
        data.frame(year = 2020,
                   log_abundance = NA)) %>%
  arrange(year) %>%
  left_join(., pred_dat)

# examine correlations!
cor(dat[,-1], use = "p")

imp <- mice(data = dat, method = "norm", m = 100)

# pull out 2020 estimates
estimated_2020 <- NA


for(i in 1:100){

estimated_2020[i] <- complete(imp, i)$log_abundance[41]

}


# plot estimate value compared with observed time series
estimated <- data.frame(year = 2020,
                        log_abundance = mean(estimated_2020),
                        LCI = quantile(estimated_2020, 0.025),
                        UCI = quantile(estimated_2020, 0.975))


abundance.plot <- abundance %>%
  mutate(LCI = NA,
         UCI = NA)

abundance.plot <- rbind(abundance.plot, estimated)

ggplot(abundance.plot, aes(year, log_abundance)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  ylab("log(immature snow crab abundance)") +
  theme(axis.title.x = element_blank())

# save plot
ggsave("./Figs/imputed_survey_abundance.png", width = 5, height = 3, units = 'in')

## process imp and pass mice imputations to brms ------------------

imputed.data <- list()

# get a df to plot y values for diagnostics
y <- data.frame()

for(i in 1:100){

  temp <- complete(imp, i)[,1:2] %>%
    mutate(log_abundance_lead1 = lead(log_abundance))

  y <- rbind(y, data.frame(
    year = temp$year[1:41],
    log_abundance_lead1 = temp$log_abundance_lead1[1:41])
  )

  # add borealization trend
  temp <- left_join(temp, trend)

  imputed.data[[i]] <- temp

  }

# set up and run brms

form <- bf(log_abundance_lead1 ~ log_abundance + s(trend))

## fit
mice_brm <- brm_multiple(form,
              data = imputed.data,
              cores = 4, chains = 4, iter = 3000,
              save_pars = save_pars(all = TRUE),
              control = list(adapt_delta = 0.99999999, max_treedepth = 12))

saveRDS(mice_brm, file = "output/mice_brm.rds")

summary(mice_brm)

# diagnostics
mice_brm <- readRDS("./output/mice_brm.rds")
check_hmc_diagnostics(mice_brm$fit)
neff_lowest(mice_brm$fit)
rhat_highest(mice_brm$fit)

bayes_R2(mice_brm)

plot(conditional_effects(mice_brm), ask = FALSE)

y <- y %>%
  group_by(year) %>%
  summarise(log_abundance_lead1 = mean(log_abundance_lead1))

y <- y$log_abundance_lead1

yrep_mice_brm  <- fitted(mice_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_mice_brm[sample(nrow(yrep_mice_brm), 25), ]) +
  ggtitle("mice_brm")

ggsave("./Figs/mice_brms_ppc.png", width = 6, height = 4, units = 'in')

trace_plot(mice_brm$fit)

# plot mice_brm

## 95% CI
ce1s_1 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(trend$trend[trend$year %in% 1980:2020]), amount = 0.01),
                          rep(NA, 100-length(trend$trend[trend$year %in% 1980:2020])))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Log abundance") +
  geom_rug(aes(x=rug.anom, y=NULL))

print(g1)

ggsave("./Figs/mice_borealization_abundance_regression.png", width = 6, height = 4, units = 'in')


# now predict for 2022 survey
new.dat <- data.frame(log_abundance = abundance$log_abundance[abundance$year == 2021],
                      trend = trend$trend[trend$year == 2021])


pred.2022 <- posterior_epred(mice_brm, newdata = new.dat)


mean(pred.2022)
LCI.80 <- exp(quantile(pred.2022, 0.1))
UCI.80 <- exp(quantile(pred.2022, 0.9))

LCI.95 <- exp(quantile(pred.2022, 0.025))
UCI.95 <- exp(quantile(pred.2022, 0.975))

overall.mean <- exp(mean(abundance$log_abundance[abundance$year %in% 1980:2019]))

LCI.80 / overall.mean
UCI.80 / overall.mean

LCI.95 / overall.mean
UCI.95 / overall.mean

# add blank year to abundance.plot
xtra.80 <- data.frame(year = 2022,
                   log_abundance = 5.5, # dummy value to make plot work
                   LCI = quantile(pred.2022, 0.1),
                   UCI = quantile(pred.2022, 0.9))

# abundance.plot <- rbind(abundance.plot, xtra)

# create data frame for 95% CI
xtra.95 <- data.frame(year = 2022,
                      log_abundance = 5.5, # dummy value to make plot work
                      LCI = quantile(pred.2022, 0.025),
                      UCI = quantile(pred.2022, 0.975))

ggplot(abundance.plot, aes(year, log_abundance)) +
  geom_line() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), color = "dark grey") +
  geom_point(size = 2) +
  geom_errorbar(data = xtra.95, aes(x = year, ymin = LCI, ymax = UCI), color = "gold") +
  geom_errorbar(data = xtra.80, aes(x = year, ymin = LCI, ymax = UCI), color = "firebrick") +

  ylab("Log abundance") +
  theme(axis.title.x = element_blank())

# save plot
ggsave("./Figs/imputed_and_predicted_survey_abundance.png", width = 5, height = 3, units = 'in')


## alternate model ---------------------------

form2 <- bf(log_abundance_lead1 ~ s(log_abundance) + s(trend))

## fit
mice_brm2 <- brm_multiple(form2,
                         data = imputed.data,
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.99999, max_treedepth = 16))

saveRDS(mice_brm2, file = "output/mice_brm2.rds")


# diagnostics
mice_brm2 <- readRDS("./output/mice_brm2.rds")
check_hmc_diagnostics(mice_brm2$fit)
neff_lowest(mice_brm2$fit)
rhat_highest(mice_brm2$fit)

bayes_R2(mice_brm2)

plot(conditional_effects(mice_brm2), ask = FALSE)
