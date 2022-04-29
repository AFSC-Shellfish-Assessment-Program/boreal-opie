## modify model to estimate missing abundance in 2020 ----------------------------

# reload data

# load borealization DFA trend 

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

# clean up and combine

abundance <- abundance %>%
  rename(year = AKFIN_SURVEY_YEAR) %>%
  mutate(log_abundance = log(ABUNDANCE_MIL), .keep = "unused")

dat <- left_join(trend, abundance)

plot.dat <- dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)


dat <- dat %>%
  mutate(trend2 = zoo::rollmean(trend, 2, fill = "NA", align = "right"),
         log_abundance_lag1 = lead(log_abundance, n = 1),
         abundance_change = log_abundance_lag1 - log_abundance) %>%
  filter(year %in% 1980:2020)


ggplot(dat, aes(abundance_change)) +
  geom_histogram(fill = "grey", color = "black")

# # set up model with mi
# 
# form <- bf(log_abundance_lag1 | mi() ~ mi(log_abundance) + s(trend)) + 
#      bf(log_abundance | mi() ~ 1) + set_rescor(FALSE)
# 
# ## fit
# mi_brm <- brm(form,
#                  data = dat,
#                  cores = 4, chains = 4, iter = 2000,
#                  save_pars = save_pars(all = TRUE),
#                  control = list(adapt_delta = 0.999, max_treedepth = 10))
# 
# saveRDS(mi_brm, file = "output/mi_brm.rds")
# 
# summary(mi_brm)
# 
# ## plot s(z)
# cs = conditional_smooths(mi_brm)
# plot(cs, ask = FALSE)
# 
# ## plot estimated missing values
# post = as.matrix(mi_brm)
# hist(post[ , "Ymi_logabundancelag1[40]"])
# hist(post[ , "Ymi_logabundance[41]"]) 
# 
# # seems that there is added uncertainty because the same value is estimated twice? (once at lag0 and once at lag1)?
# 
# # plot estimate value compared with observed time series
# estimated <- data.frame(year = 2020, 
#                         log_abundance = mean(post[ , "Ymi_logabundance[41]"]),
#                         LCI = quantile(post[ , "Ymi_logabundance[41]"], 0.025),
#                         UCI = quantile(post[ , "Ymi_logabundance[41]"], 0.975))
# 
# 
# abundance.plot <- abundance %>%
#   mutate(LCI = NA,
#          UCI = NA)
# 
# abundance.plot <- rbind(abundance.plot, estimated)
# 
# ggplot(abundance.plot, aes(year, log_abundance)) +
#   geom_line() +
#   geom_point() +
#   geom_errorbar(aes(ymin = LCI, ymax = UCI))
# 

# estimate with mice
library(mice)

# load prediction time series 
pred_dat <- read.csv("./Data/imputation_data.csv")

# clean up and remove survey data - can consider these in a lagged application if needed
pred_dat <- pred_dat[1:42,1:8] %>%
  select(-survey_female_mat_biomass, -survey_male_mat_biomass) 
  
# and log transform

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
  geom_errorbar(aes(ymin = LCI, ymax = UCI))

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
              cores = 4, chains = 4, iter = 2000,
              save_pars = save_pars(all = TRUE),
              control = list(adapt_delta = 0.9999, max_treedepth = 10))

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
  labs(x = "Borealization trend", y = "Log abundance") +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g1)

ggsave("./Figs/mice_borealization_abundance_regression.png", width = 6, height = 4, units = 'in')
