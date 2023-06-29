## estimate preindustrial and historical 
## probabilities for given 3-yr rolling mean sst anomalies

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./scripts/stan_utils.R")

theme_set(theme_bw())


# load preindustrial outcomes
preindustrial_files <- list.files("./output/preindustrial_borealization_outcomes")

preindustrial <- data.frame() 

for(i in 1:length(preindustrial_files)){
    # i <- 1
    temp <- read.csv(paste("./output/preindustrial_borealization_outcomes/", preindustrial_files[i], sep = ""))
    
    preindustrial <- rbind(preindustrial, temp)
    
}


# load historical outcomes
historical_files <- list.files("./output/historical_borealization_outcomes")

historical <- data.frame() 

for(i in 1:length(historical_files)){
    # i <- 1
    temp <- read.csv(paste("./output/historical_borealization_outcomes/", historical_files[i], sep = ""))
    
    historical <- rbind(historical, temp)
    
}

# clean up names
historical <- historical %>%
    rename(observed.year = ersst.year,
           annual.event = annual.events)

# combine
outcomes <- rbind(preindustrial, historical)

# load model weights 
weights <- read.csv("./data/normalized_CMIP6_weights.csv")

weights <- weights %>%
    filter(region == "Eastern_Bering_Sea") %>%
   select(model, normalized_weight) %>%
   rename(model_weight = normalized_weight)

outcomes <- left_join(outcomes, weights) %>%
    select(model, model_weight, period, observed.year, borealization_index, annual.event)

    # # drop NAs due to smoothing
# drop <- is.na(outcomes$annual.anomaly.3yr)
# outcomes <- outcomes[!drop,]

## brms: setup ---------------------------------------------

library(data.table)
dt <- as.data.table(outcomes)
dt[model == "ACCESS-CM2" & period == "historical", .(N = .N), by = .(period, model, observed.year)]
dt[model == "ACCESS-CM2" & period == "historical" & observed.year == 2017, ]

prop <- dt[ , .(count = sum(annual.event, na.rm = T),
                N = .N,
                borealization_index = unique(borealization_index),
                model_weight = unique(model_weight)),
           by = .(period, model, observed.year)]
prop[ , prop := count / N]
prop[ , n_comparison := if_else(period == "preindustrial", N / 10, N / 40)]
prop[ , count_use := round(prop * n_comparison)]
prop[ , model_fac := as.factor(model)]
prop[ , year_fac := as.factor(observed.year)]


head(prop[period == "preindustrial", ], n = 20)
head(prop[period == "historical", ], n = 20)

g <- ggplot(prop) +
    geom_point(aes(x = year_fac, y = prop, color = model_fac)) +
    facet_wrap( ~ period) +
    theme(legend.position="none")
print(g)

# now need to sum results for each year
prop_summary <- prop %>%
    group_by(period, model, observed.year) %>%
    summarize(borealization_index = mean(borealization_index),
              model_weight = mean(model_weight),
              count_use = mean(count_use),
              n_comparison = mean(n_comparison))

prop_summary <- prop_summary %>%
    mutate(model_fac = as.factor(model))

# # reduce to 10 best models to ease computational load
# weights <- weights %>%
#     mutate(model = reorder(model, desc(model_weight)))
# 
# ggplot(weights, aes(model, model_weight)) +
#     geom_point() +
#     geom_line() +
#     theme(axis.text.x = element_text(angle = 90))
# 
# # 10 models with weights > 1 and year >= 1972 (period for which we have borealization data)!
# 
# prop_best_models <- prop %>%
#     filter(model_weight > 1,
#            observed.year >= 1972)
# 
# unique(prop_best_models$model)
# 
# g <- ggplot(prop_best_models) +
#     geom_point(aes(x = year_fac, y = prop, color = model_fac)) +
#     facet_wrap( ~ period) +
#     theme(legend.position="none")
# print(g)

# run binomial brms model

form <-  bf(count_use | trials(n_comparison) + weights(model_weight, scale = TRUE) ~
            period + s(borealization_index, by = period, k = 6) +
            s(observed.year, by = period, k = 6) + (1 | model_fac))

# form <-  bf(count | trials(N) + weights(model_weight, scale = TRUE) ~
#                 period + s(annual.anomaly.1yr, by = period, k = 6) +
#                 s(ersst.year, by = period, k = 6) + (1 | model_fac))

far_borealization_brms <- brm(form,
                 data = prop_summary,
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 4000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 15))

# saveRDS(far_brms2, paste("./CMIP6/brms_output/", regions[i], "_binomial2.rds", sep = ""))

saveRDS(far_borealization_brms, "./output/far_borealization_brms.rds")

# diagnostics

check_hmc_diagnostics(far_borealization_brms$fit)
neff_lowest(far_borealization_brms$fit)
rhat_highest(far_borealization_brms$fit)
summary(far_borealization_brms)
bayes_R2(far_borealization_brms)
 
## plot

# load model object for predicting preindustrial and historical probabilities
mod <- readRDS("./output/far_borealization_brms.rds")

bayes_R2(mod)

far_pred <- data.frame()


# set up borealization
observed.borealization <- read.csv("./output/observed_borealization_from_posteriors.csv") %>% 
    filter(year %in% 1950:2022) %>%
    group_by(year) %>%
    summarize(borealization_index = mean(borealization_index))

## setup new data
nd <- data.frame(period = c("historical", "preindustrial"),
                 observed.year = rep(observed.borealization$year, each = 2),
                 borealization_index = rep(observed.borealization$borealization_index, each = 2),
                 n_comparison = 1000,
                 model_fac = NA)

nd_pre <- nd[nd$period == "preindustrial", ]
nd_his <- nd[nd$period == "historical", ]

## make predictions
## exclude random effects for model_fac
pre_pp <- posterior_epred(mod, newdata = nd_pre, re_formula = NA)
his_pp <- posterior_epred(mod, newdata = nd_his, re_formula = NA)

## Calc probabilities
## These are our posterior probabilities to use for FAR calculation
pre_prob <- pre_pp / unique(nd$n_comparison)
his_prob <- his_pp / unique(nd$n_comparison)


## Calc FAR
far <- 1 - (pre_prob / his_prob)
range(far, na.rm = TRUE)

far_pred <- rbind(far_pred,
                  data.frame(borealization_index = nd_pre$borealization_index,
                             year = nd_pre$observed.year,
                             far = apply(far, 2, mean),
                             lower = apply(far, 2, quantile, probs = 0.025),
                             upper = apply(far, 2, quantile, probs = 0.975)))

# far_pred <- far_pred %>%
#     filter(year %in% 1972:2022)


g2 <- ggplot(far_pred) +
    geom_line(aes(x = year, y = far), size = 0.2, color = "red3") +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
    ylab("Fraction of attributable risk") +
    xlab("Year")

print(g2)

# aside - calculate risk ratio

far_pred <- far_pred %>%
    mutate(rr = 1/(1-far))

g3 <- ggplot(far_pred) +
    geom_line(aes(x = year, y = rr), size = 0.2, color = "red3") +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
    ylab("Risk ratio") +
    xlab("Year")

print(g3)

# far_brms2 <- readRDS("./CMIP6/brms_output/Gulf_of_Alaska_binomial2.rds")

ce <- conditional_effects(far_borealization_brms)
plot(ce, ask = FALSE)

ce2 <- conditional_effects(far_borealization_brms, effect = "borealization_index:period", re_formula = NA)
plot(ce2)
head(ce2[[1]])
tail(ce2[[1]])



## setup new data
x <- ce2[[1]]
nd <- x[ , c("period", "borealization_index", "n_comparison", "model_fac")]
nd$n_comparison <- 1000
nd_pre <- nd[nd$period == "preindustrial", ]
nd_his <- nd[nd$period == "historical", ]

## make predictions
## exclude random effects for model_fac
pre_pp <- posterior_epred(far_borealization_brms, newdata = nd_pre, re_formula = NA)
his_pp <- posterior_epred(far_borealization_brms, newdata = nd_his, re_formula = NA)

## Calc probabilities
## These are our posterior probabilities to use for FAR calculation
pre_prob <- pre_pp / unique(nd$N)
his_prob <- his_pp / unique(nd$N)

range(pre_prob)
range(his_prob)
# plot(as.vector(pre_prob))
# plot(as.vector(his_prob))



## Plot conditional effects
pre_pred <- data.frame(period = "preindustrial",
                       annual.anomaly.3yr = nd_pre$annual.anomaly.3yr,
                       prob = apply(pre_prob, 2, mean),
                       lower = apply(pre_prob, 2, quantile, probs = 0.025),
                       upper = apply(pre_prob, 2, quantile, probs = 0.975))
his_pred <- data.frame(period = "historical",
                       annual.anomaly.3yr = nd_his$annual.anomaly.3yr,
                       prob = apply(his_prob, 2, mean),
                       lower = apply(his_prob, 2, quantile, probs = 0.025),
                       upper = apply(his_prob, 2, quantile, probs = 0.975))
pred <- rbind(pre_pred, his_pred)

g <- ggplot(pred) +
    geom_line(aes(x = annual.anomaly.3yr, y = prob, color = period), size = 0.8) +
    geom_ribbon(aes(x = annual.anomaly.3yr, ymin = lower, ymax = upper, fill = period), alpha = 0.15) +
    scale_color_manual(values = c("tomato", "steelblue")) +
    scale_fill_manual(values = c("tomato", "steelblue"))
print(g)



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
