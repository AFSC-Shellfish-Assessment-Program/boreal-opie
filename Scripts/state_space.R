## State-space model test for missing log_abundance


library(rstan)
library(data.table)
library(ggplot2)
theme_set(theme_bw())


## Read-in data
dat = read.csv("./Data/imm_abun.csv", row.names = 1)
names(dat) <- c("year", "abundance")
dat$log_abundance = log(dat$abundance)
trend = read.csv("./output/dfa_trend.csv")
trend = data.table(year = trend$t, trend = trend$estimate)

## Combine data
abund = c(dat$log_abundance[dat$year %in% 1980:2019], NA, dat$log_abundance[dat$year %in% 2021:2022])
opie = data.table(year = 1980:2022, log_abundance = abund)
opie$trend = trend$trend[trend$year %in% 1980:2022]

## Get data for stan
lst = list(N = nrow(opie),
           N_obs = sum(!is.na(opie$log_abundance)),
           N_mis = sum(is.na(opie$log_abundance)),
           y_obs = opie$log_abundance[!is.na(opie$log_abundance)],
           ii_obs = which(!is.na(opie$log_abundance)),
           ii_mis = as.array(which(is.na(opie$log_abundance))),
           c = opie$trend)

## This is a quick test model, probably would want a non-centered
## parameterization if we move forward with this
mod = rstan::stan_model(file = "./Scripts/stan/ssm.stan")

## Fit model
fit = rstan::sampling(mod,
                      data = lst,
                      iter = 5000,
                      cores = 4,
                      chains = 4,
                      seed = 42,
                      control = list(adapt_delta = 0.999, max_treedepth = 10))
saveRDS(fit, file = "./output/state_space.rds")
print(fit)

fit <- readRDS("output/state_space.rds")


## Look at state vector
mat = as.matrix(fit, pars = paste0("x[", 1:43, "]"))
med = apply(mat, 2, median)
lo = apply(mat, 2, quantile, probs = 0.025)
up = apply(mat, 2, quantile, probs = 0.975)

plot(opie$log_abundance)
plot(med, type = "o", col = 2)
segments(x0 = 1:43, y0 = lo, y1 = up, col = 2)


## Look at missing value estimate
mat = as.matrix(fit, pars = "y_mis[1]")
med = apply(mat, 2, median)
lo = apply(mat, 2, quantile, probs = 0.025)
up = apply(mat, 2, quantile, probs = 0.975)

plot(opie$log_abundance, type = 'o')
points(41, med, col = 2)
segments(x0 = 41, y0 = lo, y1 = up, col = 2)
