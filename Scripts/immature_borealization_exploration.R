# exploratory correlations:
# different male size classes vs borealization index

library(tidyverse)

theme_set(theme_bw())

dat <- read.csv("./output/male3059_drop10_df.csv", row.names = 1)

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t,
         bor_index = estimate)

dat <- dat %>%
  select(year, imp_log_mean) %>%
  left_join(., trend)


ccf(dat$bor_index, dat$imp_log_mean)

dat$lag_bor_index = lag(dat$bor_index)

ccf(dat$lag_bor_index[2:47], dat$imp_log_mean[2:47])

ggplot(dat, aes(bor_index, imp_log_mean)) +
  geom_point()

ggplot(dat, aes(lag_bor_index, imp_log_mean)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))

# bring in 60-95 mm males

dat2 <- read.csv("./output/male6095_drop5_df.csv", row.names = 1)

dat2 <- dat2 %>%
  mutate(lead_imp_log_mean = lead(imp_log_mean)) %>%
  select(year, imp_log_mean, lead_imp_log_mean) %>%
  rename(log_male_60_95 = imp_log_mean,
         lead_male_60_95 = lead_imp_log_mean)

dat <- dat %>%
  rename(log_male_30_59 = imp_log_mean) %>%
  left_join(., dat2)

ggplot(dat, aes(log_male_30_59, lead_male_60_95)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))

ggplot(dat, aes(bor_index, log_male_60_95)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))

ggplot(dat, aes(lag_bor_index, log_male_60_95)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))

ggplot(dat, aes(lag_bor_index, log_male_30_59)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))

ggplot(dat, aes(bor_index, log_male_30_59)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 6))
