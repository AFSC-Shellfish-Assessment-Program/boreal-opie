# all the script required to make Figs 1 & 2 for proposed ms.

library(tidyverse)
library(MARSS)
library(sf)
library(ggmap)
library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Load map layers
map_layers <- readRDS("./Data/map_layers.rda")


## Trim survey grid to survey area
map_layers$survey.grid %>%
  st_transform(crs = st_crs(map_layers$survey.area)) %>%
  st_intersection(map_layers$survey.area) -> survey.grid

## Specify plot boundary, transform to map crs
data.frame(x = c(-180, -150), 
           y = c(54.5, 67)) %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(., crs = map_layers$crs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame() -> plot.boundary


## Plot 
fig1a <- ggplot() +
  geom_sf(data = map_layers$bathymetry, color=alpha("grey70")) +
  # geom_sf(data = map_layers$survey.grid, fill=NA, color=alpha("grey70"), linewidth = 1)+
  geom_sf(data = map_layers$survey.area, fill = alpha(cb[6], alpha=0.1), size = 0) +
  geom_sf(data = map_layers$akland, fill = "grey80", size=0.1) +
  geom_sf(data = map_layers$survey.area, fill = NA) +
  scale_x_continuous(breaks = c(-180, -175, -170, -165, -160, -155, -150), labels = paste0(c(180, 175, 170, 165, 160, 155, 150), "°W")) + 
  #scale_y_continuous(breaks = c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), labels = paste0(c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), "°N")) +
  coord_sf(xlim = plot.boundary$X,
           ylim = plot.boundary$Y)+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA, color = "grey70"),
        legend.key.size = unit(0.65,'cm'),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank()) 



# fig1a <- ggplot() + theme_void()

# abundance (Fig. 1b)

abundance <- read.csv("./output/imputed_male_30-95_imm_female_abundance.csv") 

# remove SD = 0 for plotting
change <- abundance$SD == 0
abundance$SD[change] <- NA

pos_dodge = position_dodge(width = 0.6)

fig1b <- ggplot(abundance, aes(as.numeric(year), abundance, color = sex)) +
  geom_hline(yintercept = 0) +
  geom_point(position = pos_dodge) +
  geom_line(position = pos_dodge) + 
  geom_errorbar(aes(ymin = abundance - 2*SD,
                    ymax = abundance + 2*SD), width = 0.5, position = pos_dodge) +
  scale_color_manual(values = cb[c(2,4)]) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.6, 0.8)) +
  labs(y = expression(Snow~crab~abundance~(10^9))) 


# fig1b <- ggplot(abundance, aes(year, abundance)) +
#   geom_col(fill = cb[6]) +
#   geom_errorbar(aes(ymin = abundance - 2*SD,
#                     ymax = abundance + 2*SD), position = dodge, width = 0.7) +
#   theme(axis.title.x = element_blank()) +
#   geom_hline(yintercept = 0) +
#   labs(y = expression(Snow~crab~abundance~(10^9))) 

# load DFA model
mod <- readRDS("./output/DFA_model.rds")

# load DFA data 
dat <- read.csv("./output/dfa time series.csv")

dfa.dat <- dat %>%
  dplyr::select(-order) %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(year) %>%
  dplyr::select(-year) %>%
  t()

colnames(dfa.dat) <- 1972:2022

# plot loadings and trend

CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                      mean=CI$par$Z[1:12],
                      upCI=CI$par.upCI$Z[1:12],
                      lowCI=CI$par.lowCI$Z[1:12])

dodge <- position_dodge(width=0.9)


plot.CI$names[c(3,4,5,6,8,12)] <- c("Jan-Feb ice", "Mar-Apr ice", "Bloom extent", "Phytoplankton size", "Bottom temp.", "Calanus")

plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:12])

fig1c <- ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[6]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=9), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

# plot trend
trend <- data.frame(t=1972:2022,
                    estimate=as.vector(mod$states),
                    conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                    conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))


fig1d <- ggplot(trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[6]) +
  geom_hline(yintercept = 0) +
  geom_point(color=cb[6]) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[6]) + xlab("") + ylab("Borealization index")

##
library(brms)

male_brm <- readRDS("./output/fit_male.rds")

## 95% CI
ce1s_1 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(male_brm, effect = "trend2_lag1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend2_lag1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend2_lag1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend2_lag1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend2_lag1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend2_lag1[["lower__"]]

fig1e <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (lag 1-2)", y = "Log CPUE") +
  annotate("text", x = -1.5, y = 7.5, label = "Male", size = 6)

##
female_brm <- readRDS("./output/fit_female.rds")

## 95% CI
ce1s_1 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(female_brm, effect = "trend3", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend3
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend3[["lower__"]]

fig1f <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index (lag 0-2)", y = "Log CPUE") +
  annotate("text", x = -1, y = 6.5, label = "Female", size = 6)


png("./figs/fig1.png", width = 12, height = 6, units = 'in', res = 300)

ggpubr::ggarrange(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f,
                  ncol = 3, nrow = 2,
                  labels = "auto")

dev.off()


### Fig. 2 ---------------------------------

# load borealization DFA trend and sst anomaly time series

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate)

sst <- read.csv("./Data/regional_north_pacific_ersst_anomaly_time_series.csv")

sst <- sst %>%
  filter(region == "Eastern_Bering_Sea") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(sst.anomaly = annual.anomaly.unsmoothed)

dat <- left_join(trend, sst)
sst_boreal_brm <- readRDS("./output/sst_boreal_brm.rds")

# plot
## 95% CI
ce1s_1 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sst_boreal_brm, effect = "sst.anomaly", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$sst.anomaly
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$sst.anomaly[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$sst.anomaly[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$sst.anomaly[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$sst.anomaly[["lower__"]]


fig2a <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  labs(x = "SST anomaly wrt 1854-1949 (SD)", y = "Borealization index") +
  geom_text(data = dat, aes(sst.anomaly, trend, label = year), size = 3)


probs <- read.csv("./output/probabilistic_attribution_stats.csv")

fig2b <- ggplot(probs, aes(year, FAR)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = LCI_FAR,
                  ymax = UCI_FAR),
              fill = "dark grey", 
              alpha = 0.5) +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ylab("Fraction of attributable risk")

fig2c <- ggplot(probs, aes(year, RR)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab("Risk ratio")

## pdfs
resample.pdf <- read.csv("./output/resampled_borealization_pdfs.csv")

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)
 
resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))


fig2d <- ggplot(resample.pdf, aes(period, borealization_index)) +
  geom_hline(yintercept = c(-1,2), color = cb[c(6,7)], lty = 2, alpha = 0.8) +
  geom_violin(fill = cb[2], lty = 0, alpha = 0.5) +
  coord_flip(ylim = c(-4.5, 4.5)) +
  xlab("North Pacific warming") +
  ylab("Borealization index") +
  scale_x_discrete(labels = c("Preindustrial",
                              "1950 to 0.5°",
                              "0.5° to 1.0°",
                              "1.0° to 1.5°",
                              "1.5° to 2.0°")) +
  scale_y_continuous(breaks = seq(-4,4, by = 2))


fig2d

png("./figs/fig2.png", width = 8, height = 7, units = 'in', res = 300)

ggpubr::ggarrange(ggpubr::ggarrange(fig2a, fig2b, fig2c,  ncol = 1, labels = "auto"),
                  fig2d, ncol = 2, widths = c(0.45, 0.55), labels = c("", "d"))

dev.off()

## version without sst-borealization plot
png("./figs/fig2_three_panel.png", width = 7, height = 5, units = 'in', res = 300)

ggpubr::ggarrange(ggpubr::ggarrange(fig2b, fig2c,  ncol = 1, labels = "auto"),
                  fig2d, ncol = 2, widths = c(0.5, 0.5), labels = c("", "c"))

dev.off()
