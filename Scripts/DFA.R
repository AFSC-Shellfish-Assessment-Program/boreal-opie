# combine multiple indicators of borealization with 
# a Dynamic Factor Analysis model

library(tidyverse)
library(MARSS)
library(corrplot)
library(oce)
theme_set(theme_bw())

# Updated 1/10/23 with 2022 data by EJF

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## data processing -------------------------------------

#bloom timing
d1 <- read.csv("./Data/bloom_timing.csv")

d1 <- d1 %>%
  filter(north_south == "south") %>%
  rename(value = globcolour_peak_mean) %>%
  mutate(name = "Bloom timing") %>%
  dplyr::select(year, name, value)

#bloom type  
d2 <- read.csv("./Data/bloom_type.csv")

check  <- d2 %>%
  filter(north_south == "south") %>%
  mutate(name = case_when(gl_type == "ice_full" ~ "Ice-edge bloom",
                          gl_type == "ice_free" ~ "Open water bloom")) %>%
  rename(value = count) %>%
  dplyr::select(year, name, value) %>%
  group_by(name) %>%
  summarise(count = n())

check

d2  <- d2 %>%
  filter(north_south == "south") %>%
  mutate(name = case_when(gl_type == "ice_full" ~ "Ice-edge bloom",
                          gl_type == "ice_free" ~ "Open water bloom")) %>%
  rename(value = count) %>%
  dplyr::select(year, name, value)%>%
  filter(name == "Open water bloom")

#sea ice
d3 <- read.csv("./Data/ice.csv")

d3 <- d3 %>%
  rename(`Jan-Feb ice cover` = JanFeb_ice,
         `Mar-Apr ice cover` = MarApr_ice) %>%
  filter(year >= 1972) %>%
  pivot_longer(cols = -year)

#chl-a size fraction
d4 <- read.csv("./Data/summarized_chla.csv")

d4 <- d4 %>%
  rename(`Bloom size` = total_avg_chla,
         `Phytoplankton size ratio` = avg_chla_ratio) %>%
  pivot_longer(c(1,3)) 

#BCS prevalence 
d5 <- read.csv("./Data/bcs_prev.csv", row.names = 1)

d5 <- d5 %>% 
  dplyr::select(YEAR, Population) %>%
  rename(value = Population,
         year = YEAR) %>%
  mutate(name = "Hematodinium")
  
#bottom temp
d6 <- read.csv("./Data/date_corrected_bottom_temp.csv")

d6 <- d6 %>%
  rename(value = bottom.temp) %>%
  mutate(name = 'Bottom temperature')

#groundfish CPUE
d7 <- read.csv("./Data/groundfish_mean_cpue.csv", row.names = 1)

d7 <- d7 %>%
  rename(year = YEAR,
         `Pacific cod` = mean_cod_CPUE,
         `Arctic groundfish` = mean_arctic_CPUE) %>%
  pivot_longer(cols = -year)

#zooplankton abundances
d8 <- read.csv("./Data/summarized_zooplankton.csv")

d8 <- d8 %>%
  rename(name = taxa, 
         value = log_abundance)

dat <- rbind(d1, d2, d3, d4, d5, d6, d7, d8)

# add NAs for plot
plot.dat <- dat %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  pivot_longer(cols = -year)

# change Calanus for plotting
change <- plot.dat$name == "Calanus_glacialis"
plot.dat$name[change] <- "Calanus"

# reorder for plot
plot.order <- data.frame(name = unique(plot.dat$name),
                         order = c(5,4,1,2,6,7,10,3,11,12,8,9))

plot.dat <- left_join(plot.dat, plot.order)

plot.dat$name <- reorder(plot.dat$name, plot.dat$order)


# save as Extended Data Fig. 4
tiff("./Figs/Extended_Data_Fig_4.tiff", width = 10, height = 6, units = 'in', res = 300)

ggplot(plot.dat, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  theme(axis.title.x = element_blank()) +
  ylab("Value")

dev.off()




# save time series
write.csv(dat, "./output/dfa time series.csv", row.names = F)

dfa.dat <- plot.dat %>%
  dplyr::select(-order) %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(year) %>%
  dplyr::select(-year) %>%
  t()
  
colnames(dfa.dat) <- unique(d3$year)

# and plot correlations
cors <- cor(t(dfa.dat), use = "p")
diag(cors) <- 0

max(cors)
min(cors) 

plot <- as.data.frame(t(dfa.dat))

# plot correlations
corrplot(cors, method = "sq", col.lim = c(-0.865, 0.865), col = oceColorsPalette(64), tl.col = "black", cl.cex = 0.7, order = "FPC")


# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models & store results
for(R in levels.R) {
  for(m in 1:2) {  # considering either 1- or 2-trend model
    
    dfa.model = list(A="zero", R=R, m=m)
    
    kemz = MARSS(dfa.dat, model=dfa.model,
                 form="dfa", z.score=TRUE, control=cntl.list)
    
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data # diagonal and unequal is the best model but doesn't converge

# save model selection table--note that unconstrained models did not converge!
write.csv(model.data, "./output/dfa_model_selection_table.csv",
          row.names = F)

## fit the best model --------------------------------------------------
model.list = list(A="zero", m=2, R="diagonal and unequal") # best model is two-trend

# not sure that these changes to control list are needed for this best model, but using them again!
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# rotate
# get the inverse of the rotation matrix
Z.est <- coef(mod, type = "matrix")$Z

H.inv <- varimax(coef(mod, type = "matrix")$Z)$rotmat

# rotate factor loadings
Z.rot <- Z.est %*% H.inv

# rotate trends
trends.rot <- solve(H.inv) %*% mod$states

# Add CIs to marssMLE object
mod <- MARSSparamCIs(mod)

# Use coef() to get the upper and lower CIs
Z.low <- coef(mod, type = "Z", what = "par.lowCI")
Z.up <- coef(mod, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv

plot.CI <- data.frame(names=rownames(dfa.dat),
  mean = as.vector(Z.rot),
  upCI = as.vector(Z.rot.up),
  lowCI = as.vector(Z.rot.low)
)
plot.CI

dodge <- position_dodge(width=0.9)

plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)


plot.CI$trend <- rep(c("T1", "T2"), each = 12)

loadings.plot <- ggplot(plot.CI, aes(x=names, y=mean, fill = trend)) +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=9), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

# plot trend
trend <- data.frame(trend = rep(c("T1", "T2"), each = length(1972:2022)),
                    t=1972:2022,
                    estimate=as.vector(mod$states),
                    conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                    conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))


trend.plot <- ggplot(trend, aes(t, estimate, color = trend, fill = trend)) +
  theme_bw() +
  geom_line() +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=0, alpha=0.1) + xlab("") + ylab("Trend")

# save
png("./Figs/best_two_trend_DFA_loadings_trend.png", width = 9, height = 3.5, units = 'in', res = 300)

ggpubr::ggarrange(loadings.plot,
                  trend.plot,
                  ncol = 2,
                  widths = c(0.45, 0.55),
                  labels = "auto")

dev.off()
# only one loading can be distinguished from 0! 
# reject this model and fit second-best model (1 trend diagonal and unequal)

model.list = list(A="zero", m=1, R="diagonal and unequal") # second-best model - this is the borealization index

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# save 
saveRDS(mod, "./output/DFA_model.rds")

# plot fits to data
DFA_pred <- print(predict(mod))

DFA_pred <- DFA_pred %>%
  mutate(year = rep(1972:2022, 12)) 

# get R^2 for each time series
summarise <- DFA_pred %>%
  group_by(.rownames) %>%
  summarise(R_sq = cor(y, estimate, use = "pairwise")^2) %>%
  mutate(plot_label = paste(.rownames, " (", round(R_sq, 3), ")", sep = ""))

DFA_pred <- left_join(DFA_pred, summarise)

tiff("./figs/Extended_Data_Fig_5.tiff", width = 10, height = 6, units = 'in', res = 300)

ggplot(DFA_pred, aes(estimate, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~plot_label, ncol = 4, scale = "free") +
  labs(x = "Estimated", y = "Observed")

dev.off()

# process loadings and trend

CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                          mean=CI$par$Z[1:12],
                          upCI=CI$par.upCI$Z[1:12],
                          lowCI=CI$par.lowCI$Z[1:12])

dodge <- position_dodge(width=0.9)


plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:12])


# plot trend
trend <- data.frame(t=1972:2022,
                        estimate=as.vector(mod$states),
                        conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                        conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))



# and save loadings and trend
write.csv(plot.CI, "./Output/dfa_loadings.csv", row.names = F)
write.csv(trend, "./Output/dfa_trend.csv", row.names = F)
