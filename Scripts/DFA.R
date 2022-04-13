# combine multiple indicators of borealization with 
# a Dynamic Factor Analysis model

library(tidyverse)
library(MARSS)

theme_set(theme_bw())

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## data processing -------------------------------------

d1 <- read.csv("./Data/bloom_timing.csv")

d1 <- d1 %>%
  mutate(name = paste(sub_6domain, "_bloom_timing", sep = "")) %>%
  rename(value = peak_mean) %>%
  select(year, name, value)
  
d2 <- read.csv("./Data/bloom_type.csv")

d2 <- d2 %>%
  filter(!is.na(bloom_type))

ggplot(d2, aes(year, count, color = bloom_type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~sub_6domain)

d2 <- d2 %>%
  mutate(name = paste(sub_6domain, bloom_type, "bloom", sep = "_")) %>%
  rename(value = count) %>%
  select(year, name, value)

d3 <- read.csv("./Data/march_ice_cover.csv")
  
d3 <- d3 %>%
  mutate(name = "March_ice_cover") %>%
  rename(value = sic_average)

d4 <- read.csv("./Data/chl_a.csv")

d4 <- d4 %>%
  pivot_longer(cols = -year)

dat <- rbind(d1, d2, d3, d4)

ggplot(dat, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y")

dfa.dat <- dat %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(desc(year))
  select(-year) %>%
  t()
  
colnames(dfa.dat) <- d3$year

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
  for(m in 1:3) {  # allowing up to 3 trends
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
model.data

# save model selection table
write.csv(model.data, "./Results/legacy catch dfa model selection table 1956-1990.csv",
          row.names = F)

## fit the best model --------------------------------------------------
model.list = list(A="zero", m=2, R="diagonal and unequal") # best model for early era

# not sure that these changes to control list are needed for this best model, but using them again!
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# and rotate the loadings
Z.est = coef(mod, type="matrix")$Z
H.inv = varimax(coef(mod, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)


Z.rot$names <- rownames(dfa.dat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(dfa.dat)
Z.rot$plot.names <- reorder(Z.rot$names, 1:length(Z.rot$names))

loadings <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  ylab("Loading") + xlab("") + ggtitle("1956-1990 legacy catch - diagonal and unequal") + 
  scale_fill_manual(values=cb[2:3]) +
  theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) 

loadings


# plot trends
# first rotate the trends!
H.inv = varimax(coef(mod, type="matrix")$Z)$rotmat
trends.rot = solve(H.inv) %*% mod$states

trends.plot <- data.frame(year = 1956:1990,
                          trend1 = trends.rot[1,],
                          trend2 = trends.rot[2,]) %>%
  pivot_longer(cols = -year)

trend <- ggplot(trends.plot, aes(year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[2:3]) +
  ggtitle("1956-1990 legacy catch - diagonal and unequal") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8)) +
  geom_hline(yintercept = 0)

trend

# combine in one saved plot
png("./Figs/legacy_dfa_loadings_and_trend.png", 
    width=5, height=8, units='in', res=300)

ggpubr::ggarrange(loadings, trend, nrow = 2, ncol = 1,
                  heights = c(0.6, 0.4))

dev.off()

ggsave("./Figs/legacy_catch_trend_plot_diagonal_unequal.png", width = 5.5, height = 3, units = 'in')
