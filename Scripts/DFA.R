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

d5 <- read.csv("./Data/bcs_prev.csv", row.names = 1)

plot <- d5 %>%
  pivot_longer(cols = -year)

ggplot(plot, aes(year, value, color = name)) +
  geom_line() +
  geom_point()


d5$immature_bcs <- rowMeans(d5[,2:3], na.rm = T)

d5 <- d5 %>% 
  select(year, immature_bcs) %>%
  pivot_longer(cols = -year)
  

dat <- rbind(d1, d2, d3, d4, d5)

ggplot(dat, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y")

dfa.dat <- dat %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(year) %>%
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
  for(m in 1:1) {  # find best single-trend model
    
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
write.csv(model.data, "./output/dfa_model_selection_table.csv",
          row.names = F)

## fit the best model --------------------------------------------------
model.list = list(A="zero", m=2, R="diagonal and equal") # best model 

# not sure that these changes to control list are needed for this best model, but using them again!
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# plot
CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                          mean=CI$par$Z[1:12],
                          upCI=CI$par.upCI$Z[1:12],
                          lowCI=CI$par.lowCI$Z[1:12])

dodge <- position_dodge(width=0.9)


plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:11])

loadings.plot <- ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)
