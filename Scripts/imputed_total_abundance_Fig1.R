# total abundance for Fig. 1
# impute for missing stations

# load ----

library(tidyverse)
library(Hmisc)
library(patchwork)
library(moments)

theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##############################################

## data
#EBS haul data ----
sc_catch <- read.csv("./Data/crabhaul_opilio.csv")

#EBS strata data ----
sc_strata <- read_csv("./Data/STRATA_OPILIO_NEWTIMESERIES.csv")

#Standard 375 survey stations
sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  pull(STATION_ID) ->standard

#################################
#################################
# process catch data
sc_catch %>%
  mutate(YEAR = floor(CRUISE/100)) %>%
  filter(HAUL_TYPE %in% c(3,4)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(N_CRAB = sum(SAMPLING_FACTOR, na.rm = T),
            CPUE = N_CRAB / mean(AREA_SWEPT)) %>%
  #add zero catch stations
  right_join(sc_catch %>% 
               mutate(YEAR = floor(CRUISE/100)) %>%
               filter(HAUL_TYPE %in% c(3,4)) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)) %>%
  replace_na(list(CPUE = 0)) %>%
  mutate(log_cpue = log(CPUE +1)) %>%
  #filter for standard stations (1979 has 70 extra non-standard stations)
  filter(GIS_STATION %in% standard) -> cpue


# Remove stations that have never caught snow crab 
cpue %>%
  group_by(GIS_STATION) %>%
  summarise(station_mean = mean(CPUE)) %>% 
  filter(station_mean > 0) %>% #58 stations dropped
  #b) Remove stations < 5th CPUE percentile range 
  filter(station_mean > quantile(station_mean, 0.05)) %>% #16 additional stations dropped
  pull(GIS_STATION) -> keepers 

#Final dataset for imputation
cpue %>%
  filter(GIS_STATION %in% keepers) -> dat

# get data set of excluded stations for later averaging when computing total abundance
cpue %>%
  filter(!GIS_STATION %in% keepers) -> excluded_stations

excluded_stations <- unique(excluded_stations$GIS_STATION)

# plot means by station
plot_dat <- dat %>%
  group_by(YEAR) %>% 
  mutate(mean_log_cpue = mean(log_cpue))


ggplot(plot_dat, aes(mean_log_cpue)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black")


# plot means on map
plot_map <- dat %>%
  group_by(GIS_STATION) %>% 
  mutate(mean_log_cpue = mean(log_cpue),
         latitude = mean(MID_LATITUDE),
         longitude = mean(MID_LONGITUDE))

abundance_map <- ggplot(plot_map, aes(longitude, latitude, color = mean_log_cpue)) +
  geom_point(size = 3, shape = 15) +
  scale_colour_gradient(
    low = "purple",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  theme(legend.position = c(0.8, 0.8))

abundance_map 

# get into matrix form for mice
dat <- tapply(dat$log_cpue, list(dat$YEAR, dat$GIS_STATION), mean)


# examine correlations
dat1 <- dat2 <- as.matrix(dat)

r <- cor(dat1, dat2, use = "p")
r #Cross-year correlations between each station combination

# choose 15 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

# plot r used in imputation
plot_r <- data.frame(r = r[pred])

all_drop5_hist <- ggplot(plot_r, aes(r)) +
  geom_histogram(fill = "grey", color = "black", bins = 50)

all_drop5_hist

colnames(pred) <- rownames(pred) <- colnames(dat) <- str_remove_all(colnames(pred), "-")

blocks <- mice::make.blocks(dat)

imp <- mice::mice(data = dat, method = "pmm", m=100, predictorMatrix = pred, blocks = blocks)

saveRDS(imp, "./output/total_abundance_imputations.RDS")
imp_total_abundance_drop_5 <- readRDS("./output/total_abundance_imputations.RDS")


# are there NAs in complete(imp)?
check <- is.na(complete(imp_total_abundance_drop_5))

sum(check) # 0

# get weights (area of each stratum) for weighted mean abundance
# first, strata areas for weighted mean abundance - use 2022, when all stations were sampled

area <- sc_strata %>%
  filter(SURVEY_YEAR == 2022) %>%
  #STATION_ID %in% unique(drop_10th$GIS_STATION)) %>%
  mutate(STATION = str_remove_all(STATION_ID, "-")) %>%
  select(STATION, DISTRICT, TOTAL_AREA)

# sum_area <- area %>%
#   group_by(DISTRICT) %>%
#   summarise(total_count = n(),
#             total_area = mean(TOTAL_AREA))
# 
# weights <- data.frame(STATION = colnames(complete(imp_total_abundance_drop_5, action = 1))) %>%
#   left_join(., area)
# 
# # get proportion of each stratum that is included in the imputation
# summary_weights <- weights %>%
#   group_by(DISTRICT) %>%
#   summarise(count_sampled = n()) %>%
#   left_join(., sum_area) %>%
#   mutate(weight = (count_sampled / total_count) * total_area) %>%
#   select(DISTRICT, weight)
# 
# weights <- left_join(weights, summary_weights)

# set up weighted mean function
ff <- function(x) weighted.mean(x, weights$weight, na.rm = T)

imputed.dat <- imputed.back.trans.dat <- concentration.dat <-  data.frame()

# create bounds for imputed values 
# (min = 0, max = max observed)

lower_bound <- 0
upper_bound <- max(dat$log_cpue, na.rm = T) 

# data frame of 0 cpue for excluded stations
excluded <- data.frame(STATION = str_remove_all(excluded_stations, "-"),
                       cpue = 0)


imputed_abundance <- data.frame()

for(i in 1:100){
  
  # i <- 1
  temp <- complete(imp_total_abundance_drop_5, action = i)
  
  # apply bounds
  change <- temp < lower_bound
  temp[change] <- lower_bound
  
  change <- temp > upper_bound
  temp[change] <- upper_bound
                              
  
  # and get back-transformed means
  back.temp <- exp(temp)-1


  # loop through each year and calculate total estimated abundance
  # (sum of mean cpue*stratum areas)
  for(j in 1:nrow(back.temp)){
    # j <- 2
    
    temp.temp <- data.frame(STATION = colnames(back.temp),
                            cpue = as.vector(t(back.temp[j,])))
    
    # add in excluded stations as cpue = 0
    temp.temp <- rbind(temp.temp,
                       excluded)
    
    # join with area for each stratum
    temp.temp <- left_join(temp.temp, area)
    
    
    temp.sum <- temp.temp %>%
      group_by(DISTRICT) %>%
      summarise(mean_cpue = mean(cpue),
                total_area = mean(TOTAL_AREA),
                abundance = mean_cpue*total_area)
    
    imputed_abundance <- rbind(imputed_abundance,
                               data_frame(imputation = i,
                                          year = rownames(back.temp)[j],
                                          abundance = sum(temp.sum$abundance)))
    
  } # close j loop (years/rows)
  
print(paste("Imputation # ", i, sep = ""))
  
} # close i loop (imputation)


# summarize and plot
imputation_summary <- imputed_abundance %>%
  group_by(year) %>%
  summarise(abundance = mean(abundance)) %>%
  mutate(SD = tapply(imputed_abundance$abundance, imputed_abundance$year, sd))
          
# add in NAs for 2020

xtra <- data.frame(year = 2020,
                   abundance = NA,
                   SD = NA)

imputation_summary <- rbind(imputation_summary, xtra)

# and divide by 1e9 (billions of animals!)
imputation_plot <- imputation_summary %>%
  mutate(abundance = abundance / 1e9,
         SD = SD / 1e9)


ggplot(imputation_plot, aes(as.numeric(year), abundance)) +
  geom_col(fill = cb[6]) + 
  geom_errorbar(aes(ymin = abundance - 2*SD,
                    ymax = abundance + 2*SD))

# save data for combined Fig 1
write.csv(imputation_plot, "./output/imputed_total_abundance.csv", row.names = F)

# first, back-transformed data
plot.back.trans <- data.frame(year = c(1975:2019, 2021:2022),
                              mean = tapply(log(imputed.back.trans.dat$mean), imputed.back.trans.dat$year, mean),
                              sd = tapply(log(imputed.back.trans.dat$mean), imputed.back.trans.dat$year, sd))

# add in NAs
xtra <- data.frame(year = 2020,
                   mean = NA,
                   sd = NA)

plot.back.trans <- rbind(plot.back.trans,
                         xtra)


# replace sd = 0 with NA
check <- plot.back.trans$sd == 0
plot.back.trans$sd[check] <- NA

# save
write.csv(plot.back.trans, "imputed_total_abundance.csv", row.names = F)

# plot
ggplot(plot.back.trans, aes(year, mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd))