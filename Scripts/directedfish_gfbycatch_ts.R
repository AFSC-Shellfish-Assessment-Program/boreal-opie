library(tidyverse)
library(raster)
library(terra)

# DIRECTED FISHERY TS ---------------------------------------------------------------------------------------------
#1995-2020 directed fishery raw and potsum files
raw_df_9520 <- read.csv("./Data/SNOWCRAB-1995-2020_crab_dump.csv")
potsum_df_9520 <- read.csv("./Data/SNOWCRAB-1995-2020_potsum.csv")

#2021 directed fishery raw and potsum files
raw_df_21 <- read.csv("./Data/SNOWCRAB-2021_crab_dump.csv")
potsum_df_21 <- read.csv("./Data/SNOWCRAB-2021_potsum.csv")

#View(potsum_df_9520)


names(raw_df_9520)

raw_df <- rbind(raw_df_9520, raw_df_21)
potsum_df <- rbind(potsum_df_9520, potsum_df_21)

#write.csv(raw_df, "./Data/raw_df.csv")

raw_df %>%
  mutate(sampdate = lubridate::mdy(sampdate), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate)) %>%
  group_by(trip, adfg, spn, year) %>%
  summarise(female = sum(sex == 2),
            imm_male = sum(sex == 1 & size <95),
            male_3059 = sum(sex == 1 & size <60 & size >=30),
            male_6095 = sum(sex == 1 & size >=60 & size <=95),
            mat_male = sum(sex == 1 & size >=95),
            legal = sum(sex == 1 & size >=78)) %>%
  group_by(trip, adfg, year) %>%
  summarise(female = sum(female),
            imm_male = sum(imm_male),
            mat_male = sum(mat_male),
            male_3059 = sum(male_3059),
            male_6095 = sum(male_6095),
            legal = sum(legal),
            msr_potlifts = n(),
            male_3059pp = male_3059/msr_potlifts,
            male_6095pp = male_6095/msr_potlifts) -> male_sum

potsum_df  %>%
  mutate(sampdate = lubridate::mdy(sampdate), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate)) %>%
  group_by(trip, adfg, year) %>%
  summarise(msr_potlifts = sum(msr_pot == "Y"),
            tot_potlifts = max(spn)) -> total_plfs  

right_join(total_plfs, male_sum, by=c("trip", "adfg", "year", "msr_potlifts")) %>%
  group_by(year) %>%
  summarise(msr_potlifts = sum(msr_potlifts),
            tot_imm_male = sum(imm_male),
            tot_mat_male = sum(mat_male),
            tot_male3059 = sum(male_3059),
            tot_male6095 = sum(male_6095),
            logmean_male3059pp = log(mean(male_3059pp + 1)),
            logmean_male6095pp = log(mean(male_6095pp + 1)),
            imm_male_pp = tot_imm_male/msr_potlifts,
            mat_male_pp = tot_mat_male/msr_potlifts) -> opilio_df_males_pp

write.csv(opilio_df_males_pp, "./Data/opilio_male_directfish_ts.csv")



ggplot() +
  geom_line(opilio_df_males_pp, mapping = aes(x = year, y = imm_male_pp, col = "red"), linewidth = 1)+
  geom_line(opilio_df_males_pp, mapping = aes(x = year, y = mat_male_pp, col = "blue"), linewidth = 1) +
  ylab("Total per Pot") +
  xlab("Year")+
  scale_x_continuous(breaks= seq(min(opilio_df_males_pp$year), max(opilio_df_males_pp$year), by=5),minor_breaks = NULL)+
  scale_color_manual(name = NULL, values = c("blue", "red"),
                    labels = c("Mature (>=95mm)", "Immature (<95mm)"))+
  theme_bw()


## GF BYCATCH -------------------------------------------------------------------------------------------------------------------
# Read in data
length_opilio <- read.csv("./Data/DEBRIEFED_CRAB_LENGTH.txt") %>%
  filter(SPECIES_NAME == "OPILIO TANNER CRAB")
    
haul_opilio <- read.csv("./Data/DEBRIEFED_HAUL.txt") %>%
  filter(SPECIES_NAME == "OPILIO TANNER CRAB")
    
# Pull out months and year
haul_opilio %>%
  mutate(sampdate = lubridate::mdy(HAUL_DATE), 
         year = lubridate::year(sampdate),
        month = lubridate::month(sampdate),
        day = lubridate::day(sampdate)) -> haul_opilio
  #filter(month %in% 6:8, GEAR_TYPE == "NON_PELAGIC_TRAWL") -> summer_npt # unhash if you just want NPT data in summer
    
  
#Calculate total crab measured, total immature, and total mature males and the proportion of each by haul in core area
length_opilio %>%
  filter(HAUL_JOIN %in% haul_opilio$HAUL_JOIN) %>%
  group_by(HAUL_JOIN) %>%
  summarise(msr_crab = sum(FREQUENCY),
            imm_total = sum(SEX == "M" & LENGTH < 95),
            mat_total = sum(SEX == "M" & LENGTH >=95),
            male_3059 = sum(SEX == "M" & LENGTH >=30 & LENGTH <60),
            male_6095 = sum(SEX == "M" & LENGTH >=60 & LENGTH <=95),
            imm_prop = imm_total/msr_crab,
            mat_prop = mat_total/msr_crab,
            male_3059_prop = male_3059/msr_crab,
            male_6095_prop = male_6095/msr_crab) -> length_sum
  
  #Calculate total number of crab extrapolated up to the haul by haul join, haul, and year
  haul_opilio %>%
    group_by(HAUL_JOIN, VESSEL_ID, HAUL_SET_NUMBER, year) %>%
    summarise(catch_extrap = sum(NUMBER_OF_SPECIES_EXTRAPOLATED_UP_TO_THE_HAUL)) -> haul_sum
  
  #Join haul data with length data
  right_join(haul_sum, length_sum, by = "HAUL_JOIN") %>%
    summarise(year = year,
              catch_extrap = catch_extrap,
              imm_extrap = imm_prop * catch_extrap,
              mat_extrap = mat_prop * catch_extrap,
              male_3059_extrap = male_3059_prop * catch_extrap,
              male_6095_extrap = male_6095_prop * catch_extrap) %>%
    group_by(year) %>%
    summarise(catch_total = sum(catch_extrap),
              imm_total = sum(imm_extrap),
              mat_total = sum(mat_extrap),
              logmean_male3059 = log(mean(male_3059_extrap + 1)),
              logmean_male6095 = log(mean(male_6095_extrap + 1))) -> opilio_bycatch_year_total
  
  write.csv(opilio_bycatch_year_total, "./Data/opilio_male_bycatch_ts.csv")
  
  #Plot by year
  ggplot() +
    geom_line(bycatch_year_total, mapping = aes(x = year, y = imm_total, col = "red"), size = 1)+
    geom_line(bycatch_year_total, mapping = aes(x = year, y = mat_total, col = "blue"), size = 1) +
    ylab("Total catch per year") +
    xlab("Year")+
    scale_x_continuous(breaks= seq(min(bycatch_year_total$year), max(bycatch_year_total$year), by=3),minor_breaks = NULL)+
    scale_color_manual(name = NULL, values = c("blue", "red"),
                       labels = c("Mature (>=95mm)", "Immature (<95mm)"))+
    theme_bw() -> total_catch_per_year

# save just bycatch / directed catch data for imputation -------------------------
  opie_bycatch <- read.csv("./Data/opilio_male_bycatch_ts.csv", row.names = 1) %>%
    dplyr::select(year, logmean_male3059, logmean_male6095) %>%
    filter(year >= 1988) # truncating first two years with very low catch
  
  opie_directfish <- read.csv("./Data/opilio_male_directfish_ts.csv", row.names = 1) %>%
    dplyr::select(year, logmean_male3059pp, logmean_male6095pp)
  
# combine
out <- left_join(opie_bycatch, opie_directfish) %>%
  dplyr::rename(small_log_cpue_bycatch = logmean_male3059,
         large_log_cpue_bycatch = logmean_male6095,
         small_log_cpue_directed = logmean_male3059pp,
         large_log_cpue_directed = logmean_male6095pp)

write.csv(out, "./output/directed_bycatch_all_years.csv")
    
# CROSS-CORRELATION BETWEEN DIRECTED FISHERY and BYCATCH DATA for OPILIO and YFS, POLLOCK, and PLAICE RECRUIT AND SSB ---------------------------------------------------------------------------------------------

  years <- 1995:2019 #set years for data
    
  #Read in data, filter by years
    opie_bycatch <- read.csv("./Data/opilio_male_bycatch_ts.csv") %>%
      filter(year %in% years) #opie bycatch
    
    opie_directfish <- read.csv("./Data/opilio_male_directfish_ts.csv") %>%
      filter(year %in% years) #opie directed fishery catch
    
    yfs_ssbrecruit <- read.csv("./Data/BSAIyfin_SSBrecruit.csv") %>%
      filter(year %in% years) #yfs recruitment and SSB
             
    pollock_ssbrecruit <- read.csv("./Data/BSpollock_SSBrecruit.csv") %>%
      filter(year %in% years) #pollock recruitment and SSB
    
    plaice_ssbrecruit <- read.csv("./Data/BSplaice_SSBrecruit.csv") %>%
      filter(year %in% years) #plaice recruitment and SSB
    
  # Log transform ssb and recruitment timeseries for groundfish
    yfs_ssbrecruit %>%
      mutate(logSSB_yfs = log(SSB + 1),
             logrecruit_yfs = log(recruitment +1)) -> yfs_ssbrecruit
    
    pollock_ssbrecruit %>%
      mutate(logSSB_pol = log(SSB + 1),
             logrecruit_pol = log(recruitment +1)) -> pollock_ssbrecruit
    
    plaice_ssbrecruit %>%
      mutate(logSSB_plaice = log(SSB +1),
             logrecruit_plaice = log(recruitment +1)) -> plaice_ssbrecruit
    
  # Create one dataframe with opie directed fishery catch per pot (both sizes), opie bycatch (both sizes), yfs recruit and SSB, 
  # plaice recruit and SSB, and pollock recruit and SSB, all log transformed
  comb_df <- data.frame(year = years,
                   male3059.directfish = opie_directfish$logmean_male3059, male3059.bycatch = opie_bycatch$logmean_male3059,
                   male6095.directfish = opie_directfish$logmean_male6095, male6095.bycatch = opie_bycatch$logmean_male6095,
                   yfs.recruit = yfs_ssbrecruit$logrecruit_yfs, yfs.ssb = yfs_ssbrecruit$logSSB_yfs,
                   pollock.recruit = pollock_ssbrecruit$logrecruit_pol, pollock.ssb = pollock_ssbrecruit$logSSB_pol,
                   plaice.recruit = plaice_ssbrecruit$logrecruit_plaice, plaice.ssb = plaice_ssbrecruit$logSSB_plaice)
  
  # Create a dataframe of all possible combinations of variable names (different timeseries)
  var_names <- expand.grid(colnames(comb_df[-1]), colnames(comb_df[-1]))
  
  # Create empty dataframe to store ccf output by pairwise combination
  ccf_df <- as.data.frame(matrix(ncol = nrow(var_names), nrow = 21))
  
  # Run for loop to calculate ccf by pairwise combination, store in df
  for(ii in 1:nrow(var_names)){
    var1 = var_names[ii, 1]
    var2 = var_names[ii, 2]
    
    ccf(comb_df[grepl(var1, colnames(comb_df))], comb_df[grepl(var2, colnames(comb_df))], 
                      plot = FALSE)$acf[,,1] -> ccf_df[,ii]
    
    colnames(ccf_df)[ii] <- paste(var1, var2, sep = "-")
    
  }

  # Save output to csv
  write.csv(ccf_df, "./Output/ccf_opiliogroundfish.csv")
  write.csv(comb_df, "./Output/bycatch_directfish_groundfish_ts.csv")
