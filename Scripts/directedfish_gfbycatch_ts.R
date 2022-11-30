#1995-2020 directed fishery raw and potsum files
raw_df_9520 <- read.csv("./Data/SNOWCRAB-1995-2020_crab_dump.csv")
potsum_df_9520 <- read.csv("./Data/SNOWCRAB-1995-2020_potsum.csv")

#2021 directed fishery raw and potsum files
raw_df_21 <- read.csv("./Data/SNOWCRAB-2021_crab_dump.csv")
potsum_df_21 <- read.csv("./Data/SNOWCRAB-2021_potsum.csv")

View(potsum_df_9520)


names(raw_df_9520)

raw_df <- rbind(raw_df_9520, raw_df_21)
potsum_df <- rbind(potsum_df_9520, potsum_df_21)

write.csv(raw_df, "./Data/raw_df.csv")

raw_df %>%
  mutate(sampdate = lubridate::mdy(sampdate), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate)) %>%
  group_by(trip, adfg, spn, year) %>%
  summarise(female = sum(sex == 2),
            imm_male = sum(sex == 1 & size <95),
            mat_male = sum(sex == 1 & size >=95),
            legal = sum(sex == 1 & size >=78)) %>%
  group_by(trip, adfg, year) %>%
  summarise(female = sum(female),
            imm_male = sum(imm_male),
            mat_male = sum(mat_male),
            legal = sum(legal),
            msr_potlifts = n()) -> male_sum

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
            imm_male_pp = tot_imm_male/msr_potlifts,
            mat_male_pp = tot_mat_male/msr_potlifts) -> opilio_df_males_pp


ggplot() +
  geom_line(opilio_df_males_pp, mapping = aes(x = year, y = imm_male_pp, col = "red"), size = 1)+
  geom_line(opilio_df_males_pp, mapping = aes(x = year, y = mat_male_pp, col = "blue"), size = 1) +
  ylab("Total per Pot") +
  xlab("Year")+
  scale_x_continuous(breaks= seq(min(opilio_df_males_pp$year), max(opilio_df_males_pp$year), by=5),minor_breaks = NULL)+
  scale_color_manual(name = NULL, values = c("blue", "red"),
                    labels = c("Mature (>=95mm)", "Immature (<95mm)"))+
  theme_bw()



raw_df %>%
  mutate(sampdate = lubridate::mdy(sampdate), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         matsex = case_when((sex == 1 & size < 95) ~ "Immature Male",
                           (sex == 1 & size >= 95) ~ "Mature Male",
                           (sex == 2) ~ "Female")) %>%
  group_by(trip, adfg, year, matsex) %>%
  summarise(total = n())




  summarise(female = sum(sex == 2),
            imm_male = sum(sex == 1 & size <95),
            mat_male = sum(sex == 1 & size >=95),
            legal = sum(sex == 1 & size >=78)) %>%
  group_by(trip, adfg, year) %>%
  summarise(female = sum(female),
            imm_male = sum(imm_male),
            mat_male = sum(mat_male),
            legal = sum(legal),
            msr_potlifts = n()) -> male_sum
  
###################################################################################################################### 
## Snow gf bycatch
  
#---------------------
#Filtering data by core area stations
  #Define core stations
  core_sta <- c("P-32", "Q-31", "R-31", "S-31", "T-30", "U-29", "V-28", "V-27", "V-26", "V-25", "U-25", "T-25", "S-25",
                "S-24", "S-23", "S-22", "R-22", "Q-22", "Q-21", "Q-20", "Q-19", "Q-18", "P-18", "O-18", "N-18", "N-01",
                "M-01", "L-01", "K-01", "K-02", "L-03", "K-04", "K-05", "K-06", "J-06", "H-06", "H-05", "H-04", "H-03",
                "H-02", "G-02", "G-01", "E-18", "F-19", "GF_2019", "GF-2120", "H-22", "I-24", "J-25", "L-27", "N-28", 
                "O-29", "O-30", "P-30", "P-31")
    
  #Filter GIS stations in core area df by core stations, pull out core lat/lon 
  perc50_core %>%
    filter(GIS_STATION %in% core_sta) -> xx
    
  xx[match(core_sta, xx$GIS_STATION), ] %>%
    select(LATITUDE, LONGITUDE) %>%
    rename(y = LATITUDE, x = LONGITUDE) %>%
    na.omit()-> core_crds
  
  #Create polygon around core stations
  library(sf)
  core_crds %>%
    st_as_sf(coords = c("x", "y"), crs = 4326) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON") %>%
    vect() -> core_poly
  
  
  length_opilio <- read.csv("./Data/DEBRIEFED_CRAB_LENGTH.txt") %>%
    filter(SPECIES_NAME == "OPILIO TANNER CRAB")
  
  haul_opilio <- read.csv("./Data/DEBRIEFED_HAUL.txt") %>%
    filter(SPECIES_NAME == "OPILIO TANNER CRAB")
  
  #Filter haul data by summer months and NPT
  haul_opilio %>%
    mutate(sampdate = lubridate::mdy(HAUL_DATE), 
           year = lubridate::year(sampdate),
           month = lubridate::month(sampdate),
           day = lubridate::day(sampdate)) %>%
    filter(month %in% 6:8, GEAR_TYPE == "NON_PELAGIC_TRAWL") -> summer_npt
  
  #Define CRS, transform summer npt by crs, create spat vect=
  cc <- "EPSG:4326"
  
  summer_npt %>%
    dplyr::rename("x" = "GEAR_RETRIEVAL_LONDD", "y" = "GEAR_RETRIEVAL_LATDD") %>%
    akgfmaps::transform_data_frame_crs(out.crs = cc)-> test.points
  
  vect(test.points, geom = c("x", "y"), crs = cc) %>%
    project(cc) -> test.vect
  
  #Extract summer NPT bycatch data that falls within core area polygon
  extract(core_poly, test.vect) %>%
    cbind(test.points) %>%
    filter(is.na(id.x) == FALSE) -> summer_npt_core
  
  #Plot to make sure the extraction worked
  plot(core_poly)
  points(summer_npt_core$x, summer_npt_core$y)
  
  #Calculate total crab measured, total immature, and total mature males and the proportion of each by haul in core area
  length_opilio %>%
    filter(HAUL_JOIN %in% summer_npt_core$HAUL_JOIN) %>%
    group_by(HAUL_JOIN) %>%
    summarise(msr_crab = sum(FREQUENCY),
              imm_total = sum(SEX == "M" & LENGTH < 95),
              mat_total = sum(SEX == "M" & LENGTH >=95),
              imm_prop = imm_total/msr_crab,
              mat_prop = mat_total/msr_crab) -> length_sum
  
  #Calculate total number of crab extrapolated up to the haul by haul join, haul, and year
  summer_npt_core %>%
    group_by(HAUL_JOIN, VESSEL_ID, HAUL_SET_NUMBER, year) %>%
    summarise(catch_extrap = sum(NUMBER_OF_SPECIES_EXTRAPOLATED_UP_TO_THE_HAUL)) -> haul_sum
  
  #Join haul data with length data
  right_join(haul_sum, length_sum, by = "HAUL_JOIN") %>%
    summarise(year = year,
              catch_extrap = catch_extrap,
              imm_extrap = imm_prop * catch_extrap,
              mat_extrap = mat_prop * catch_extrap) %>%
    group_by(year) %>%
    summarise(catch_total = sum(catch_extrap),
              imm_total = sum(imm_extrap),
              mat_total = sum(mat_extrap)) -> bycatch_year_total
  
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