library(DBI)
library(ROracle)

drv <- dbDriver("Oracle")
connect.string <- 'AFSC'

## Input username and password (every user needs own) as objects. User name and password must be in parentheses

un<-"USERNAME"
pw<-"PASSWORD"

## Establish connection. NOTE: YOU MUST BE LOGGED INTO VPN TO CONNECT TO ORACLE

con <- dbConnect(drv, username =un, password=pw,dbname = connect.string)

## Submit query to Oracle database using standard SQL code language within parenthese
## This example targets haul table while excluding haul_type = 17 (retow hauls)
## This exclusion is appropriate in all cases except for female BBRKC
## A list of major tables to query from is included at the end of the code

rs <- dbSendQuery(con, "select * from CRAB.HAUL_NEWTIMESERIES where HAUL_TYPE <>17")

data <- fetch(rs)

View(data)

## Write CSV file with data

write.csv(data,"./Data/Haultable_Time_series.csv")


# List of critical tables to query from:
# 1.) HAUL_NEWTIMESERIES ; Haul table for EBS
# 2.) HAUL_NEWTIMESERIES_NBS ; Haul table for NBS
# 3.) EBSCRAB ; Specimen table for EBS...very, VERY, large
# 4.) EBSCRAB_NBS ; Specimen table for EBS
# 5.) STRATA_RKC_NEWTIMESERIES ; RKC strata table for EBS
# 6.) STRATA_BAIRDI_NEWTIMESERIES ; Bairdi strata table for EBS
# 7.) STRATA_OPILIO_NEWTIMESERIES ; Opilio strata table for EBS
# 8.) STRATA_BKC_NEWTIMESERIES ; BKC strata table for EBS
# 9.) STRATA_EI_NEWTIMESERIES ; BKC strata table for EBS
# 10.) STRATA_NBS_NEWTIMESERIES_STANDARDIZED_AREA ; Standardized strata table for NBS

