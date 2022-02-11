install.packages('RODBC')

# Load RODBC package
library(RODBC)

# Create a connection to the CRAB database called "dc"
pw <- readline("Input Password: ")

dc <- odbcConnect("AFSC", uid="crab",  pwd=pw)

# Query the desired datatable and place the query output into the data frame "df". 
# See list at end of code for major tables of interest.

 df <- sqlQuery(dc, 
 "SELECT *
 FROM
 CRAB.HAUL_NEWTIMESERIES where HAUL_TYPE<>17")

# When finished, it's a good idea to close the connection
odbcClose(channel)


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
