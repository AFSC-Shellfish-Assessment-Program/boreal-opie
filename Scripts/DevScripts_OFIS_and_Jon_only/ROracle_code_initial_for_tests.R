R.home(component = "home")
Sys.getenv("PATH")
Sys.getenv("OCI_INC")
Sys.getenv("OCI_LIB64")
Sys.getenv("ORACLE_HOME")
Sys.getenv("TNS_ADMIN")
install.packages("DBI")
library(DBI)

library(ROracle)

drv <- dbDriver("Oracle")

connect.string <- 'AFSC'

con <- dbConnect(drv, username ="<your username>", password="<your password>",dbname = connect.string)

rs <- dbSendQuery(con, "select * from global_name")

data <- fetch(rs)

View(data)
