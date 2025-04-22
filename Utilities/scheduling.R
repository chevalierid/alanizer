rm(list=ls())
library(suncalc)

LATITUDE <- +49.250811969024724
LONGITUDE <- -123.2388863
MAX_FUTURE_DAYS <- 200

all_dates <- Sys.Date()-1 + 1:MAX_FUTURE_DAYS

times <- getSunlightTimes(date = all_dates,lat=LATITUDE,lon=LONGITUDE,tz="PST")
sunrises <- as.data.frame(as.numeric(as.POSIXct(times$sunrise)))
sunsets <- as.data.frame(as.numeric(as.POSIXct(times$sunset)))

alltimes <- c(t(cbind(sunrises,sunsets)))
alltimes

delta_t <- diff(alltimes)

delta_t

delta_t_str <- paste(toString(delta_t), sep = ",")

define <- sprintf("#define START_TIME %i;\nuint32_t delta_ts[] {%s};\nuint32_t sunrise_0 {%i}\n", sunrises[1,1], paste(delta_t_str, collapse=','), as.numeric(as.POSIXct(times$sunrise[1])))
cat(define,file = "schedule_data_v4.csv")
