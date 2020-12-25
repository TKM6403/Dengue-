rm(list=ls())
setwd("~/Desktop/stats")
library(spatialpred)
library(dengueThailand)
library(foreach)
library(ggplot2)
library(dplyr)
library(grid)
library(readr)
library(gridExtra)
library(cowplot)
library(reshape2)
library(rjags)
library(wesanderson)
library(cowplot)
library(scales)
library(pROC)
library(ResourceSelection)
library(readxl)
library(relaimpo)
library(ridge)
library(car)
library(mgcv)
library(glmnet)
library(caret)
library(tidyverse)
library(dynlm)
library(dyn)
library(readxl) 
library(astsa)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


#SOI data- Southern Oscillation Index
#Air Pressure Index
# <0 means El Nino
# >0 means La Nina
#averaging the index values over months or seasons helps 
#to isolate more sustained deviations from the average, 
#like those associated with ENSO.
SOI <-as.data.frame(read.csv("SOIdata.csv"))
SOI$SOI<- as.numeric(as.character(SOI$Southern.Oscillation.Index..SOI.))
str(SOI$SOI)
soibymonth <- c(SOI$SOI[734:838], NA, NA, NA)

singaporefull$SOI <- soibymonth

#ONI 2012-2020
ONI <- c(1.0,-0.8,	-0.6,	-0.5,	-0.4,	-0.2,	0.1,	0.3	,0.3,	0.3	,0.2,	0.0	,-0.2,
	-0.4,	-0.3,	-0.2,	-0.2,	-0.3,	-0.3,	-0.4,	-0.4,	-0.3,	-0.2,	-0.2,	-0.3,
	-0.4,	-0.4,	-0.2,	0.1,	0.3,	0.2,	0.1,	0.0,	0.2,	0.4,	0.6,	0.7,
	0.6,	0.6, 0.6,	0.8,	1.0,	1.2,	1.5,	1.8,	2.1,	2.4,	2.5,	2.6,
	2.5,	2.2,	1.7,	1.0,0.5,	0.0,	-0.3,	-0.6,	-0.7,	-0.7,	-0.7,	-0.6,
-0.3,	-0.1,	0.1,	0.3,	0.4,	0.4,	0.2,	-0.1,	-0.4,	-0.7,	-0.9,	-1.0,
	-0.9,	-0.8,	-0.6,	-0.4,	-0.1,	0.1,	0.1,	0.2,	0.4,	0.7,	0.9,	0.8,
	0.8,	0.8,	0.8,	0.7,	0.6,	0.5,	0.3,	0.1,	0.1,	0.3,	0.5,	0.5,0.5,	0.6,	0.5,	0.3	,0.0,	-0.2,	-0.4,	-0.6, NA, NA, NA)


#NDVI data
library(raster)
library(ncdf4)
library(rgdal)
library(ggplot2)
nc_data <- nc_open('LAI_mean_monthly_1981-2015.nc4')
# Save the print(nc) dump to a text file
{
  sink('monthlyLAI.txt')
  print(nc_data)
  sink()
}
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

head(lon)

ndvi.array <- ncvar_get(nc_data, "LAI")
dim(ndvi.array) 

fillvalue <- ncatt_get(nc_data, "LAI", "_FillValue")
fillvalue
nc_close(nc_data) 
ndvi.array[ndvi.array == fillvalue$value] <- NA
ndvi.slice <- ndvi.array[, , 1] 
dim(ndvi.slice)

#SAVE INTO A RASTER
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)

writeRaster(r, "LAI_1981.tif", "GTiff", overwrite=TRUE)

r_brick <- brick(ndvi.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r_brick <- flip(t(r_brick), direction='y')

toolik_lon <- -149.5975
toolik_lat <- 68.6275
toolik_series <- extract(r_brick, SpatialPoints(cbind(toolik_lon,toolik_lat)), method='simple')

toolik_df <- data.frame(year= seq(from=1981, to=2015, by=1), NDVI=t(toolik_series))

closeAllConnections()

fullsinglist <- read_excel_allsheets("Singapore2012-2020.xlsx")
singaporefull <- data.frame("Infected"=c(fullsinglist$`2012`$...12[2:53], fullsinglist$`2013`$...12[2:53], fullsinglist$`2014`$...12[2:53], fullsinglist$`2015`$...12[2:53], fullsinglist$`2016`$...12[2:53], fullsinglist$`2017`$...12[2:53], fullsinglist$`2018`$...12[2:53],fullsinglist$`2019`$...12[2:53],fullsinglist$`2020`$...12[2:53]), 
                            "Year"=c(seq(2012, 2012, length.out=52),seq(2013, 2013, length.out=52), seq(2014, 2014, length.out=52),seq(2015, 2015, length.out=52),seq(2016, 2016, length.out=52), seq(2017, 2017, length.out=52), seq(2018, 2018, length.out=52), seq(2019, 2019, length.out=52),seq(2020, 2020, length.out=52)),
                            "Cumweek"=c(seq(1,52*9, by=1)),
                            "Weekinyear"=c(seq(1, 52, length.out=52),seq(1, 52, length.out=52),seq(1, 52, length.out=52),seq(1, 52, length.out=52),seq(1, 52, length.out=52), seq(1, 52, length.out=52), seq(1, 52, length.out=52), seq(1, 52, length.out=52),seq(1, 52, length.out=52)),
                            "Population"=c(seq(5312000, 5850342, length.out=52*9)))
singaporefull$Population <- as.numeric(as.character(singaporefull$Population))
singaporefull$Infected <- as.numeric(as.character(singaporefull$Infected))
singaporefull$Year <- as.factor(singaporefull$Year)
singaporefull$incidence <- 1000000*singaporefull$Infected/singaporefull$Population

##Singapore Weekly Case Plot
p <- ggplot(data=singaporefull, aes(x=Cumweek, y=Infected, color=Year)) +
  geom_line()+
  geom_point()
p + xlab("Epi Week") + ylab("Cases") + ggtitle("Singapore Weekly Cases 2012-2020")+scale_color_grey()+labs(y= "Cases", x = "Epidemiological Week")




month <- c("January ", "February ", "March ", "April", "May", "June", "July", "August", "September", "October", "November", "December")
singaporemonthly <- data.frame("Month"=rep(month, times=9),
                               "Year"=c(seq(2012, 2012, length.out=12),seq(2013, 2013, length.out=12), seq(2014, 2014, length.out=12),seq(2015, 2015, length.out=12),seq(2016, 2016, length.out=12), seq(2017, 2017, length.out=12), seq(2018, 2018, length.out=12), seq(2019, 2019, length.out=12),seq(2020, 2020, length.out=12)),
                               "ContMonth"=c(seq(1, 12, length.out=12),seq(1, 12, length.out=12), seq(1, 12, length.out=12),seq(1, 12, length.out=12),seq(1, 12, length.out=12), seq(1, 12, length.out=12), seq(1, 12, length.out=12), seq(1, 12, length.out=12),seq(1, 12, length.out=12)),
                               "Population"=c(seq(5312000, 5850342, length.out=108)))


#WEEK TO MONTH
for (k in 0:8){
singaporemonthly$Infected[1+12*k] <- sum(singaporefull$Infected[(52*k+1):(52*k+5)],na.rm = TRUE)
singaporemonthly$Infected[2+12*k] <- sum(singaporefull$Infected[(52*k+6):(52*k+9)],na.rm = TRUE)
singaporemonthly$Infected[3+12*k] <- sum(singaporefull$Infected[(52*k+10):(52*k+13)],na.rm = TRUE)
singaporemonthly$Infected[4+12*k] <- sum(singaporefull$Infected[(52*k+14):(52*k+17)],na.rm = TRUE)
singaporemonthly$Infected[5+12*k] <- sum(singaporefull$Infected[(52*k+18):(52*k+21)],na.rm = TRUE)
singaporemonthly$Infected[6+12*k] <- sum(singaporefull$Infected[(52*k+22):(52*k+25)],na.rm = TRUE)
singaporemonthly$Infected[7+12*k] <- sum(singaporefull$Infected[(52*k+26):(52*k+30)],na.rm = TRUE)
singaporemonthly$Infected[8+12*k] <- sum(singaporefull$Infected[(52*k+31):(52*k+35)],na.rm = TRUE)
singaporemonthly$Infected[9+12*k] <- sum(singaporefull$Infected[(52*k+36):(52*k+39)],na.rm = TRUE)
singaporemonthly$Infected[10+12*k] <- sum(singaporefull$Infected[(52*k+40):(52*k+43)],na.rm = TRUE)
singaporemonthly$Infected[11+12*k] <- sum(singaporefull$Infected[(52*k+44):(52*k+47)],na.rm = TRUE)
singaporemonthly$Infected[12+12*k] <- sum(singaporefull$Infected[(52*k+48):(52*k+52)],na.rm = TRUE)

}
singaporemonthly$SOI <- soibymonth
singaporemonthly$Year <- as.factor(singaporemonthly$Year)
singaporemonthly$factormonth <- factor(singaporemonthly$Month, levels=c("January ", "February ", "March ", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
singaporemonthly$cummonth <- factor(c(seq(1:108)))
singaporemonthly$incidence <- 1000000*singaporemonthly$Infected/singaporemonthly$Population
singaporemonthly$ONI <- ONI
singaporemonthly <- singaporemonthly[1:105,]

#MONTHLY INCIDENCE BY YEAR
p<-ggplot(data=singaporemonthly, aes(x=factormonth, y=incidence, color=Year)) +
  geom_bar(stat="identity") +
  geom_line()+facet_wrap(~Year)
p

#SOI BY MONTH
q<-ggplot(data=singaporemonthly, aes(x=cummonth, y=SOI, color=Year)) +
  geom_line()+facet_wrap(~Year)
q

barplot(singaporemonthly$Infected, main="2013-2019",
     xlab="Cumulative Week", ylab="Infected",
     xlim=c(1, 108), ylim=c(-8000, 8000))
lines(singaporemonthly$SOI*4000)
#ENSO WITH LAG
library(dyn)

lag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  supsmu <- supsmu(x,y)
  fit.lm <- lm(y~x)
  summary1 <- summary(fit.lm)
  plot(x,y, main=sprintf('%s %s, %s %s', "Lag", time, "p-value=",  round(summary1$coefficients[2,4], 4) ), sub=sprintf('%s %s', "p-value", summary1$coefficients[2,4]), xlab = sprintf('%s', "variable"), ylab="incidence")
  lines(x,summary1$coefficients[1,1]+x*summary1$coefficients[2,1], col=2)
  return(summary1$coefficients[2,4]) #return p-value
}
# RUN, SEE PLOT, CHECK P-Values
sapply(timelags,lag,variable=singaporemonthly$SOI, incidence=singaporemonthly$incidence)

#CCF in R
lag2.plot <- function (series1, series2, max.lag = 0, corr = TRUE, smooth = TRUE, 
          col = gray(0.1), ...) 
{
  as.ts = stats::as.ts
  par = graphics::par
  plot = graphics::plot
  lines = graphics::lines
  ts.intersect = stats::ts.intersect
  legend = graphics::legend
  name1 = paste(deparse(substitute(series1)), "(t-", sep = "")
  name2 = paste(deparse(substitute(series2)), "(t)", sep = "")
  series1 = as.ts(series1)
  series2 = as.ts(series2)
  max.lag = as.integer(max.lag)
  m1 = max.lag + 1
  prow = ceiling(sqrt(m1))
  pcol = ceiling(m1/prow)
  a = stats::ccf(series1, series2, max.lag, plot = FALSE)$acf
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(prow, pcol), mar = c(2.5, 2.5, 1, 1), mgp = c(1.6, 
                                                              0.6, 0), cex.main = 1, font.main = 1)
  for (h in 0:max.lag) {
    plot(stats::lag(series1, -h), series2, xy.labels = FALSE, 
         main = paste(name1, h, ")", sep = ""), ylab = name2, 
         xlab = "", col = col, panel.first = Grid(), ...)
    if (smooth == TRUE) 
      lines(stats::lowess(ts.intersect(stats::lag(series1, 
                                                  -h), series2)[, 1], ts.intersect(stats::lag(series1, 
                                                                                              -h), series2)[, 2]), col = "red")
    if (corr == TRUE) 
      legend("topright", legend = round(a[m1 - h], digits = 3), 
             text.col = "blue", bg = "white", adj = 0.25, 
             cex = 0.85)
    on.exit(par(old.par))
  }
}


ONI_singapore <- ts(singaporemonthly$ONI)
SOI_singapore <- ts(singaporemonthly$SOI)
singapore_incidence <- ts(singaporemonthly$incidence)
x = ccf(ONI_singapore, singapore_incidence)
lag2.plot(ONI_singapore, incidence, 12)
lag2.plot(SOI_singapore, incidence, 12)

ONI_honduras <- ts(hondurasmonthly$ONI)
SOI_honduras <- ts(hondurasmonthly$SOI)
honduras_incidence <- ts(hondurasmonthly$incidence)
x = ccf(ONI_honduras, honduras_incidence)
lag2.plot(ONI_honduras, honduras_incidence, 12)
lag2.plot(SOI_honduras, honduras_incidence, 12)

sapply(timelags,lag,variable=singaporemonthly$ONI, incidence=singaporemonthly$incidence)




####HONDURAS FULL DATASET
x <- c(944.934,876.64,755.382,691.064,709.674,728.283,678.766,842.765,691.228,660.306,574.523,373.952,1009.276,1108.559,1207.842,1307.126,1388.796,1419.104,1561.047,1709.217,1611.954,1542.391,1646.387,1750.384,1498.452,1516.82,1292.716,1215.732,1068.435,909.389,871.542,766.308,661.074,555.839,450.605,375.381,324.869,274.356,312.885,343.612,343.612,388.512,479.52,496.357,513.195,471.577,421.065,492.545,576.732,478.484,576.538,674.592,662.785,624.158,610.592,615.184,619.776,572.055,502.302,568.949,499.494,430.04,232.015,413.86,574.247,568.635,563.022,589.952,627.171,687.824,742.375,691.799,717.056,626.779,690.762,721.987,525.549,519.774,498.867,453.406,407.944,395.028,386.008,376.988,356.556,297.962,239.367,180.773,249.85,221.437,193.024,276.044,239.65,177.347,162.725,148.103,133.481,118.859,109.844,150.254,142.84,122.996,103.151,100)
honduras2017 <- c(108.271,112.281,126.316,100.251,114.286,98.246,118.296,78.195,90.226,78.195,96.241,74.185,106.266,96.241,68.17,108.271,98.246,128.321,106.266,118.296,144.361,130.326,178.446,214.536,218.546,226.566,162.406,158.396,134.336,106.266,92.231,118.296,98.246,126.316,102.256,104.261,82.206,68.17,60.15,60.15,80.201,120.301,48.12,54.135,72.18,32.08,36.09,46.115,46.115,42.105, 38,40)
honduras2018 <- c(44.332,58.438,58.438,44.332,68.514,96.725,100.756,104.786,128.967,153.149,163.224,199.496,173.3,213.602,231.738,213.602,263.98,183.375,211.587,225.693,183.375,177.33,223.678,201.511,183.375,141.058,155.164,118.892,92.695,92.695,116.877,78.589,94.71,100.756,116.877,106.801,102.771,102.771,137.028,70.529,139.043,169.27,185.39,219.647,278.086,191.436,213.602,221.662,276.071,259.95,221.662,122.922)
honduras2019 <- read_excel("honduras.xlsx")
honduras2019inf <- subset(honduras2019$Infected, !is.na(honduras2019$Infected))
honduras2019infall <- c(honduras2019inf, seq(1000,20, by=-800/8) )
honduras2020 <- 0.009588*c(28429.752,47603.306,66776.86,80000,89256.198,109090.909,117024.793,104462.81,99834.711,101818.182,89256.198,59504.132,52892.562,50909.091,49586.777,47603.306,44958.678,38347.107,33719.008,23140.496,8595.041)
honduras2020 <- c(honduras2020, rep(NA, 52-21))

hondurasfull <- data.frame("Infected"=c(x,honduras2017,honduras2018,honduras2019infall,honduras2020), 
                            "Year"=c(seq(2015, 2015, length.out=52),seq(2016, 2016, length.out=52), seq(2017, 2017, length.out=52), seq(2018, 2018, length.out=52), seq(2019, 2019, length.out=52),seq(2020, 2020, length.out=52)),
                            "Cumweek"=c(seq(1,52*6, by=1)))
hondurasfull$Population <- seq(9113000, 9952671, length.out = 52*6)
hondurasfull$Year <- as.factor(hondurasfull$Year)
hondurasfull$incidence <- 1000000*hondurasfull$Infected/hondurasfull$Population
hondurasfull$Cum
p <- ggplot(data=hondurasfull, aes(x=Cumweek, y=Infected, color=Year)) +
  geom_line()+
  geom_point()
p + xlab("Epi Week") + ylab("Cases")# +facet_wrap(~Year)


month <- c("January ", "February ", "March ", "April", "May", "June", "July", "August", "September", "October", "November", "December")
hondurasmonthly <- data.frame("Month"=rep(month, times=6),
                               "Year"=c(seq(2015, 2015, length.out=12),seq(2016, 2016, length.out=12), seq(2017, 2017, length.out=12), seq(2018, 2018, length.out=12), seq(2019, 2019, length.out=12),seq(2020, 2020, length.out=12)),
                               "ContMonth"=c(seq(1, 12, length.out=12),seq(1, 12, length.out=12), seq(1, 12, length.out=12),seq(1, 12, length.out=12),seq(1, 12, length.out=12), seq(1, 12, length.out=12)),
                               "Population"=c(seq(9113000, 9952671, length.out=12*6)))

for (k in 0:5){
  hondurasmonthly$Infected[1+12*k] <- sum(hondurasfull$Infected[(52*k+1):(52*k+5)],na.rm = TRUE)
  hondurasmonthly$Infected[2+12*k] <- sum(hondurasfull$Infected[(52*k+6):(52*k+9)],na.rm = TRUE)
  hondurasmonthly$Infected[3+12*k] <- sum(hondurasfull$Infected[(52*k+10):(52*k+13)],na.rm = TRUE)
  hondurasmonthly$Infected[4+12*k] <- sum(hondurasfull$Infected[(52*k+14):(52*k+17)],na.rm = TRUE)
  hondurasmonthly$Infected[5+12*k] <- sum(hondurasfull$Infected[(52*k+18):(52*k+21)],na.rm = TRUE)
  hondurasmonthly$Infected[6+12*k] <- sum(hondurasfull$Infected[(52*k+22):(52*k+25)],na.rm = TRUE)
  hondurasmonthly$Infected[7+12*k] <- sum(hondurasfull$Infected[(52*k+26):(52*k+30)],na.rm = TRUE)
  hondurasmonthly$Infected[8+12*k] <- sum(hondurasfull$Infected[(52*k+31):(52*k+35)],na.rm = TRUE)
  hondurasmonthly$Infected[9+12*k] <- sum(hondurasfull$Infected[(52*k+36):(52*k+39)],na.rm = TRUE)
  hondurasmonthly$Infected[10+12*k] <- sum(hondurasfull$Infected[(52*k+40):(52*k+43)],na.rm = TRUE)
  hondurasmonthly$Infected[11+12*k] <- sum(hondurasfull$Infected[(52*k+44):(52*k+47)],na.rm = TRUE)
  hondurasmonthly$Infected[12+12*k] <- sum(hondurasfull$Infected[(52*k+48):(52*k+52)],na.rm = TRUE)
  
}
hondurasmonthly$SOI <- soibymonth[37:108]
hondurasmonthly$Year <- as.factor(hondurasmonthly$Year)
hondurasmonthly$factormonth <- factor(hondurasmonthly$Month, levels=c("January ", "February ", "March ", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
hondurasmonthly$cummonth <- factor(c(seq(1:72)))
hondurasmonthly$incidence <- 1000000*hondurasmonthly$Infected/hondurasmonthly$Population
hondurasmonthly$ONI <- ONI[37:108]
hondurasmonthly <- hondurasmonthly[1:65,]

#MONTHLY INCIDENCE BY YEAR
p<-ggplot(data=hondurasmonthly, aes(x=factormonth, y=incidence, color=Year)) +
  geom_bar(stat="identity") +
  geom_line()+facet_wrap(~Year)
p

#ENSO WITH LAG
library(dyn)

# RUN, SEE PLOT, CHECK P-Values
timelags <- c(1:12)
par(mfrow = c(3, 4))
par(mar=c(3,3,3,3))
#HONDURAS SOI VS INCIDENCE
sapply(timelags,lag,variable=hondurasmonthly$SOI, incidence=hondurasmonthly$incidence)
x <- (218.07-59.22*singaporemonthly$SOI[1:95])
y <- singaporemonthly$incidence[11:105]
ggplot <- data.frame("x"=x, "y"=y, "week"=c(1:95), "year"=singaporemonthly$Year[11:105])
p <- ggplot(ggplot, aes(x=week, y=y, col=year)) +geom_bar(stat="identity", alpha=0.75) 
p + geom_line(data=ggplot, aes(x=week, y=x), colour="red")+ggtitle("SOI vs Singapore Incidence")+xlab("Month")+ylab("incidence per 1000000")




sapply(timelags,lag,variable=hondurasmonthly$ONI, incidence=hondurasmonthly$incidence)



#SOI data- Southern Oscillation Index
#Air Pressure Index
# <0 means El Nino
# >0 means La Nina
#averaging the index values over months or seasons helps 
#to isolate more sustained deviations from the average, 
#like those associated with ENSO.


#ONI 


#PREPARE HONDURAS DATA



hweather <- NULL
hweather <- read_excel("HondurasWeather12-20.xlsx")
hweather$P <- hweather$`Precipitation integer in hundredths of a millimeter - liquid equivalent - "0" is used for trace amounts and "-1" is used for no precipitation`
hweather$P[hweather$P==0] <-0.1
hweather$P[hweather$P==-1] <-0 

hweather <- hweather[,2:13]
hweather <- as.matrix(hweather)

hondurasweather<-matrix(0, nrow=459, ncol=12)

for (k in 1:nrow(hondurasweather)){
for (i in 1:12){
  hondurasweather[k,i] <- (hweather[(7*k-6),i]+hweather[(7*k-5),i]+hweather[(7*k-4),i]+hweather[(7*k-3),i]+hweather[(7*k-2),i]+hweather[(7*k-1),i]+hweather[(7*k),i])/7
}
}
hondurasweather <- as.data.frame(hondurasweather)
library(tidyverse)
hondurasweather <- hondurasweather %>% 
  rename(
    avemaxtemp = V1,
    avemintemp = V2,
    avetemp = V3,
    heatingdays=V4,
    coolingdays=V5,
    DROP=V6,
    avehumidity=V7,
    avewind=V8,
    avedew=V9,
    avevisibility=V10,
    avesealevel=V11,
    aveprecipitation=V12
  )
honduras$aveprecipitation <- honduras$aveprecipitation/100
hondurasweather$cumrain <- hondurasweather$aveprecipitation*7

avehumid <- vector()


#PREPARE SINGAPORE DATA
sweather <- NULL
sweather <- read_excel("SingaporeWeather12-20.xlsx")
sweather$P <- sweather$`Precipitation integer in hundredths of a millimeter - liquid equivalent - "0" is used for trace amounts and "-1" is used for no precipitation`
sweather$P[sweather$P==0] <-0.1
sweather$P[sweather$P==-1] <-0 

sweather <- sweather[,2:13]
sweather <- as.matrix(sweather)

singaporeweather<-matrix(0, nrow=459, ncol=12)

for (k in 1:nrow(singaporeweather)){
  for (i in 1:12){
    singaporeweather[k,i] <- (sweather[(7*k-6),i]+sweather[(7*k-5),i]+sweather[(7*k-4),i]+sweather[(7*k-3),i]+sweather[(7*k-2),i]+sweather[(7*k-1),i]+sweather[(7*k),i])/7
  }
}
singaporeweather <- as.data.frame(singaporeweather)
library(tidyverse)
singaporeweather <- singaporeweather %>% 
  rename(
    avemaxtemp = V1,
    avemintemp = V2,
    avetemp = V3,
    heatingdays=V4,
    coolingdays=V5,
    DROP=V6,
    avehumidity=V7,
    avewind=V8,
    avedew=V9,
    avevisibility=V10,
    avesealevel=V11,
    aveprecipitation=V12
  )

singaporeweather$aveprecipitation <- singaporeweather$aveprecipitation/100
singaporeweather$cumrain <- singaporeweather$aveprecipitation*7

avehumid <- vector()


#CREATE CRUCIAL DATASETS
str(singaporefull)
str(singaporeweather)
singapore <- cbind(singaporefull[1:nrow(singaporeweather),], singaporeweather )
singapore$avevisibility[singapore$avevisibility<0] <- 9
singapore<- singapore[1:457,]
singapore$period <- seq(1:457)
singapore$period[1:60] <- "1"
singapore$period[61:193] <- "2"
singapore$period[194:416] <- "3"
singapore$period[417:457] <- "4"
decomposedRessing <- decompose(ts(singapore$incidence, frequency = 52), type="mult")
plot(decomposedRessing)
singapore$seasonal <- unclass(decomposedRessing$seasonal) 
#singapore$seasonal <- (singapore$seasonal/0.2958333)-2.318874 #scale proportionally to R0
ggplot(singapore, aes(x=Cumweek, y=seasonal, color=Year))+geom_point()+geom_line()+ggtitle("Singapore Seasonal Decomposition")+xlab("Cumulative Epi. Week")+ylab("Beta Value")

p <- ggplot(data=singapore, aes(x=Cumweek, y=Infected, color=Year)) +
  geom_line()+
  geom_point()
p + xlab("Cumulative epidemiological Week") + ylab("Cases") + ggtitle("Singapore Reported Cases")



str(hondurasfull)
str(hondurasweather)
hondurasweatheralt <- hondurasweather[157:459,]
hondurasfullalt <- hondurasfull[1:nrow(hondurasweatheralt),]
honduras <- cbind(hondurasfullalt, hondurasweatheralt)
honduras <- honduras[1:281,]
honduras$avevisibility[honduras$avevisibility<0] <- 9
honduras$period <- 1:281
honduras$period[1:208] <- "1"
honduras$period[209:281] <-"2"
p <- ggplot(data=honduras, aes(x=Cumweek, y=Infected, color=period)) +
  geom_line()+
  geom_point()
p + xlab("Cumulative epidemiological Week") + ylab("Cases") + ggtitle("Honduras Cases by Period")
decomposedReshond <- decompose(ts(honduras$incidence, frequency = 52), type="mult")
plot(decomposedReshond)
honduras$seasonal <- decomposedReshond$seasonal
honduras$seasonal <- unclass(honduras$seasonal)
#honduras$seasonal <- (honduras$seasonal*2.560601)-1.243248 #scale p
summary(honduras$seasonal)
ggplot(honduras, aes(x=Cumweek, y=seasonal, color=Year))+geom_point()+geom_line()+ggtitle("Honduras Seasonal Decomposition")+xlab("Cumulative Epi. Week")+ylab("Beta Value")




plot(honduras$seasonal)

for (t in 1:length(honduras$seasonal)){
  honduras$beta[t] <- (1.30768+0.32924*cos(2*3.14159262535*t/26.36067)+0.09073*sin(2*3.14159262535*t/26.31303))
}
lines(honduras$beta, lwd=3, col=2)
xc <- 1:length(honduras$beta)
xs <- 1:length(honduras$beta)
for (t in 1:length(honduras$beta)){
  xc[t] <- cos(2*pi*t/52)
  xs[t] <- sin(2*pi*t/52)
  fit <- lm(honduras$seasonal~xc+xs)
}

honduras$beta <- (1.31454-1.14983*xc+0.10990*xs)
plot(honduras$seasonal, main="Honduras Seasonal Decomposition", xlab="Cumulative Week", ylab="Beta Estimate")
lines(honduras$beta, lwd=3, col=2)

#honduras plots
colors <- c("Infected" = "magenta", "Temperature" = "Blue", "Precipitation" = " orange", "Cumulative Precipitation"= "green", "Humidity"="Red")


ggplot(honduras, aes(x=Cumweek)) + 
  geom_line(aes(y = Infected), color = "magenta1") + 
  geom_line(aes(y = avetemp*60), color="steelblue", linetype="solid") +
  geom_line(aes(y = aveprecipitation*60), color="darkorange1")+
  geom_line(aes(y = avedew*100), color="springgreen4" )+
  geom_line(aes(y = avehumidity*60), color= "red3")+
  labs(color="Legend")+
  scale_color_manual(name="Key",values = colors) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Cases",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./60, name="mm/%/C")
  )

ggplot(singapore, aes(x=cumweek)) + 
  geom_line(aes(y = Infected), color = "magenta1") + 
  geom_line(aes(y = avetemp*60), color="steelblue", linetype="solid") +
  geom_line(aes(y = aveprecipitation*60), color="darkorange1")+
  geom_line(aes(y = avedew*100), color="springgreen4" )+
  geom_line(aes(y = avehumidity*60), color= "red3")+
  labs(color="Legend")+
  scale_color_manual(name="Key",values = colors) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Cases",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./60, name="mm/%/C")
  )




#create seaprate plots by outbreak year
#p<-ggplot(allinfected, aes(x=Week, y=Infected, group=Year)) +
#  geom_line(aes(color=Year))+
#  geom_point(aes(color=Year))
#p + xlab("Epi Week") + ylab("Cases")

fitquant2020seasaonal <- read_csv("singapore2020quantiles-seasonal.csv")

credibleplot <- function(df, year){
  fitcases <- data.frame(mean=1:52, low=1:52, high=1:52, real=singapore$Infected[(52*year-104623):(52*year-104572)])
  for (t in 1:52){
    fitcases$mean[t] <- rnorm(1, df[t+52,3]/6, (df[t,3]))
    fitcases$low[t] <- rnorm(1, df[t+52,1]/6, (df[t,1]))
    fitcases$high[t] <- rnorm(1, df[t+52,5]/6, (df[t,5]))
    fitcases$week[t] <- t
  }
  #plot(fitcases$real, xlim=c(1,52),ylim=c(0,2000), main=sprintf('%s',year), xlab="Epi. Week", ylab="Reported Cases")
  
  #polygon(c(1:52, rev(1:52)),c(fitcases$low, rev(fitcases$high)), col='grey80')
  #lines(fitcases$real, lwd=3, col=2)
  #write.csv(df, sprintf('%s %s %s',"singapore", year, ".csv"))
  return(print(df[211,]/(singapore$Population[(52*year-104623)])))
}
#source('SimDengue.R')
singaporeJAGS <- function (year){
JAGS.sim.data<-NULL
JAGS.sim.data$C <- (singapore$Infected[(52*year-104623):(52*year-104572)])
#JAGS.sim.data$C <- singapore$Infected[365:416]
JAGS.sim.data$T<- length(JAGS.sim.data$C)
JAGS.sim.data$N<-singapore$Population[(52*year-104623):(52*year-104572)]
JAGS.sim.data$B<- 0.0001692308*JAGS.sim.data$N
JAGS.sim.data$rho1 <- 1/6
JAGS.sim.data$tauI <- 9.535602e-01
JAGS.sim.data$tauO <- 8.643903e-02
#JAGS.sim.data$tauR <- 4.981290e-01
#JAGS.sim.data$tauS<- 5.085361e-01
JAGS.sim.data$I_0 <- (singapore$Infected[(52*year-104624)])*6
JAGS.sim.data$seasonal <- singapore$seasonal[1:52]

#JAGS.sim.data$theta <- singapore$seasonal
#JAGS.sim.data$CASES<-as.numeric(CASES)
#JAGS.sim.data$J<-1
params <- c("Sp0", "I_0", "R_0", "Sp", "I", "R", "tauS", "tauI", "tauR","tauO","Sstep", "gamma", "remove", "theta", "CVar")

jmod <- jags.model("seasonaltheta.R", 
                   data=JAGS.sim.data,
                   n.chains=4,
                   n.adapt=10000)


#update(jmod, n.iter=7500)
 
jsamp <- coda.samples(jmod,
                      variable.names=params,
                      n.iter=50000,
                      thin=10)

fit <- summary(as.mcmc.list(jsamp))
fitquant <- as.data.frame(fit$quantiles)
return(credibleplot(fitquant,year))
}
par(mfrow = c(2, 4))
lapply(c(2012:2018), singaporeJAGS)
x <- 2
system("say Just finished darling. I am from Malaysia and think everyone is the enemy ! ")




#fitquantsingapore2020$est <- c( singapore$Infected[417:457], singapore$Infected[1:173])
singaporeSOest <- data.frame("2012"= c(0.5805082,0.6115396,0.6258038,0.6386093,0.6607451), "2013"=c(0.6226488,0.6384096,0.6458848,0.6530281,0.6662257), "2014"=c(0.6425418,0.6580222,0.6661782,0.674524,0.6885837), "2015"=c(0.6512513,0.6650233,0.6727587,0.6791138,0.6909084), "2016"=c(0.6802246,0.6941641,0.7005729,0.7076188,0.7193353), "2017"=c(0.6059835,0.6252681,0.6348157,0.6440285,0.6603158), "2018"=c(0.6225916,0.6392203,0.6476792,0.6556523,0.669757), "2019"=c(0.6383349,0.6529303,0.6599208,0.6667067,0.6790214), "2020"= c(0.6287421,0.650892,0.6602687,0.6689984,0.6839346))

plot(unlist(singaporeSOest[3,])~c(2012:2020),ylim=c(0.5,0.7), main="Initial Susceptible estimations", xlab="Year", ylab="% of populations")
polygon(c(2012:2020, rev(2012:2020)),c(singaporeSOest[1,], rev(singaporeSOest[5,])), col='grey80')
lines(unlist(singaporeSOest[3,])~c(2012:2020), lwd=3, col=2)



credibleplot <- function(df, year){
  fitcases <- data.frame(mean=1:52, low=1:52, high=1:52, real=honduras$Infected[(52*year-104779):(52*year-104728)])
  for (t in 1:52){
    fitcases$mean[t] <- rnorm(1, df[t+52,3]/6, (df[t,3]))
    fitcases$low[t] <- rnorm(1, df[t+52,1]/6, (df[t,1]))
    fitcases$high[t] <- rnorm(1, df[t+52,5]/6, (df[t,5]))
    fitcases$week[t] <- t
  }
  dev.off()
  plot(fitcases$real, xlim=c(1,52),ylim=c(0,2000), main=sprintf('%s',year), xlab="Epi. Week", ylab="Reported Cases")
  
  polygon(c(1:52, rev(1:52)),c(fitcases$low, rev(fitcases$high)), col='grey80')
  lines(fitcases$real, lwd=3, col=2)
  #write.csv(df, sprintf('%s %s %s',"honduras", year, ".csv"))
  print(df[211,]/(honduras$Population[(52*year-104779)]))
}

hondurasJAGS <- function (year){
  JAGS.sim.data<-NULL
  JAGS.sim.data$C <- (honduras$Infected[(52*2019-104779):(52*2019-104728)])
  #JAGS.sim.data$C <- singapore$Infected[365:416]
  JAGS.sim.data$T<- length(JAGS.sim.data$C)
  JAGS.sim.data$N<-honduras$Population[(52*2019-104779):(52*2019-104728)]
  JAGS.sim.data$B<- 0.0004230769*JAGS.sim.data$N
  JAGS.sim.data$rho1 <- 1/6
  JAGS.sim.data$theta <-1.3
  JAGS.sim.data$I_0 <- (honduras$Infected[(52*2019-104780)])*6
  
  JAGS.sim.data$tauI <- 9.890516e-01
  JAGS.sim.data$tauO <- 7.737229e-02
  JAGS.sim.data$tauR <- 4.981290e-01
  JAGS.sim.data$tauS<- 5.085361e-01
  
  #JAGS.sim.data$theta <- 1.3
  #JAGS.sim.data$theta <- singapore$seasonal
  #JAGS.sim.data$CASES<-as.numeric(CASES)
  #JAGS.sim.data$J<-1
  params <- c("Sp0", "I_0", "R_0", "Sp", "I", "R", "tauS", "tauI", "tauR","tauO","Sstep", "gamma", "remove", "theta", "CVar")
  
  jmod <- jags.model("constanttheta.jags.R", 
                     data=JAGS.sim.data,
                     n.chains=4,
                     n.adapt=10000)
  
  
  #update(jmod, n.iter=7500)
  
  jsamp <- coda.samples(jmod,
                        variable.names=params,
                        n.iter=50000,
                        thin=10)
  
  fit <- summary(as.mcmc.list(jsamp))
  fitquant <- as.data.frame(fit$quantiles)
  return(credibleplot(fitquant,year))
}
lapply(c(2015:2018), hondurasJAGS)

singaporeyears <- c(2013:2019)
lapply(singaporeyears, singaporeJAGS)
singaporetheta <-(fitquantsingapore[1836:2291,3])
singaporelambda<-fitquantsingapore[1374:1830,3]
#For now singapore values cover 2017-2019
#2012-2016 is 52*5=260 values
plot(fitquanthonduras[1:281,3])




par(mfrow = c(4, 4))
par(mar=c(3,3,3,3))
timelags <- c(1:16)


lag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  supsmu <- supsmu(x,y)
  fit.lm <- lm(y~x)
  summary1 <- summary(fit.lm)
  plot(x,y, main=sprintf('%s %s, %s %s', "Lag", time, "p-value=",  round(summary1$coefficients[2,4], 4) ), sub=sprintf('%s %s', "p-value", summary1$coefficients[2,4]), xlab = sprintf('%s', "variable"), ylab="incidence")
  lines(x,summary1$coefficients[1,1]+x*summary1$coefficients[2,1], col=2)
  return(summary1$coefficients[2,4]) #return p-value
}

quadlag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  supsmu <- supsmu(x,y)
  x2 <- x^2
  fit.lm <- lm(y~x2+x)
  summary1 <- summary(fit.lm)
  plot(x,y, main=sprintf('%s %s', "Lag", time ), sub=sprintf('%s %s', "R-squared=", summary1$r.squared), xlab = sprintf('%s', "variable"), ylab="incidence")
  lines(x,summary1$coefficients[1,1]+x2*summary1$coefficients[2,1]+x*summary1$coefficients[3,1], col=2)
  return(pf(summary1$fstatistic[1],summary1$fstatistic[2], summary1$fstatistic[3], lower.tail = F))
}

expolag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  supsmu <- supsmu(x,y)
  fit.lm <- lm(log(y)~x)
  summary1 <- summary(fit.lm)
  plot(x,y, main=sprintf('%s %s, %s %s', "Lag", time, "p-value=",  round(summary1$coefficients[2,4], 4) ), sub=sprintf('%s %s', "p-value", summary1$coefficients[2,4]), xlab = sprintf('%s', "variable"), ylab="incidence")
  lines(x,exp(summary1$coefficients[1,1]+x*summary1$coefficients[2,1]), col=2)
  return(summary1$coefficients[2,4])
}

library(mgcv)
####Generalized additive models
additivelag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  df <- data.frame(x=x, y=y)
  fit.lm <- gam(y~s(x))
  summary1 <- summary(fit.lm)
  plot <-  ggplot(df, aes(x,y))+
    geom_point()+
    stat_smooth(method=gam, formula=y~s(x))
return(summary1$s.pv)
}


timelags <- c(1:16)

sapply(timelags,lag,variable=singapore$avemaxtemp[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avemaxtemp[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avemaxtemp[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avehumidity[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avehumidity[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avehumidity[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singaporeweather$avevisibility[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singaporeweather$avevisibility[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singaporeweather$avevisibility[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$cumrain[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avevisibility[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avevisibility[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avemintemp[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avemintemp[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avemintemp[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avewind[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avewind[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avewind[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avesealevel[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avesealevel[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avesealevel[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avetemp[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avetemp[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avetemp[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$avedew[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$avedew[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$avedew[1:457], incidence=singaporelambda)

sapply(timelags,lag,variable=singapore$aveprecipitation[1:457], incidence=singaporelambda)
sapply(timelags,quadlag,variable=singapore$aveprecipitation[1:457], incidence=singaporelambda)
sapply(timelags,expolag,variable=singapore$aveprecipitation[1:457], incidence=singaporelambda)


sapply(timelags,lag,variable=singapore$avemaxtemp, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$avehumidity, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$avevisibility, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$cumrain, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$avemintemp, incidence=singapore$incidence) #may be quadratic?
sapply(timelags,lag,variable=singapore$avewind, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$avesealevel, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$avetemp, incidence=singapore$incidence) #may also be quadratic
sapply(timelags,lag,variable=singapore$avedew, incidence=singapore$incidence)
sapply(timelags,lag,variable=singapore$aveprecipitation, incidence=singapore$incidence)



###estimate function for theta
s <- supsmu(singapore$Cumweek, singapore$Infected, bass=7)
plot(singapore$Cumweek, singapore$Infected)
lines(s, lwd=2, col=2)
closeAllConnections()
####HONDURAS JAGS
hondurasjags <- function(year){
JAGS.sim.data<-NULL
JAGS.sim.data$C <- honduras$Infected
JAGS.sim.data$T<- length(JAGS.sim.data$C)
JAGS.sim.data$N<-seq(honduras$Population[1], honduras$Population[length(honduras$Population)], length.out=length(JAGS.sim.data$C))
JAGS.sim.data$B<- 0.0216*JAGS.sim.data$N/52
JAGS.sim.data$rho1 <- 0.20
#JAGS.sim.data$CASES<-as.numeric(CASES)
#JAGS.sim.data$J<-1
params <- c("Sp", "tauS", "tauO", "tauP","theta", "gamma", "lambda", "Sp0", "I" ,"theta", "C")

jmod <- jags.model("simdengue.jags", 
                   data=JAGS.sim.data,
                   n.chains=4,
                   inits = list(tauS=0.1),
                   n.adapt=100000)


#update(jmod, n.iter=7500)

jsamp <- coda.samples(jmod,
                      variable.names=params,
                      n.iter=100000,
                      thin=10)

fit <- summary(as.mcmc.list(jsamp))
fitquant<- as.data.frame(fit$quantiles)
write_xlsx(fitquanthonduras, sprintf('%s %s %s', "HondurasYear", "2020", ".xlsx"))

}
lapply(c(2015:2019), hondurasjags)
View(fitquanthonduras)
honduraslambda <- fitquanthonduras[846:1126,3]
#1-104 is 2015-2016 in honduras
#theta spans from 1133-1411 
#2017-2019 should be 1237-1392


hondurastheta <-(fitquanthonduras[1133:1411,3])
hondurasthetafornow <- (fitquanthonduras[1237:1392,3])
hondurasinfectedfornow <- fitquanthonduras[105:260,3]
par(mfrow = c(1, 1))
par(mar=c(3,3,3,3))
timelags <- c(1:16)

plot(hondurastheta, hondurasfull$incidence[1:279], main="Theta v Incidence- Honduras", xlab = "theta from JAGS", ylab="Honduras Incidence")
#TO TEST:[1] "avemaxtemp"       "avemintemp"       "avetemp"         
#[4] "heatingdays"      "coolingdays"      "DROP"            
#[7] "avehumidity"      "avewind"          "avedew"          
#[10] "avevisibility"    "avesealevel"      "aveprecipitation"
#[13] "cumrain"          "Infected"         "cumweek"         
#[16] "year"    

#READ in FITQUANTS
singaporefitquant <- read_excel("Singaporefitquant.xlsx")
x <- c(which(singaporefitquant$parameter=="theta"))
singaporetheta <- vector()
for (i in 1:length(x)){
singaporetheta[i]<-singaporefitquant$`50%`[x[i]]
}
plot(singaporetheta)

x <- c(which(singaporefitquant$parameter=="S0"))
singaporeinitsusc <- vector()
for (i in 1:length(x)){
  singaporeinitsusc[i]<-singaporefitquant$`50%`[x[i]]
}

x <- c(which(singaporefitquant$parameter=="lambda"))
singaporelambda <- vector()
for (i in 1:length(x)){
  singaporelambda[i]<-singaporefitquant$`50%`[x[i]]
}

hondurasfitquant <- read_excel("HondurasFitquant.xlsx")
x <- c(which(hondurasfitquant$parameter=="lambda"))
honduraslambda <- vector()
for (i in 1:length(x)){
  honduraslambda[i]<-hondurasfitquant$`50%`[x[i]]
}
plot(honduraslambda)

x <- c(which(hondurasfitquant$parameter=="theta"))
hondurastheta <- vector()
for (i in 1:length(x)){
  hondurastheta[i]<-hondurasfitquant$`50%`[x[i]]
}
plot(hondurastheta)
#REGRESS CLIMATE VARIABLES VS THETA APPROXIMATION
sapply(timelags,lag,variable=hondurasweather$avemaxtemp, incidence=hondurasthetafornow)
sapply(timelags,expolag,variable=honduras$avehumidity, incidence=honduraslambda)
sapply(timelags,lag,variable=hondurasweather$avevisibility, incidence=honduraslambda)
sapply(timelags,lag,variable=honduras$avevisibility, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$cumrain, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$avemintemp, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$avewind, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$avesealevel, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$avetemp, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$avedew, incidence=hondurasthetafornow)
sapply(timelags,lag,variable=honduras$aveprecipitation, incidence=hondurasthetafornow)


sapply(timelags,lag,variable=honduras$avemaxtemp, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avehumidity, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avevisibility, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$cumrain, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avemintemp, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avewind, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avesealevel, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avetemp, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$avedew, incidence=honduras$incidence)
sapply(timelags,lag,variable=honduras$aveprecipitation, incidence=honduras$incidence)

#get best climate fits for each varible
regressiontype <- function(variable){
 x<-  sapply(timelags,lag,variable=variable[1:457], incidence=singaporelambda)
 bestx <- which(x==min(x))
 y <- sapply(timelags,quadlag,variable=variable[1:457], incidence=singaporelambda)
 besty <- which(y==min(y))
 z <- sapply(timelags,expolag,variable=variable[1:457], incidence=singaporelambda)
bestz<-which(z==min(z))

  return(c(bestx,min(x),besty,min(y),bestz,min(z)))
 }

vectorlist <- list(singapore$avemaxtemp, singapore$avehumidity, singapore$avevisibility,singapore$cumrain,singapore$avemintemp, singapore$avewind,singapore$avesealevel,singapore$avetemp,singapore$avedew,singapore$aveprecipitation)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
library(writexl)
write_xlsx(x, "singaporeclimate~lambdanew.xlsx")


regressiontype <- function(variable){
  x<-  sapply(timelags,lag,variable=variable[1:457], incidence=singapore$incidence[1:457])
  bestx <- which(x==min(x))
  y <- sapply(timelags,quadlag,variable=variable[1:457], incidence=singapore$incidence[1:457])
  besty <- which(y==min(y))
  z <- sapply(timelags,expolag,variable=variable[1:457], incidence=singapore$incidence[1:457])
  bestz<-which(z==min(z))
  w<- sapply(timelags,additivelag,variable=variable[1:457], incidence=singapore$incidence[1:457])
  bestw<-which(w==min(w))
  
  return(c(bestx,scientific(min(x)),besty,scientific(min(y)),bestz,scientific(min(z)), bestw, scientific((min(w)))))
}

vectorlist <- list(singapore$avemaxtemp, singapore$avehumidity, singapore$avevisibility,singapore$cumrain,singapore$avemintemp, singapore$avewind,singapore$avesealevel,singapore$avetemp,singapore$avedew,singapore$aveprecipitation)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
write_xlsx(x, "singaporeclimate~incidenceJOURNAL.xlsx")

regressiontype <- function(variable){
  x<-  sapply(timelags,lag,variable=variable, incidence=singaporemonthly$incidence)
  bestx <- which(x==min(x))
  y <- sapply(timelags,quadlag,variable=variable, incidence=singaporemonthly$incidence)
  besty <- which(y==min(y))
  z <- sapply(timelags,expolag,variable=variable, incidence=singaporemonthly$incidence)
  bestz<-which(z==min(z))
  w<- sapply(timelags,additivelag,variable=variable, incidence=singaporemonthly$incidence)
  bestw<-which(w==min(w))
  
  return(c(bestx,scientific(min(x)),besty,scientific(min(y)),bestz,scientific(min(z)), bestw, scientific((min(w)))))
}

vectorlist <- list(singaporemonthly$ONI, singaporemonthly$SOI)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
write_xlsx(x, "singaporeclimate~incidenceJOURNAL.xlsx")



### DO THE SAME FOR HONDURAS
regressiontype <- function(variable){
  x<-  sapply(timelags,lag,variable=variable, incidence=honduraslambda)
  bestx <- which(x==min(x))
  y <- sapply(timelags,quadlag,variable=variable, incidence=honduraslambda)
  besty <- which(y==min(y))
  z <- sapply(timelags,expolag,variable=variable, incidence=honduraslambda)
  bestz<-which(z==min(z))
  w <- sapply(timelags,additivelag,variable=variable, incidence=honduraslambda)
  bestw<-which(w==min(w))
  
  return(c(bestx,min(x),besty,min(y),bestz,min(z), bestw, min(w)))
}

vectorlist <- list(honduras$avemaxtemp, honduras$avehumidity, honduras$avevisibility,honduras$cumrain,honduras$avemintemp, honduras$avewind,honduras$avesealevel,honduras$avetemp,honduras$avedew,honduras$aveprecipitation)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
library(writexl)
write_xlsx(x, "hondurasclimate~lambda.xlsx")


regressiontype <- function(variable){
  x<-  sapply(timelags,lag,variable=variable, incidence=honduras$incidence)
  bestx <- which(x==min(x))
  y <- sapply(timelags,quadlag,variable=variable, incidence=honduras$incidence)
  besty <- which(y==min(y))
  z <- sapply(timelags,expolag,variable=variable, incidence=honduras$incidence)
  bestz<-which(z==min(z))
  w <- sapply(timelags,additivelag,variable=variable, incidence=honduras$incidence)
  bestw<-which(w==min(w))
  
  
  return(c(bestx,scientific(min(x)),besty,scientific(min(y)),bestz,scientific(min(z)), bestw, scientific((min(w)))))
}

vectorlist <- list(honduras$avemaxtemp, honduras$avehumidity, honduras$avevisibility,honduras$cumrain,honduras$avemintemp, honduras$avewind,honduras$avesealevel,honduras$avetemp,honduras$avedew,honduras$aveprecipitation)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
write_xlsx(x, "hondurasclimate~incidenceforpaper.xlsx")




regressiontype <- function(variable){
  x<-  sapply(timelags,lag,variable=variable, incidence=hondurasmonthly$incidence)
  bestx <- which(x==min(x))
  y <- sapply(timelags,quadlag,variable=variable, incidence=hondurasmonthly$incidence)
  besty <- which(y==min(y))
  z <- sapply(timelags,expolag,variable=variable, incidence=hondurasmonthly$incidence)
  bestz<-which(z==min(z))
  w<- sapply(timelags,additivelag,variable=variable, incidence=hondurasmonthly$incidence)
  bestw<-which(w==min(w))
  
  return(c(bestx,scientific(min(x)),besty,scientific(min(y)),bestz,scientific(min(z)), bestw, scientific((min(w)))))
}

vectorlist <- list(hondurasmonthly$ONI, hondurasmonthly$SOI)
x <- lapply(vectorlist, regressiontype)
x <- as.data.frame(x)
write_xlsx(x, "hondurasmonthlyclimate~incidenceJOURNAL.xlsx")
########
#########




bestlag <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  supsmu <- supsmu(x,y)
  fit.lm <- lm(y~x)
  summary1 <- summary(fit.lm)
  return(summary1$coefficients[2,4])
  }
whichlag <- function (country, variable){
  x <- sapply(timelags,bestlag,variable=variable, incidence=country$incidence) #get all p-values
  y <- which(x==min(x)) #capture min p-value
  return(sprintf('%s %s %s, %s %s', "Lag", y, "week(s)",  "p=", formatC(x[y],format="e"))) 
  
}

library(htmlTable)
#tbl <- matrix(c(whichlag(singapore, singapore$avemaxtemp), whichlag(singapore, singapore$avehumidity), whichlag(singapore, singapore$avevisibility), whichlag(singapore, singapore$cumrain), whichlag(singapore, singapore$avemintemp),whichlag(singapore, singapore$avewind), whichlag(singapore, singapore$avesealevel), whichlag(singapore, singapore$avetemp), whichlag(singapore, singapore$avedew), whichlag(singapore, singapore$aveprecipitation),whichlag(singapore, singapore$avemaxtemp), whichlag(singapore, singapore$avehumidity), whichlag(singapore, singapore$avevisibility), whichlag(singapore, singapore$cumrain), whichlag(singapore, honduras$avemintemp),whichlag(honduras, honduras$avewind), whichlag(honduras, honduras$avesealevel), whichlag(honduras, honduras$avetemp), whichlag(honduras, honduras$avedew), whichlag(honduras, honduras$aveprecipitation) ),
#              ncol = 2,
#              dimnames = list(rows=c("Max Temp", "Humidity", "Visibility", "Cumulative Rain", "Minimum Temperature", "Wind", "Sea Level Pressure", "Average Temperature", "Dew Point ", "Average Precipitation"),
#                              cols= c("Singapore", "Honduras"))) %>% 
#  htmlTable
#(tbl)



#APPLICATION OF CLIMATE TO INCIDENCE
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


climateprediction <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  df <- data.frame("Variable" =x,
                   "Incidence"=y)
  fit.lm <- lm(y~x)
  s <- summary(fit.lm)
  forplot <- s$coefficients[1,1]+s$coefficients[2,1]*x
  plot(y, main=sprintf('%s %s', "R^2=", s$r.squared))
  return(lines(forplot, col=2, lwd=2))
}
quadclimateprediction <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  df <- data.frame("Variable" =x,
                   "Incidence"=y)
  x2 <- x^2
  fit.lm <- lm(y~x2+x)
  s <- summary(fit.lm)
  forplot <- s$coefficients[1,1]+s$coefficients[2,1]*x2+s$coefficients[3,1]*x
  plot(y, main=sprintf('%s %s', "R^2=", s$r.squared))
  return(lines(forplot, col=2, lwd=2))
}
expoclimateprediction <- function(variable,incidence, time){
  x <- c(variable[1: (length(variable)-time)])
  y <- c(incidence[(1+time):length(incidence)]) 
  df <- data.frame("Variable" =x,
                   "Incidence"=y)
  x2 <- x^2
  fit.lm <- lm(log(y)~x)
  s <- summary(fit.lm)
  forplot <- s$coefficients[1,1]+s$coefficients[2,1]*x2+s$coefficients[3,1]*x
  plot(y, main=sprintf('%s %s', "R^2=", s$r.squared))
  return(lines(forplot, col=2, lwd=2))
}
#gelcman.diag(jsamp)
par(mfrow = c(1, 1))
par(mar=c(3,3,3,3))
quadclimateprediction(singapore$avedew, singapore$incidence, 7)
quadclimateprediction(singapore$avemaxtemp, singapore$incidence, 14)
a<- (singapore$avemaxtemp[3:443])^2
b<- (singapore$avedew[10:450])^2
c <- (singapore$avetemp[3:445])^2
d <- exp(singapore$avesealevel[1:441])

fit.lm <- lm(singapore$incidence[17:457]~b +singapore$avedew[10:450])
x <- summary(fit.lm)

fit.lm <- lm(log(honduras$incidence[12:281])~(honduras$avehumidity[1:270]))
x <- summary(fit.lm)
plot(honduras$incidence, main=sprintf('%s %s', "R^2=", x$r.squared))
lines(exp(7.2253722-0.049934*honduras$avehumidity[1:270]), col=2, lwd=2)

fit.lm <- lm(log(honduras$incidence[13:281])~(honduras$avehumidity[2:270])+honduras$avedew[1:269])
x <- summary(fit.lm)


fit.lm <- lm(singapore$incidence[17:457]~a +singapore$avemaxtemp[3:443] + b +singapore$avedew[10:450]+singapore$avewind[8:448])
x <- summary(fit.lm)
plot(singapore$incidence[17:457], main=sprintf('%s %s', "R^2=", x$r.squared))
lines(9223.728+4.870*a-302.297*singapore$avemaxtemp[3:443]+8.644*b-394.599*singapore$avedew[10:450]-1.882*singapore$avewind[8:448], col=2, lwd=3)
fit.lm <- lm(honduras$incidence[17:281]~honduras$avemaxtemp[16:280]+honduras$avevisibility[4:268]+honduras$cumrain[1:265])

plot(honduras$incidence)
lines(3.64*honduras$avemaxtemp[16:280]-41.78*honduras$avevisibility[4:268]+0.006774*honduras$cumrain[1:265], lwd=3, col=2)

#coda::traceplot(jsamp)
#summary(jsamp)

#a <- as.mcmc.list(jsamp)
#plot(a)
#summary(a)
install.packages("ncdf4_1.15.tar.gz", repos = NULL, type="source")

library(raster)
library(sp)
library(MODIS)
library(gdalUtils)
gdalinfo("MOD17A3H.A2000001.h21v09.006.2015141183401.hdf")


  x <- c(singapore$avedew[1: (length(singapore$avedew)-7)])
  y <- c(singapore$incidence[(1+7):length(singapore$incidence)]) 
  supsmu <- supsmu(x,y)
  mean <- mean(x)
  fit.lm <- lm(y~ (x-mean)^2 )
  summary1 <- summary(fit.lm)
  plot(x,y, main=sprintf('%s %s, %s %s', "Lag", time, "p-value=",  round(summary1$coefficients[2,4], 3) ), sub=sprintf('%s %s', "p-value", summary1$coefficients[2,4]), xlab = sprintf('%s', "variable"), ylab="incidence")
  lines(x,summary1$coefficients[1,1]+x*summary1$coefficients[2,1], col=2)
  return(summary(fit.lm))

