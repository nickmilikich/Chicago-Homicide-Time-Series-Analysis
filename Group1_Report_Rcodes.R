############################
# Data import and processing
############################

setwd("~/Google\ Drive\ (nmilikic@nd.edu)/Fall\ 2019/ACMS\ 40842\ Time\ Series\ Analysis/Time\ Series\ Project/")
hom = read.csv(file = "Homicides (1).csv", header = TRUE, sep = ",")
hom.dummy = numeric(0)
currMonth = 1
currYear = 2001
currVal = 0
for (i in seq(length(hom$Date),1,-1))
{
  month = as.numeric(substr(toString(hom[i,"Date"]),0,2))
  year = as.numeric(substr(toString(hom[i,"Date"]),7,10))
  if (month == currMonth & year == currYear)
    currVal = currVal + 1
  else
  {
    hom.dummy = c(hom.dummy, currVal)
    currVal = 0
    if (currMonth == 12)
    {
      currMonth = 1
      currYear = currYear + 1
    }
    else
      currMonth = currMonth + 1
  }
}

hom.ts = ts(hom.dummy, frequency = 12, start = c(2001,1))

  
###############
# Visualization
###############

library(fpp2)
autoplot(hom.ts) + xlab("Month/Year") + ylab("Monthly Homicides") + ggtitle("Monthly Chicago Homicides")
ggseasonplot(hom.ts, year.labels = TRUE) + ylab("Monthly Homicides") + ggtitle("Seasonal Plot: Monthly Chicago Homicides")
ggsubseriesplot(hom.ts) + ylab("Monthly Homicides") + ggtitle("Subseries Plot: Monthly Chicago Homicides")
ggAcf(hom.ts) + ggtitle("ACF: Monthly Chicago Homicides")
pacf(hom.ts)

hom.train = window(hom.ts, end = c(2017,09))
hom.test = window(hom.ts, start = c(2017,10))

##############################
# Time Series Regression Model
##############################

hom.tslm <- tslm(hom.train~trend+season)
summary(hom.tslm)
# The season coefficients show a peak at season 7 (July), suggesting that season impacts the number of homicides
#   The seasons where temperature is normally higher (May through September/October) are positive
#   meaning that the number of homicides is expected to increase during those months.
#Therefore, our initial attempt to analyze the linear model with external Temperature data was unnecessary

autoplot(forecast(hom.tslm, h = 24), PI = FALSE) + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
checkresiduals(hom.tslm)
# Acf plot shows strong trend

#Breusch Godfrey Test
# p-value: 2.062e-12 -> there is some correlation of the lags

accuracy(forecast(hom.tslm, h = 24), hom.test)
#RMSE (test): 8.001407
library(fpp)
CV(hom.tslm) #1016.756

#############
# TSLM ARIMA ERRORS
#############
temp = read.csv(file = "Temperature.txt", header = TRUE, skip = 4, sep = ",")
temp$Value[143] <- 0.5*(temp$Value[142]+temp$Value[144])
temp$Anomaly[143] <- 0.5*(temp$Anomaly[142]+temp$Anomaly[144])
temp.ts = ts(temp$Value, frequency = 12, start = c(2001,1))
# temp.ts[143] is missing data that has been replaced as the average between temp.ts[142] and temp.ts[144]
temp.train = window(temp.ts, end = c(2017,09))
temp.test = window(temp.ts, start = c(2017,10))
temp = temp[-226,]

hom.dumm2 <- cbind(hom.dummy, temp)
hom.ts2 = ts(hom.dumm2, frequency = 12, start = c(2001,1))
hom.train2 <- window(hom.ts2, end = c(2017,09))
hom.test2 = window(hom.ts2, start = c(2017,10))

hom.arimaerrors <- auto.arima(hom.train2[,"hom.dummy"], xreg = hom.train2[,"Value"])
summary(hom.arimaerrors)
#ARIMA(0,1,1)(2,0,0)[12]
#AICc = 1445.78
#RMSE = 8.705785
autoplot(forecast(hom.arimaerrors, xreg = hom.test2[,"Value"], h =24),PI = FALSE) + autolayer(hom.test2[,"Value"]) + xlab("Month/Year")
ggtsdisplay(residuals(hom.arimaerrors, type="response"), main = "Regression Errors")
checkresiduals(hom.arimaerrors)
CV(hom.arimaerrors)

#############
# ARIMA Model
#############

hom.arima = auto.arima(hom.train, lambda = BoxCox.lambda(hom.ts))
summary(hom.arima)
#ARIMA(2,1,1)(2,1,2)[12]
# This model has a first order differencing with a seasonal lag of m = 12
autoplot(forecast(hom.arima, h = 24), PI = FALSE) + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
checkresiduals(hom.arima)
# Significant lag at ~35
#LJung-Box test
# p-value = 0.8883 -> fail to reject that the time series is not autocorrelated

accuracy(forecast(hom.arima, h = 24), hom.test)
#RMSE (test): 14.793524

###########
# Differencing Graphs
###########

hom.ts.differenced <- diff(hom.ts, lag = 12)
autoplot(hom.ts.differenced) + ggtitle("Seasonally Differenced Homicide Data")
hom.ts.differenced <- diff(hom.ts.differenced)
autoplot(hom.ts.differenced) + ggtitle("1st Difference of Seasonally Differenced Homicide Data")
autoplot(diff(hom.ts)) + ggtitle("1st Differenced Homicide Data")
kpss.test(hom.ts.differenced) #p-val = 0.1

adf.test(hom.ts.differenced) #p-val = 0.04299
#ADF test states stationarity

###########
# ETS Model
###########

#The first model is Holt Winters Additive Seasonality
hom.fit1 = hw(hom.train, seasonal="additive")
summary(hom.fit1)
autoplot(forecast(hom.fit1, h = 24), PI = FALSE) + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
accuracy(hom.fit1)

#hom.fit2 = hw(hom.train, seasonal="multiplicative")
#summary(hom.fit2)
#autoplot(forecast(hom.fit2, h = 24), PI = FALSE) + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
#accuracy(hom.fit2)
#   The multiplicative seasonality was not needed. So we only fit the additive method (which was the same as the ets())

#ETS Model
hom.ets <- ets(hom.train)
summary(hom.ets) 
#ETS(M,N,A)
autoplot(forecast(hom.ets, h=24), PI = FALSE) + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
checkresiduals(hom.ets)
accuracy(forecast(hom.ets, h = 24), hom.test)


#######
# Neural Network Method
#######

library(forecast)
hom.neuron <- nnetar(hom.train, P = 12)
#NNAR(12,1,6)[12]
hom.neuron %>% forecast(h = 24, PI = FALSE) %>% autoplot() + autolayer(hom.test) + xlab("Month/Year") + ylab("Monthly Homicides")
checkresiduals(hom.neuron)
accuracy(forecast(hom.neuron,h=24), hom.test)





