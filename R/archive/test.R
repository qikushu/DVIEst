################
#  User parameters
################

# Load data
myClimateFile = "./inst/extdata/climateData.txt"
myHeadingFile = "./inst/extdata/DVICalc.txt"
myClimateDf = read.table(myClimateFile, head=T)
myHeadingDf = read.table(myHeadingFile, head=T)

myFixedParams <- list(
  tb = 8, # base temperature in Celsius degree
  tc = 42, # ceiling temperature in Celsius degree
  to = 30,  # optimum temperature in Celsius degree
  pb = 0,  #  base photoperiod in hour (decimal)
  po = 10, # optimum  photoperiod in hour (decimal)
  pc = 24 # ceiling photoperiod in hour (decimal)
)


mySowingDateGroups = list(
  c("2019-02-01","2019-03-01","2019-04-01", "2019-05-01"),
  c("2019-03-01","2019-04-01","2019-05-01", "2019-06-01"),
  c("2019-04-01","2019-05-01","2019-06-01", "2019-07-01"),
  c("2019-05-01","2019-06-01","2019-07-01", "2019-08-01"),
  c("2019-06-01","2019-07-01","2019-08-01", "2019-09-01"),
  c("2019-07-01","2019-08-01","2019-09-01", "2019-10-01"),
  c("2019-08-01","2019-09-01","2019-10-01", "2019-11-01"),
  c("2019-09-01","2019-10-01","2019-11-01", "2019-12-01"),
  c("2019-10-01","2019-11-01","2019-12-01", "2020-01-01")
  #c("2019-02-01","2019-03-01","2019-04-01", "2019-05-01","2019-06-01", "2019-07-01")
  #c("2019-07-01","2019-08-01","2019-09-01", "2019-10-01","2019-11-01", "2019-12-01")
)
mySowingDateGroups = list("all")

myAccessions = c("NKTP","IR24", "STK")
myAccessions = c("PSH", "NKG", "IMY","NKTP","MSMKK", "IR24", "STK")

myAccessions = c("NKTP")
# Fixed parameters
myCriticalDaylength <- list(PSH=11.86667, NKG=11.90000, IMY=12.16667, NKTP=12.75000, MSMKK=13.08333, SBPS=13.33333, IR24=13.31667, STK=13.31667)


################
#  Let's estimate
################
#Use library GA
library(GA)
myMaxiter = 100
myHeadingDf <-
estDVIparam(climateDf = myClimateDf, headingDf = myHeadingDf, accessions = myAccessions, criticalDaylength = myCriticalDaylength, sowingDateGroups = mySowingDateGroups, fixedParams = myFixedParams, maxiter = myMaxiter)
# NKTP   all
# 3.415394 6.938856 96.63959 0.1177181
# cor:  0.6001106

source("R/01_functions.R")
set.seed(123)
heading <- na.omit(subset(myHeadingDf, subset = Name == "NKTP"))
model <- buildDVImodel(climate = myClimateDf,
                       sowing = heading$Sowing,
                       heading = heading$Heading,
                       critical = myCriticalDaylength[["NKTP"]],
                       acc = "NKTP")
out <- estDVIparam(object = model, ga_param = gaFixedParam(parallel = 2))
out <- predict(object = out)
summary(out)
plotDviEst(out)
