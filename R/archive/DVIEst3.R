################
#  Functions
################
funcT <- function(tb, tc, to, t, alpha) {
    if (tb <= t && t <= tc) {
        log_funcT = alpha * log((t-tb)/(to- tb)) + alpha * (tc - to) / (to- tb) * log((tc-t)/(tc - to))
        output = exp(log_funcT)
    } else {
        output = 0
    }
    return(output)
}

funcP <- function(pb, pc, po, p, beta, pCritical) {

    # When pCritical is NULL
    if (is.null(pCritical)) {
        # Nakagawa et al. (1990)
        if (po <= p) {
            log_funcP = beta *  log((p-pb)/(po-pb)) + beta  * (pc - po ) / (po-pb)  * log((pc-p)/(pc-po))
            output = exp(log_funcP)
        } else {
            output = 1
        }
        # else
    } else {
        # This study using critical daylength (photoperiod)
        if (pCritical < p) {
            output = 0
        } else if (po <= p && p <= pCritical) {
            log_funcP = beta *  log((p-pb)/(po-pb)) + beta  * (pc - po) / (po-pb)  * log((pc-p)/(pc-po))
            output = exp(log_funcP)
        } else {
            output = 1
        }
    }
    return(output)
}

# Calculation of DVS at the given sowing and heading date.
calcDVS <- function(fixedParams, alpha, beta, g, sowingDate, headingDate, climateDf) {

    ######################################################################
    # initalization
    tb <- fixedParams$tb # base temperature in Celsius degree
    tc <- fixedParams$tc # ceiling temperature in Celsius degree
    to <- fixedParams$to # optimum temperature in Celsius degree
    po <- fixedParams$po # optimum photoperiod in hour (decimal)
    pb <- fixedParams$pb # base photoperiod in hour (decimal)
    pc <- fixedParams$pc # ceiling photoperiod in hour (decimal)
    pCritical <- fixedParams$pCritical # critical daylength in hour (decimal)
    dvs1 = 0.145 + 0.005 * g
    dvs2 = 0.345 + 0.005 * g
    # alpha is a parameter of temperature sensitivity to be estimated
    # beta is a parameter of photoperiod sensitivity to be estimated
    # g is a parameter of the earliness of flowering to be estimated
    # (scaled by the number of days from seedling emergence to flowering)
    #  under optimal photoperiod and temperature.

    # # Get the row indices for the sowing and heading date
    # sowingDate <- as.Date("2023-08-17")
    sowingDateRowIndex <- which(climateDf$Date == sowingDate)
    headingDateRowIndex <- which(climateDf$Date == headingDate)


    ###########################
    # At sowing date, observed DVS set as 0 and DVR is defined as f(t)/ g
    ###########################
    # temprature at sowing data
    t = climateDf[sowingDateRowIndex,"Temp"]
    DVSsowing = DVR = funcT(tb, tc, to, t, alpha) / g
    errorSowing = 0 - DVSsowing

    ###########################
    # At heading date DVS = 1
    ###########################
    DVSheading = 0
    for (i in sowingDateRowIndex:headingDateRowIndex) {
        date = climateDf[i,"Date"]
        t = climateDf[i,"Temp"]
        p = climateDf[i,"DayLength"]

        if (dvs1 <= DVSheading && DVSheading <= dvs2) {
            DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) * funcP(pb=pb, pc=pc, po=po, p=p, beta=beta, pCritical=pCritical) / g
        } else {
            DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) / g
        }

        DVSheading = DVSheading + DVR
        # for debug
        #print(paste0("Date: ", date, " Temp: ", t, " DayLength: ", p, " DVR :", DVR, " DVSheading : ", DVSheading ))
    }

    errorHeading = 1 - DVSheading
    return(c(errorSowing, errorHeading))
}

calcRMSE <- function(fixedParams, alpha, beta, g, headingDf, climateDf) {

    # initialization
    allErrors = c()
    # cat(sprintf("%s\t%s\t%s\t%s", "alpha", "beta", "g", "rmse"))

    for (i in 1:nrow(headingDf)) {
        sowingDate = headingDf[i, "Sowing"]
        headingDate = headingDf[i, "Heading"]
        errors = calcDVS(fixedParams=fixedParams, alpha=alpha, beta=beta, g=g, sowingDate=sowingDate, headingDate=headingDate, climateDf=climateDf)

        # Record errors as the vectors
        allErrors = c(allErrors, errors)
    }

    if (length(allErrors) > 0) {
        rmse <- sqrt(mean(allErrors^2))
    } else {
        rmse <- NULL
    }

    return(rmse)

}

# construnctoin of RMSE function conditioned by fixed parameters, and observed heading and climate data
makeRmseLossFunc = function(fixedParams, alpha, beta, g, headingDf, climateDf) {
    return(
        function(alpha, beta, g) {
            return(calcRMSE(fixedParams, alpha, beta, g, headingDf, climateDf))
        }
    )
}

getGAEstimate = function(ga_result) {
    solution = c(as.numeric(ga_result@solution[1,]), ga_result@fitnessValue * (-1))
    names(solution) = c("Estimated Alpha", "Estimated Beta", "Estimated G", "Minimum RSME")
    return(solution)
}


# For prediction
# Calculation of DVS at the given sowing and heading date.
predictHeadingDate <- function(fixedParams, alpha, beta, g, sowingDate, climateDf) {

    #######################################################################
    # To perform debug set debug = 1. Not to perform debug set debug = 0
    debug = 0

    ######################################################################
    # initalization
    tb <- fixedParams$tb # base temperature in Celsius degree
    tc <- fixedParams$tc # ceiling temperature in Celsius degree
    to <- fixedParams$to # optimum temperature in Celsius degree
    po <- fixedParams$po # optimum photoperiod in hour (decimal)
    pb <- fixedParams$pb # base photoperiod in hour (decimal)
    pc <- fixedParams$pc # ceiling photoperiod in hour (decimal)
    pCritical <- fixedParams$pCritical # critical daylength in hour (decimal)
    dvs1 = 0.145 + 0.005 * g
    dvs2 = 0.345 + 0.005 * g

    # # Get the row indices for the sowing date
    sowingDateRowIndex <- which(climateDf$Date == sowingDate)
    lastIndex <- length(climateDf$Date)

    #
    DVSheading = 0
    for (i in sowingDateRowIndex:lastIndex) {
        date = climateDf[i,"Date"]
        t = climateDf[i,"Temp"]
        p = climateDf[i,"DayLength"]

        if (dvs1 <= DVSheading && DVSheading <= dvs2) {
            DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) * funcP(pb=pb, pc=pc, po=po, p=p, beta=beta, pCritical=pCritical) / g
        } else {
            DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) / g
        }

        DVSheading = DVSheading + DVR

        if (DVSheading > 1) {
            #print(DVSheading)
            return(date)
        }
    }
    return(NA)
}

runPredictHeadingDate= function(estParam, headingDf2, fixedParams, accessionAnalyzed, climateDf, showDetail = T) {
    # prediction
    est_alpha = estParam[1]
    est_beta = estParam[2]
    est_g = estParam[3]

    # Apply the predictHeadingDate() function to each row and store the results in the variable predictedDates.
    predictedDates <- vector("numeric", nrow(headingDf2))
    for (i in 1:nrow(headingDf2)) {
        sowingDate <- headingDf2[i, "Sowing"]  # SowingDate 列を取得し、適切な変数に代入
        predictedDate <- predictHeadingDate(fixedParams=fixedParams, alpha=est_alpha, beta=est_beta, g=est_g, sowingDate=sowingDate, climateDf=climateDf)
        predictedDates[i] <- predictedDate
    }

    # Add a new column to the data frame."
    headingDf3=headingDf2
    headingDf3$PredictedHeadingDate <- predictedDates

    observedHeadingDates = as.Date(headingDf3[,"Heading"])
    sowingDates = as.Date(headingDf3[,"Sowing"])
    diffObs = as.integer(sowingDates - observedHeadingDates ) * (-1)
    predictedHeadingDates = as.Date(headingDf3[,"PredictedHeadingDate"])
    diffPred = as.integer(sowingDates - predictedHeadingDates ) * (-1)
    headingDf4 = cbind(headingDf3,diffObs, diffPred)

    if (showDetail == T) {
        print(headingDf4)
    }
    cat("cor: ",cor(diffObs, diffPred,use = "complete.obs"), "\n")
    plot(diffObs, diffPred, main = accessionAnalyzed)
}

estDVIparam = function(climateDf,headingDf, accessions, criticalDaylength, sowingDateGroups, fixedParams, maxiter) {

    for (accessionAnalyzed in accessions) {

        # Fixed parameters
        pCritical = NULL
        if (exists("criticalDaylength") && exists("accessionAnalyzed") && !is.null(criticalDaylength[[accessionAnalyzed]])) {
            pCritical = criticalDaylength[[accessionAnalyzed]]
        }
        fixedParams <- list(
            tb = fixedParams$tb,
            tc = fixedParams$tc,
            to = fixedParams$to,
            pb = fixedParams$pb,
            po = fixedParams$po,
            pc = fixedParams$pc,
            pCritical = pCritical
        )


        for (sowingDateForAnalysis in sowingDateGroups) {

            # Obtain the subset dataframe of the accessions.
            headingDf2 = na.omit(headingDf[headingDf$Name == accessionAnalyzed,])

            # using only subset data
            if (sowingDateForAnalysis == "all") {
                headingDf2 = headingDf2
            } else {
                condition <- headingDf2$Sowing %in% sowingDateForAnalysis
                headingDf2 <- headingDf2[condition, ]
            }


            RmseLossFunc = makeRmseLossFunc(alpha, beta, g, fixedParams=fixedParams, headingDf=headingDf2, climateDf=climateDf)

            ga_result <- ga(type = "real-valued", fitness = function(x) -RmseLossFunc(x[1], x[2], x[3]), lower = c(0,0,30), upper = c(40, 40, 150), maxiter=maxiter, popSize=50, monitor=F, parallel=10)

            # Result summary
            # summary(ga_result)
            cat(accessionAnalyzed, " ", sowingDateForAnalysis, "\n")
            estParam = getGAEstimate(ga_result)
            cat(getGAEstimate(ga_result))
            cat("\n")

            runPredictHeadingDate(estParam=estParam, headingDf2=headingDf2, fixedParams=fixedParams, accessionAnalyzed=accessionAnalyzed,climateDf=climateDf, showDetail = F )
        }
    }
}



################
#  User parameters
################

# Load data
myClimateFile = "C:/Users/yamagata/OneDrive - Kyushu University/Research/publication/論文作成中/continuousTransplanting/climateData.txt"
myHeadingFile = "C:/Users/yamagata/OneDrive - Kyushu University/Research/publication/論文作成中/continuousTransplanting/DVICalc.txt"
myClimateFile = "~/hdd6/satreps/dvi_model_R/climateData.txt"
myHeadingFile = "~/hdd6/satreps/dvi_model_R/DVICalc.txt"
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
set.seed(123)
system.time({
  estDVIparam(climateDf = myClimateDf, headingDf = myHeadingDf, accessions = myAccessions, criticalDaylength = myCriticalDaylength,
              sowingDateGroups = mySowingDateGroups, fixedParams = myFixedParams, maxiter = myMaxiter)
})

##################################################
# Do not use below
##################################################

for (accessionAnalyzed in accessions) {

    # Fixed parameters
    pCritical = criticalDaylength[[accessionAnalyzed]]
    fixedParams <- list(
        tb = 8, # base temperature in Celsius degree
        tc = 42, # ceiling temperature in Celsius degree
        to = 30,  # optimum temperature in Celsius degree
        pb = 0,  #  base photoperiod in hour (decimal)
        po = 10, # optimum  photoperiod in hour (decimal)
        pc = 24, # ceiling photoperiod in hour (decimal)
        pCritical = pCritical
    )

    for (sowingDateForAnalysis in sowingDateGroups) {

        # Obtain the subset dataframe of the accessions.
        headingDf2 = na.omit(headingDf[headingDf$Name == accessionAnalyzed,])

        # using only subset data
        if (sowingDateForAnalysis == "all") {
            headingDf2 = headingDf2
        } else {
            condition <- headingDf2$Sowing %in% sowingDateForAnalysis
            headingDf2 <- headingDf2[condition, ]
        }

        RmseLossFunc = makeRmseLossFunc(alpha, beta, g, fixedParams=fixedParams, headingDf=headingDf2, climateDf=climateDf)


        ga_result <- ga(type = "real-valued", fitness = function(x) -RmseLossFunc(x[1], x[2], x[3]), lower = c(0,0,30), upper = c(40, 40, 150), maxiter=500, popSize=50, monitor=F, parallel=T)

        # Result summary
        # summary(ga_result)
        cat(accessionAnalyzed, " ", sowingDateForAnalysis, "\n")
        estParam = getGAEstimate(ga_result)
        cat(getGAEstimate(ga_result))
        cat("\n")

        runPredictHeadingDate(estParam=estParam, headingDf2=headingDf2, fixedParams=fixedParams, accessionAnalyzed=accessionAnalyzed, showDetail = F )
    }
}


######################
##  練習
library(GA)
# 目的関数
objective_function <- function(x1, x2) {
    z = 2 * x1^2 + 2 * x1 - 4 + x2
  return(z)
}
# 遺伝的アルゴリズムの実行, library(GA)のga()関数はfitness関数の値が最大化するときのパラメータを推定する
# ため、注意する。目的関数に-1をかけて最小化問題にする。
ga_result <- ga(type = "real-valued", fitness = function(x) {-objective_function(x[1], x[2])},
 lower = c(-10,3.2), upper = c(10, 3.2))
# 最適解（最小値）のx
summary(ga_result)

## 高速化検討
tb = 8
tc = 42
to = 30
t = 21.3
pb = 8
po = 10
pc = 24
alpha = 10
beta = 10
p=16.3

output1 = (((t-tb)/(to-tb))^1 * (((tc-t)/(tc-to))^((tc-to)/(to-tb))) )^alpha
print(output1)

output2 = alpha * log((t-tb)/(to- tb)) + alpha * (tc - to) / (to- tb) * log((tc-t)/(tc - to))
exp(output2)

f1 = function(tb, tc, t0, t, alpha) {
output1 = (((t-tb)/(to-tb))^1 * (((tc-t)/(tc-to))^((tc-to)/(to-tb))) )^alpha
}
microbenchmark(f1(tb, tc, t0, t, alpha))

f2 = function(tb, tc, t0, t, alpha) {
output2 = alpha * log((t-tb)/(to- tb)) + alpha * (tc - to) / (to- tb) * log((tc-t)/(tc - to))
exp(output2)
}
microbenchmark(f2(tb, tc, t0, t, alpha))

# 結論  logをとった方が1.27倍高速だった

f3 = function(tb, tc, t0, t, alpha) {
output1 = (((t-tb)/(to-tb))^1 * (((tc-t)/(tc-to))^((tc-to)/(to-tb))) )^alpha
}
microbenchmark(f3(tb, tc, t0, t, alpha))

f3 = function(pb, pc, p0, p, beta) {
output3 = (((p-pb)/(po-pb))^1 * (((pc-p)/(pc-po))^((pc-po)/(po-pb))) )^beta
}
microbenchmark(f3(pb, pc, p0, p, beta))

f4 = function(pb, pc, p0, p, beta) {
output4 = beta *  log((p-pb)/(po-pb)) + beta  * (pc - po ) / (po-pb)  * log((pc-p)/(pc-po))
exp(output4)
}
microbenchmark(f4(pb, pc, p0, p, beta))

# 結論  logをとった方が1.25倍高速だった






