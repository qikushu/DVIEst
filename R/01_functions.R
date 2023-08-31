#'
#'
#'
#'
#'
#'
estDVIparam <- function(climateDf,
                        headingDf,
                        accessions,
                        criticalDaylength,
                        sowingDateGroups,
                        fixedParams,
                        maxiter) {

    for (accessionAnalyzed in accessions) {

        # Fixed parameters
        if (exists("criticalDaylength") && exists("accessionAnalyzed") && !is.null(criticalDaylength[[accessionAnalyzed]])) {
            fixedParams$pCritical <- criticalDaylength[[accessionAnalyzed]]

        } else {
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # See the comment in the code of funcP() to get why I set Inf here
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            fixedParams$pCritical <- Inf
        }

        for (sowingDateForAnalysis in sowingDateGroups) {

            # Obtain the subset dataframe of the accessions.
            headingDf2 = na.omit(headingDf[headingDf$Name == accessionAnalyzed, ])

            # using only subset data
            if (sowingDateForAnalysis == "all") {
                headingDf2 = headingDf2
            } else {
                condition <- headingDf2$Sowing %in% sowingDateForAnalysis
                headingDf2 <- headingDf2[condition, ]
            }


            RmseLossFunc = makeRmseLossFunc(alpha,
                                            beta,
                                            g,
                                            fixedParams = fixedParams,
                                            headingDf = headingDf2,
                                            climateDf = climateDf)

            ga_result <- ga(type = "real-valued",
                            fitness = function(x) -RmseLossFunc(x[1], x[2], x[3]),
                            lower = c(0, 0, 30),
                            upper = c(40, 40, 150),
                            maxiter = maxiter,
                            popSize = 50,
                            monitor = FALSE,
                            parallel = TRUE)

            cat(accessionAnalyzed, " ", sowingDateForAnalysis, "\n")
            estParam = getGAEstimate(ga_result)
            cat(getGAEstimate(ga_result))
            cat("\n")

            runPredictHeadingDate(estParam = estParam,
                                  headingDf2 = headingDf2,
                                  fixedParams = fixedParams,
                                  accessionAnalyzed = accessionAnalyzed,
                                  climateDf = climateDf,
                                  showDetail = FALSE)
        }
    }
}

# The beta function to describe the response to temperature

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The fixedParamter object should be given to the funcT function as is,
# not as separate scalars, for visibility of the code and also to avoid
# unexpected misspecification and modification of the values.
# Better to keep the capsulated object capsulated.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
funcT <- function(fp, t, alpha) {
    if (fp$tb <= t && t <= fp$tc) {
        log_funcT = alpha * log((t - fp$tb) / (fp$to - fp$tb)) +
            alpha * {(fp$tc - fp$to) / (fp$to - fp$tb)} *
            log((fp$tc - t) / (fp$tc - fp$to))
        output = exp(log_funcT)

    } else {
        output = 0
    }

    return(output)
}

# The beta function to describe the response to photoperiod

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The fixedParamter object should be given to the funcP function as is,
# not as separate scalars, for visibility of the code and also to avoid
# unexpected misspecification and modification of the values.
# Better to keep the capsulated object capsulated.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
funcP <- function(fp, p, beta) {

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # If you set pCritical to Inf when pCritical can be ignored,
    # fp$pCritical < p never be TRUE and p <= fp$pCritical is always TRUE.
    # Then, fp$po <= p && p <= fp$pCritical is the same as fp$po <= p.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (is.null(fp$pCritical)) {
        # When pCritical is NULL,
        # use the equation described in Nakagawa et al. (1990)
        if (fp$po <= p) {
            log_funcP = beta * log((p - fp$pb) / (fp$po - fp$pb)) +
                beta * (fp$pc - fp$po) / (fp$po - fp$pb) *
                log((fp$pc - p) / (fp$pc - fp$po))
            output = exp(log_funcP)

        } else {
            output = 1
        }

    } else {
        if (fp$pCritical < p) {
            output = 0

        } else if (fp$po <= p && p <= fp$pCritical) {
            log_funcP = beta * log((p - fp$pb) / (fp$po - fp$pb)) +
                beta * {(fp$pc - fp$po) / (fp$po - fp$pb)} *
                log((fp$pc - p) / (fp$pc - fp$po))
            output = exp(log_funcP)

        } else {
            output = 1
        }
    }
    return(output)
}

# Calculate DVS at the given sowing and heading date.
calcDVS <- function(fixedParams,
                    alpha,
                    beta,
                    g,
                    sowingDate,
                    headingDate,
                    climateDf) {
    # alpha is a parameter of temperature sensitivity to be estimated
    # beta is a parameter of photoperiod sensitivity to be estimated
    # g is a parameter of the earliness of flowering to be estimated
    # (scaled by the number of days from seedling emergence to flowering)
    #  under optimal photoperiod and temperature.

    # Initalization
    # tb <- fixedParams$tb # base temperature in Celsius degree
    # tc <- fixedParams$tc # ceiling temperature in Celsius degree
    # to <- fixedParams$to # optimum temperature in Celsius degree
    # po <- fixedParams$po # optimum photoperiod in hour (decimal)
    # pb <- fixedParams$pb # base photoperiod in hour (decimal)
    # pc <- fixedParams$pc # ceiling photoperiod in hour (decimal)
    # pCritical <- fixedParams$pCritical # critical daylength in hour (decimal)
    dvs1 <- 0.145 + 0.005 * g
    dvs2 <- 0.345 + 0.005 * g

    # # Get the row indices for the sowing and heading date
    # sowingDate <- as.Date("2023-08-17")
    s_date <- which(climateDf$Date == sowingDate)
    h_date <- which(climateDf$Date == headingDate)

    # At sowing date, observed DVS set as 0 and DVR is defined as f(t)/ g
    # temprature at sowing data
    t <- climateDf$Temp[s_date]
    dvs_sowing <- funcT(fp = fixedParams,
                       t = climateDf$Temp[s_date],
                       alpha = alpha)
    error_sowing <- 0 - dvs_sowing / g

    # At heading date DVS = 1
    dvs_heading <- 0
    for (i in s_date:h_date) {
        t <- climateDf$Temp[i]
        p <- climateDf$DayLength[i]

        DVR <- funcT(fp = fixedParams, cl = climateDf, alpha = alpha)

        if (dvs1 <= dvs_heading && dvs_heading <= dvs2) {
            DVR <- DVR * funcP(fp = fixedParams,
                               p = p,
                               beta = beta)
        }

        dvs_heading = dvs_heading + DVR / g
        # for debug
        # print(paste0("Date: ", climateDf$Date[i], " Temp: ", t,
        # " DayLength: ", p, " DVR :", DVR, " DVSheading : ", DVSheading ))
    }

    error_heading = 1 - dvs_heading
    return(c(error_sowing, error_heading))
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

