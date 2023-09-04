gaFixedParam <- function(lower = c(0, 0, 30),
                         upper = c(40, 40, 150),
                         maxiter = 100,
                         popSize = 50,
                         parallel = 1){
    out <- list(lower = lower, 
                upper = upper, 
                maxiter = maxiter,
                popSize = popSize,
                parallel = parallel)
    class(out) <- c(class(out), "GaParam")
    return(out)
}

dviFixedParam <- function(tb = 8,
                          tc = 42,
                          to = 30,
                          pb = 0,  
                          po = 10, 
                          pc = 24){
    out <- list(tb = tb, tc = tc, to = to, pb = pb, po = po, pc = pc)
    class(out) <- c(class(out), "DviParam")
    return(out)
}

buildDVImodel <- function(date,
                          day_len,
                          temp,
                          sowing,
                          heading,
                          critical,
                          acc){
    climate <- list(date = date, day_len = day_len, temp = temp)
    
    out <- list(climate = climate, 
                sowing = sowing,
                heading = heading,
                critical = critical, 
                acc = acc)
    class(out) <- c(class(out), "DviModel")
    return(out)
}

#'
#'
#'
#'
#'
#'

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Object names were simplified to avoid typo and keep better visibility 
# of codes throughout the entire script.
# Lowercase letters connected with underscore(s) are preferred for object 
# names in the naming convention of R.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# To make arguments simple, functions to make input objects for 
# estDVIparam() were defined as shown above.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# estDVIparam = function(climateDf,headingDf, accessions, criticalDaylength, sowingDateGroups, fixedParams, maxiter) {
# Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
estDVIparam <- function(object,
                        fixed_param = dviFixedParam(),
                        ga_param = gaFixedParam()) {
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # Validate input objects
    stopifnot(inherits(x = object, "DviModel"))
    stopifnot(inherits(x = fixed_param, "DviParam"))
    stopifnot(inherits(x = ga_param, "GaParam"))
    
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Loops over accessions and sowing date groups should be 
    # handled outside of this function.
    # The estDVIparam() funcrion takes one dataset including climate information, 
    # sowing, and heading dates for one accession and estimate parameters 
    # for a DVI model for the accession.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # for (accessionAnalyzed in accessions) {
    
    # Fixed parameters
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The argument accessions and criticalDaylength are not allowed to 
    # be empty. Thus, exists("accessionAnalyzed") and 
    # exists("criticalDaylength") should be always TRUE.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Original >>>>>>>>>>>
    # if (exists("criticalDaylength") && exists("accessionAnalyzed") && !is.null(criticalDaylength[[accessionAnalyzed]])) {
    # modified >>>>>>>>>>>
    if (!is.null(object$critical)) {
        fixed_param$pCritical <- object$critical
        # <<<<<<<<<<<<<<<<<<<<
        
    } else {
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # See the comment in the code of funcP() to get why I set Inf here
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fixed_param$pCritical <- Inf
    }
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # No need to redefine fixedParams
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fixedParams <- list(
    #     tb = fixedParams$tb,
    #     tc = fixedParams$tc,
    #     to = fixedParams$to,
    #     pb = fixedParams$pb,
    #     po = fixedParams$po,
    #     pc = fixedParams$pc,
    #     pCritical = pCritical
    # )
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Delete an unnecessary loop inside this function
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # for (sowingDateForAnalysis in sowingDateGroups) {
    
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Subsetting of input data should be done before inputting it here.
    # Prepare a function to subset the input if necessary.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # # Obtain the subset dataframe of the accessions.
    # heading_subset <- na.omit(heading[heading$Name == acc, ])
    
    # using only subset data
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # # using only subset data
    # if (sowingDateForAnalysis == "all") {
    #     headingDf2 = headingDf2
    # } else {
    #     condition <- headingDf2$Sowing %in% sowingDateForAnalysis
    #     headingDf2 <- headingDf2[condition, ]
    # }
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>s>>>>>>>>>>>>>>
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # RmseLossFunc = makeRmseLossFunc(alpha, beta, g, fixedParams=fixedParams, headingDf=headingDf2, climateDf=climateDf)
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    RmseLossFunc <- .makeRmseLossFunc(alpha, 
                                      beta,
                                      g,
                                      fixed_param = fixed_param,
                                      model = object)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    ga_result <- ga(type = "real-valued",
                    fitness = function(x) -RmseLossFunc(x[1], x[2], x[3]),
                    lower = ga_param$lower,
                    upper = ga_param$upper,
                    maxiter = ga_param$maxiter,
                    popSize = ga_param$popSize,
                    monitor = FALSE,
                    parallel = ga_param$parallel)
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Let this function return an DviParam class object that contains
    # the output of ga(), the fixed parameters and input data used in 
    # parameter estimation to allow users predict heading date without 
    # specify the fixed parameters and input data again for 
    # heading-date prediction.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # cat(accessionAnalyzed, " ", sowingDateForAnalysis, "\n")
    # estParam = getGAEstimate(ga_result)
    # cat(getGAEstimate(ga_result))
    # cat("\n")
    # runPredictHeadingDate(estParam=estParam, headingDf2=headingDf2, fixedParams=fixedParams, accessionAnalyzed=accessionAnalyzed,climateDf=climateDf, showDetail = F )
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    out <- list(ga = ga_result, 
                model = object,
                fixed_param = fixed_param,
                ga_param = ga_param)
    
    class(out) <- c(class(out), "DviEst")
    return(out)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

# The beta function to describe the response to temperature

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The fixedParamter object should be given to the funcT function as is,
# not as separate scalars, for visibility of the code and also to avoid
# unexpected misspecification and modification of the values.
# Better to keep the capsulated object capsulated.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Original >>>>>>>>>>>>>>>>>
# funcT <- function(tb, tc, to, t, alpha) {
#     if (tb <= t && t <= tc) {
#         log_funcT = alpha * log((t-tb)/(to- tb)) + alpha * (tc - to) / (to- tb) * log((tc-t)/(tc - to))
# Modified >>>>>>>>>>>>>>>>>
.funcT <- function(fp, t, alpha) {
    if (fp$tb <= t && t <= fp$tc) {
        log_funcT <- alpha * log((t - fp$tb) / (fp$to - fp$tb)) +
            alpha * {(fp$tc - fp$to) / (fp$to - fp$tb)} *
            log((fp$tc - t) / (fp$tc - fp$to))
        # <<<<<<<<<<<<<<<<<<<<<<<<<<
        output <- exp(log_funcT)
        
    } else {
        output <- 0
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

# Original >>>>>>>>>>>>>
# funcP <- function(pb, pc, po, p, beta, pCritical) {
# Modified >>>>>>>>>>>>>
.funcP <- function(fp, p, beta) {
    # <<<<<<<<<<<<<<<<<<<<<<
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # If you set pCritical to Inf when pCritical can be ignored,
    # fp$pCritical < p never be TRUE and p <= fp$pCritical is always TRUE.
    # Then, fp$po <= p && p <= fp$pCritical is the same as fp$po <= p.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # # When pCritical is NULL
    # if (is.null(pCritical)) {
    #     # Nakagawa et al. (1990)
    #     if (po <= p) {
    #         log_funcP = beta *  log((p-pb)/(po-pb)) + beta  * (pc - po ) / (po-pb)  * log((pc-p)/(pc-po))
    #         output <- exp(log_funcP)
    #         
    #     } else {
    #         output <- 1
    #     }
    
    # } else {
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # if (pCritical < p) {
    #     output = 0
    #     
    # } else if (po <= p && p <= pCritical) {
    #     log_funcP = beta *  log((p-pb)/(po-pb)) + beta  * (pc - po) / (po-pb)  * log((pc-p)/(pc-po))
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (fp$pCritical < p) {
        output <- 0
        
    } else if (fp$po <= p && p <= fp$pCritical) {
        log_funcP <- beta * log((p - fp$pb) / (fp$po - fp$pb)) +
            beta * {(fp$pc - fp$po) / (fp$po - fp$pb)} *
            log((fp$pc - p) / (fp$pc - fp$po))
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        output <- exp(log_funcP)
        
    } else {
        output <- 1
    }
    # }
    return(output)
}

# Calculate DVS at the given sowing and heading date.
.calcDVS <- function(index,
                     fixed_param,
                     alpha,
                     beta,
                     g,
                     climate,
                     sowing,
                     heading) {
    # alpha is a parameter of temperature sensitivity to be estimated
    # beta is a parameter of photoperiod sensitivity to be estimated
    # g is a parameter of the earliness of flowering to be estimated
    # (scaled by the number of days from seedling emergence to flowering)
    #  under optimal photoperiod and temperature.
    
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The fixedParamter object should be given to the funcP function as is,
    # not as separate scalars, for visibility of the code and also to avoid
    # unexpected misspecification and modification of the values.
    # Better to keep the capsulated object capsulated.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # sowingDateRowIndex <- which(climateDf$Date == sowingDate)
    # headingDateRowIndex <- which(climateDf$Date == headingDate)
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    s_date_i <- which(climate$date == sowing[index])
    h_date_i <- which(climate$date == heading[index])
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # At sowing date, observed DVS set as 0 and DVR is defined as f(t)/ g
    # temprature at sowing data
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # t = climateDf[sowingDateRowIndex,"Temp"]
    # DVSsowing = DVR = funcT(tb, tc, to, t, alpha) / g
    # errorSowing = 0 - DVSsowing
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    dvs_sowing <- .funcT(fp = fixed_param,
                         t = climate$temp[s_date_i],
                         alpha = alpha)
    error_sowing <- 0 - dvs_sowing / g
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # At heading date DVS = 1
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # DVSheading = 0
    # for (i in sowingDateRowIndex:headingDateRowIndex) {
    #     date = climateDf[i,"Date"]
    #     t = climateDf[i,"Temp"]
    #     p = climateDf[i,"DayLength"]
    #     
    #     if (dvs1 <= DVSheading && DVSheading <= dvs2) {
    #         DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) * funcP(pb=pb, pc=pc, po=po, p=p, beta=beta, pCritical=pCritical) / g
    #     } else {
    #         DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) / g
    #     }
    #     
    #     DVSheading = DVSheading + DVR
    #     # for debug
    #     #print(paste0("Date: ", date, " Temp: ", t, " DayLength: ", p, " DVR :", DVR, " DVSheading : ", DVSheading ))
    # }
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    dvs_heading <- 0
    for (i in seq(s_date_i, h_date_i)) {
        dvr <- .funcT(fp = fixed_param,
                      t = climate$temp[i],
                      alpha = alpha)
        
        if (dvs1 <= dvs_heading && dvs_heading <= dvs2) {
            dvr <- dvr * .funcP(fp = fixed_param,
                                p = climate$day_len[i],
                                beta = beta)
        }
        
        dvs_heading = dvs_heading + dvr / g
    }
    
    error_heading <- 1 - dvs_heading
    return(c(error_sowing, error_heading))
}

.calcRMSE <- function(fixed_param, alpha, beta, g, model) {

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Somtimes vapply() is faster than for loop and improve the code visibility
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    # # initialization
    # allErrors = c()
    # # cat(sprintf("%s\t%s\t%s\t%s", "alpha", "beta", "g", "rmse"))
    # 
    # for (i in 1:nrow(headingDf)) {
    #     sowingDate = headingDf[i, "Sowing"]
    #     headingDate = headingDf[i, "Heading"]
    #     errors = calcDVS(fixedParams=fixedParams, alpha=alpha, beta=beta, g=g, sowingDate=sowingDate, headingDate=headingDate, climateDf=climateDf)
    #     
    #     # Record errors as the vectors
    #     allErrors = c(allErrors, errors)
    # }
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    out <- vapply(X = seq_along(model$sowing),
                  FUN.VALUE = numeric(2), 
                  FUN = .calcDVS, 
                  fixed_param = fixed_param,
                  alpha = alpha,
                  beta = beta,
                  g = g, 
                  climate = model$climate,
                  heading = model$heading,
                  sowing = model$sowing)
    out <- as.vector(out)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    if (length(out) > 0) {
        out <- sqrt(mean(out^2))
    } else {
        out <- NULL
    }
    
    return(out)
}

# construnctoin of RMSE function conditioned by fixed parameters, and observed heading and climate data
.makeRmseLossFunc = function(fixed_param, alpha, beta, g, model) {
    return(
        function(alpha, beta, g) {
            return(.calcRMSE(fixed_param, alpha, beta, g, model))
        }
    )
}

.getGAEstimate <- function(ga_result) {
    solution <- c(as.numeric(ga_result@solution[1, ]), -ga_result@fitnessValue)
    names(solution) <- c("Estimated Alpha",
                         "Estimated Beta",
                         "Estimated G",
                         "Minimum RSME")
    return(solution)
}


# For prediction
# Calculation of DVS at the given sowing and heading date.
.predictHeadingDate <- function(sowing_date,
                                fixed_param,
                                alpha,
                                beta,
                                g, 
                                climate) {
    
    #######################################################################
    # To perform debug set debug = 1. Not to perform debug set debug = 0
    # debug <- 0
    
    ######################################################################
    # # initalization
    # tb <- fixedParams$tb # base temperature in Celsius degree
    # tc <- fixedParams$tc # ceiling temperature in Celsius degree
    # to <- fixedParams$to # optimum temperature in Celsius degree
    # po <- fixedParams$po # optimum photoperiod in hour (decimal)
    # pb <- fixedParams$pb # base photoperiod in hour (decimal)
    # pc <- fixedParams$pc # ceiling photoperiod in hour (decimal)
    # pCritical <- fixedParams$pCritical # critical daylength in hour (decimal)
    dvs1 <- 0.145 + 0.005 * g
    dvs2 <- 0.345 + 0.005 * g
    
    # # Get the row indices for the sowing date
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # sowingDateRowIndex <- which(climateDf$Date == sowingDate)
    # lastIndex <- length(climateDf$Date)
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    s_date_i <- which(climate$date == sowing_date)
    last_i <- length(climate$date)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # DVSheading = 0
    # for (i in sowingDateRowIndex:lastIndex) {
    #     date = climateDf[i,"Date"]
    #     t = climateDf[i,"Temp"]
    #     p = climateDf[i,"DayLength"]
    #     
    #     if (dvs1 <= DVSheading && DVSheading <= dvs2) {
    #         DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) * funcP(pb=pb, pc=pc, po=po, p=p, beta=beta, pCritical=pCritical) / g
    #     } else {
    #         DVR = funcT(tb=tb, tc=tc, to=to, t=t, alpha=alpha) / g
    #     }
    #     
    #     DVSheading = DVSheading + DVR
    #     
    #     if (DVSheading > 1) {
    #         #print(DVSheading)
    #         return(date)
    #     }
    # }
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    dvs_heading <- 0
    for (i in seq(s_date_i, last_i)) {
        dvr <- .funcT(fp = fixed_param,
                      t = climate$temp[i],
                      alpha = alpha)
        
        if (dvs1 <= dvs_heading && dvs_heading <= dvs2) {
            dvr <- dvr * .funcP(fp = fixed_param,
                                p = climate$day_len[i],
                                beta = beta)
        }
        
        dvs_heading <- dvs_heading + dvr / g
        
        if (dvs_heading > 1) {
            return(climate$date[i])
        }
    }
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return(NA)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Define a function to predict the heading date based on the given model and 
# inputs.
# This function is defined as a method of the S3 generic function predict().
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
predict.DviEst <- function(object, new_model = NULL) {
    
    # Validate the input
    stopifnot(inherits(x = object, "DviEst"))
    
    # Get the parameters estimated using GA
    est_param <- .getGAEstimate(ga_result = object$ga)
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Allow this function to predict the heading date based on 
    # newly supplied climate information and a sowing date, but not
    # based on the input data used for the parameter estimation using GA.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(is.null(new_model)){
        model <- object$model
    } else {
        stopifnot(inherits(x = new_model, "DviModel"))
        model <- new_model
    }
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # predictedDates <- vector("numeric", nrow(headingDf2))
    # for (i in 1:nrow(headingDf2)) {
    #     sowingDate <- headingDf2[i, "Sowing"]  # SowingDate 列を取得し、適切な変数に代入
    #     predictedDate <- predictHeadingDate(fixedParams=fixedParams, alpha=est_alpha, beta=est_beta, g=est_g, sowingDate=sowingDate, climateDf=climateDf)
    #     predictedDates[i] <- predictedDate
    # }
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    hd_pred <- vapply(X = model$sowing,
                      FUN.VALUE = character(1),
                      FUN = .predictHeadingDate,
                      fixed_param = object$fixed_param,
                      alpha = est_param[1],
                      beta = est_param[2], 
                      g = est_param[3],
                      climate = model$climate)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Let this function output a predicted heading date and the number of 
    # predicted days to heading.
    # Only if NOT new_model was specified, this function returns 
    # the observed headng date, the number of observed days to heading, and
    # the correlation between the prediction and observation.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     # Add a new column to the data frame."
#     headingDf3=headingDf2
#     headingDf3$PredictedHeadingDate <- predictedDates
#     
#     observedHeadingDates = as.Date(headingDf3[,"Heading"])
#     sowingDates = as.Date(headingDf3[,"Sowing"])
#     diffObs = as.integer(sowingDates - observedHeadingDates ) * (-1)
#     predictedHeadingDates = as.Date(headingDf3[,"PredictedHeadingDate"])
#     diffPred = as.integer(sowingDates - predictedHeadingDates ) * (-1)
#     headingDf4 = cbind(headingDf3,diffObs, diffPred)
#     
#     if (showDetail == T) {
#         print(headingDf4)
#     }
#     cat("cor: ",cor(diffObs, diffPred,use = "complete.obs"), "\n")
#     plot(diffObs, diffPred, main = accessionAnalyzed)
    # Modified >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sowing_dates <- as.Date(model$sowing)
    if(is.null(model$heading)){
        dth_obs <- NA
    } else {
        hd_obs <- as.Date(model$heading)
        dth_obs <- as.integer(hd_obs - sowing_dates)
    }
    hd_pred <- as.Date(hd_pred)
    dth_pred <- as.integer(hd_pred - sowing_dates)
    
    object$estimated_param <- est_param
    object$prediction <- data.frame(HD_Pred = hd_pred, 
                                    DTH_Obs = dth_obs, 
                                    DTH_Pred = dth_pred)
    
    if(is.na(dth_obs)){
        object$correlation <- NA
        
    } else {
        object$correlation <- cor(dth_obs, dth_pred, use = "complete.obs")
    }
    
    class(object) <- c(class(object), "DviEst")
    return(object)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Define a function to show simple summary of the DviEst object
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
summary.DviEst <- function(x){
    cat("Accession: ", x$model$acc, "\n")
    cat("Estimated parameters:\n")
    print(x$estimated_param)
    cat("Correlation: ", x$correlation)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Define a function to draw a scatter plot to show correlation between 
# the predicted and observed values.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plotDviEst<- function(object){
    # Validate the input.
    stopifnot(inherits(x = object, "DviEst"))
    
    # Check if predicted values exist.
    if(is.null(object$prediction)){
        stop("No predicted values. Run predict.DviEst().")
    }
    if(all(is.na(object$prediction$DTH_Obs))){
        stop("Observed data was not provided.")
    }
    
    # Make a ggplot object with minimum aesthetics and output the object
    # to allow users additional modification of the plot.
    p <- ggplot(data.frame(Predicted = object$prediction$DTH_Pred,
                           Observed = object$prediction$DTH_Obs)) +
        geom_point(aes(x = Observed, y = Predicted)) +
        labs(title = object$model$acc)
    return(p)
}