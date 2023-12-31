#########################################################
# Functions for the parameter estimation of a DVI model #
#########################################################


#' Set genetic algorithm parameters
#'
#' Create a parameter object to control the parameter estimation by a genetic
#' algorithm.
#'
#' @param lower a vector of length equal to the decision variables providing the
#'  lower bounds of the search space in case of real-valued or permutation
#'  encoded optimizations.
#' @param upper a vector of length equal to the decision variables providing the
#'  upper bounds of the search space in case of real-valued or permutation
#'   encoded optimizations.
#' @param maxiter the maximum number of iterations to run before the GA search
#' is halted.
#' @param popSize the population size.
#' @param parallel An optional argument which allows to specify if the Genetic
#' Algorithm should be run sequentially or in parallel. See [GA::ga()] for more
#' details.
#'
#' @return A GaParam object that can be specified to the `ga_param` argument of
#' [estDVIparam()]
#'
#' @seealso GA::ga()
#'
#' @export
#'
gaParam <- function(lower = c(0, 0, 30),
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


#' Set fixed parameters of DVI models
#'
#' Create a parameter object to specify fixed parameters of a DVI model.
#'
#' @param tb a numeric value to specify the base temperature.
#' @param to a numeric value to specify the optimal temperature.
#' @param tc a numeric value to specify the ceiling temperature.
#' @param pb a numeric value to specify the base photoperiod.
#' @param po a numeric value to specify the optimal photoperiod.
#' @param pc a numeric value to specify the ceiling photoperiod.
#'
#' @return A DviParam object that can be specified to the `fixed_param`
#' argument of [estDVIparam()]
#'
#' @export
#'
dviFixedParam <- function(tb = 8,
                          to = 30,
                          tc = 42,
                          pb = 0,
                          po = 10,
                          pc = 24){
    out <- list(tb = tb, to = to, tc = tc, pb = pb, po = po, pc = pc)
    class(out) <- c(class(out), "DviParam")
    return(out)
}


#' Create an input object
#'
#' Create an input object for the parameter estimation of the DVI model using
#' a genetic algorithm.
#'
#' @param date a character vector to indicate the dates on which day lengths
#' and temperature were measured. The length should match the lengths of
#' `day_len` and `temp`.
#' @param day_len a numeric vector to indicate the dates on which day lengths
#' and temperature were measured. The length should match the lengths of `date`
#' and `temp`.
#' @param temp a numeric vector to indicate the dates on which day lengths and
#' temperature were measured. The length should match the lengths of `date`
#' and `day_len`.
#' @param sowing a numeric value to indicate the sowing date.
#' @param heading a numeric value to indicate the heading date.
#' @param critical a numeric value to indicate the critical day length.
#' @param acc a character indicates sample ID.
#' @param verbose No message if FALSE.
#'
#' @details Dates should be given in the 'YYYY-MM-DD' format.
#'
#' @return A DviModel object that can be specified to the `model` argument of
#' [estDVIparam()]
#'
#' @export
#'
buildDVImodel <- function(date,
                          day_len,
                          temp,
                          sowing,
                          heading = NULL,
                          critical = NULL,
                          acc = "",
                          verbose = FALSE){

    # Validate input values
    if(length(date) != length(day_len)){
        stop("The lengths of date and day_len should be the same.")
    }
    if(length(date) != length(temp)){
        stop("The lengths of date and temp should be the same.")
    }
    is_date <- try(as.Date(date))
    if(inherits(x = is_date, what = "try-error")){
        stop("date should be a character vector of dates in a valid format")
    }

    # Check NA values in input vectors and omit them
    incl_na <- is.na(date) | is.na(day_len) | is.na(temp)
    if(any(incl_na)){
        if(verbose){
            message("NA was found in the given climate information. \n",
                    "The data point with NA was removed.")
        }
    }

    climate <- list(date = date[!incl_na],
                    day_len = day_len[!incl_na],
                    temp = temp[!incl_na])

    # Validate sowing date values
    is_date <- try(as.Date(sowing))
    if(inherits(x = is_date, what = "try-error")){
        stop("sowing should be a character vector of dates in a valid format")
    }

    # Validate heading date values if it is specified.
    # Check NA values in input vectors and omit them
    if(is.null(heading)){
        incl_na <- is.na(sowing)

    } else {
        is_date <- try(as.Date(heading))
        if(inherits(x = is_date, what = "try-error")){
            stop("heading should be a character vector of dates in a valid format")
        }

        incl_na <- is.na(sowing) | is.na(heading)
        heading <- heading[!incl_na]
    }
    sowing <- sowing[!incl_na]

    if(any(incl_na)){
        if(verbose){
            message("NA was found in the given sowing and/or heading dates. \n",
                    "The data point with NA was removed.")
        }
    }

    # Check the given sowing and heading dates are in the range of the dates
    # given as `date` of climate information.
    sowing_valid <- sapply(sowing, function(x){
        return(any(date <= x))
    })

    if(!is.null(heading)){
        heading_valid <- sapply(heading, function(x){
            return(any(date >= x))
        })
    } else {
        heading_valid <- TRUE
    }

    if(any(!sowing_valid & !heading_valid)){
        stop("Sowing and/or heading dates contains date(s) out of the range of",
             " the dates on which climate information is available.")
    }

    out <- list(climate = climate,
                sowing = sowing,
                heading = heading,
                critical = critical,
                acc = acc)
    class(out) <- c(class(out), "DviModel")
    return(out)
}


#' DVI model parameter estimation
#'
#' Estimate the parameters of the DVI model using a genetic algorithm.
#'
#' @param object a DviModel object created by [buildDVImodel()].
#' @param fixed_param a DviParam object created by [dviFixedParam()].
#' @param ga_parma a GaParam object created by [gaParam()].
#'
#' @details Dates should be given in the 'YYYY-MM-DD' format.
#'
#' @return A DviEst object.
#'
#' @seealso [predict.DviEst()], [summary.DviEst()] and [plotDviEst]
#'
#' @importFrom GA ga
#'
#' @export
#'
#' @examples
#' fixed_param <- dviFixedParam()
#' ga_param <- gaParam()
#' date <- seq.Date(from = as.Date("2023-06-01"),
#'                  to = as.Date("2023-08-31"),
#'                  by = "day")
#' day_len <- c(seq(12, 13.5, length.out = 30),
#'              seq(13.5, 12, length.out = length(date) - 30))
#' temp = c(sample(x = 25:30, size = 30, replace = TRUE),
#'          sample(x = 25:35, size = 31, replace = TRUE),
#'          sample(x = 30:35, size = 31, replace = TRUE))
#' dvi_model <- buildDVImodel(date = date,
#'                            day_len = day_len,
#'                            temp = temp,
#'                            sowing = "2023-06-05",
#'                            heading = "2023-08-20",
#'                            critical = 13,
#'                            acc = "test")
#'
#' dvi_est <- estDVIparam(object = dvi_model,
#'                        fixed_param = fixed_param,
#'                        ga_param = ga_param)
#'

estDVIparam <- function(object,
                        fixed_param = dviFixedParam(),
                        ga_param = gaParam()) {

    # Validate input objects
    stopifnot(inherits(x = object, "DviModel"))
    stopifnot(inherits(x = fixed_param, "DviParam"))
    stopifnot(inherits(x = ga_param, "GaParam"))

    if(is.null(object$heading)){
        stop("The parameter estimation requires heading date information in",
             " the DviModel object. Remake a DviModel object with specifying ",
             "dates to the `heading` arugment in buildDVImode().")
    }
    if (!is.null(object$critical)) {
        fixed_param$critical <- object$critical

    } else {
        fixed_param$critical <- Inf
    }
    RmseLossFunc <- .makeRmseLossFunc(alpha,
                                      beta,
                                      g,
                                      fixed_param = fixed_param,
                                      model = object)

    ga_result <- ga(type = "real-valued",
                    fitness = function(x) -RmseLossFunc(x[1], x[2], x[3]),
                    lower = ga_param$lower,
                    upper = ga_param$upper,
                    maxiter = ga_param$maxiter,
                    popSize = ga_param$popSize,
                    monitor = FALSE,
                    parallel = ga_param$parallel)

    # Get the parameters estimated using GA
    est_param <- .getGAEstimate(ga_result = ga_result)

    out <- list(ga = ga_result,
                model = object,
                fixed_param = fixed_param,
                ga_param = ga_param,
                estimated_param = est_param)

    class(out) <- c(class(out), "DviEst")
    return(out)
}

# The beta function to describe the response to temperature
.funcT <- function(fp, t, alpha) {
    output <- rep(0, length(t))
    in_range <- fp$tb <= t & t <= fp$tc
    log_funcT <- alpha * log((t - fp$tb) / (fp$to - fp$tb)) +
        alpha * {(fp$tc - fp$to) / (fp$to - fp$tb)} *
        log((fp$tc - t) / (fp$tc - fp$to))
    output[in_range] <- exp(log_funcT[in_range])
    # <<<<<<<<<<<<<<<<<<<<<<<<<<
    return(output)
}

# The beta function to describe the response to photoperiod
.funcP <- function(fp, p, beta) {
    output <- rep(1, length(p))
    out_range <- fp$critical < p
    in_range <- fp$po <= p & p <= fp$critical
    log_func_p <- beta * log((p - fp$pb) / (fp$po - fp$pb)) +
        beta * {(fp$pc - fp$po) / (fp$po - fp$pb)} *
        log((fp$pc - p) / (fp$pc - fp$po))
    output[in_range] <- exp(log_func_p[in_range])
    output[out_range] <- 0
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

    dvs1 <- 0.145 + 0.005 * g
    dvs2 <- 0.345 + 0.005 * g

    # Get the row indices for the sowing and heading date
    s_date_i <- which(climate$date == sowing[index])
    h_date_i <- which(climate$date == heading[index])

    # At sowing date, observed DVS set as 0 and DVR is defined as f(t)/ g
    # temprature at sowing data

    dvs_sowing <- .funcT(fp = fixed_param,
                         t = climate$temp[s_date_i],
                         alpha = alpha)
    error_sowing <- 0 - dvs_sowing / g

    func_t_out <- .funcT(fp = fixed_param,
                         t = climate$temp[s_date_i:h_date_i],
                         alpha = alpha)
    dvr <- func_t_out / g

    dvs_heading <- cumsum(dvr)

    in_range <- dvs1 <= dvs_heading
    if(sum(in_range) > 1){
        in_range[which(in_range)[1]] <- FALSE
        func_p_out <- .funcP(fp = fixed_param,
                             p = climate$day_len[s_date_i:h_date_i][in_range],
                             beta = beta)
        dvr[in_range] <- func_t_out[in_range] * func_p_out / g
    }

    dvs_heading <- cumsum(dvr)

    out_range <- dvs2 < dvs_heading
    if(sum(out_range) > 1){
        out_range[which(out_range)[1]] <- FALSE
        dvr[out_range] <- func_t_out[out_range] / g
    }

    dvs_heading <- cumsum(dvr)

    error_heading <- 1 - tail(dvs_heading, 1)
    return(c(error_sowing, error_heading))
}

.calcRMSE <- function(fixed_param, alpha, beta, g, model) {
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


#' DVI model parameter estimation
#'
#' Estimate the parameters of the DVI model using a genetic algorithm.
#'
#' @param object a DviModel object created by [buildDVImodel()].
#' @param fixed_param a DviParam object created by [dviFixedParam()].
#' @param ga_parma a GaParam object created by [gaParam()].
#'
#' @details Dates should be given in the 'YYYY-MM-DD' format.
#'
#' @return A DviEst object.
#'
#' @seealso [predict.DviEst()], [summary.DviEst()] and [plotDviEst]
#'
#' @export
#'
#' @examples
#' fixed_param <- dviFixedParam()
#' ga_param <- gaParam()
#' date <- seq.Date(from = as.Date("2023-06-01"),
#'                  to = as.Date("2023-08-31"),
#'                  by = "day")
#' day_len <- c(seq(12, 13.5, length.out = 30),
#'              seq(13.5, 12, length.out = length(date) - 30))
#' temp = c(sample(x = 25:30, size = 30, replace = TRUE),
#'          sample(x = 25:35, size = 31, replace = TRUE),
#'          sample(x = 30:35, size = 31, replace = TRUE))
#' dvi_model <- buildDVImodel(date = date,
#'                            day_len = day_len,
#'                            temp = temp,
#'                            sowing = "2023-06-05",
#'                            heading = "2023-08-20",
#'                            critical = 13,
#'                            acc = "test")
#'
#' dvi_est <- estDVIparam(object = dvi_model,
#'                        fixed_param = fixed_param,
#'                        ga_param = ga_param)
#'

################################################################################

##############################################
# Functions for the heading date predication #
##############################################


#' Heading date prediction
#'
#' Predict heading dates based on sowing dates and the DVI model with the
#' estimated parameters.
#'
#' @param object a DviEst object created by [estDVIparam()].
#' @param new_mode a DviModel object created by [buildDVImodel()].
#'
#' @details The predicated heading date (HD_Pred) and the days to heading
#' (DTH_Pred) will be returned as a data.frame in the prediction slot of the
#' output DviEst object. If the observed heading dates are available in the
#' input DviEst object, these values also returned at the DTH_Obs column of
#' the data.frame in the prediction slot. The predicted heading dates will be
#' compared with the observed heading dates recorded in the given DviEst object,
#' and the correlation coefficient between these will be returned in the
#' correlation slot of the DviEst object.
#' If a DviModle object was specified to `new_model`, only the predicated
#' heading date (HD_Pred) and days to heading (DTH_Pred) will be returned in the
#' prediction slot.
#'
#' @return A DviEst object.
#'
#' @seealso [predict.DviEst()], [summary.DviEst()] and [plotDviEst]
#'
#' @export
#'
#' @examples
#' fixed_param <- dviFixedParam()
#' ga_param <- gaParam()
#' date <- seq.Date(from = as.Date("2023-06-01"),
#'                  to = as.Date("2023-08-31"),
#'                  by = "day")
#' day_len <- c(seq(12, 13.5, length.out = 30),
#'              seq(13.5, 12, length.out = length(date) - 30))
#' temp = c(sample(x = 25:30, size = 30, replace = TRUE),
#'          sample(x = 25:35, size = 31, replace = TRUE),
#'          sample(x = 30:35, size = 31, replace = TRUE))
#' dvi_model <- buildDVImodel(date = date,
#'                            day_len = day_len,
#'                            temp = temp,
#'                            sowing = "2023-06-05",
#'                            heading = "2023-08-20",
#'                            critical = 13,
#'                            acc = "test")
#'
#' dvi_est <- estDVIparam(object = dvi_model,
#'                        fixed_param = fixed_param,
#'                        ga_param = ga_param)
#'
#' dvi_est <- predict(object = dvi_est)
#'

predict.DviEst <- function(object, new_model = NULL) {

    # Validate the input
    stopifnot(inherits(x = object, "DviEst"))

    if(is.null(new_model)){
        model <- object$model
    } else {
        stopifnot(inherits(x = new_model, "DviModel"))
        model <- new_model
    }

    if (!is.null(object$model$critical)) {
        object$fixed_param$critical <- object$model$critical
    } else {
        object$fixed_param$critical <- Inf
    }

    hd_pred <- vapply(X = model$sowing,
                      FUN.VALUE = character(1),
                      FUN = .predictHeadingDate,
                      fixed_param = object$fixed_param,
                      alpha = object$estimated_param[1],
                      beta = object$estimated_param[2],
                      g = object$estimated_param[3],
                      climate = model$climate)

    sowing_dates <- as.Date(model$sowing)
    if(is.null(model$heading)){
        dth_obs <- NA
    } else {
        hd_obs <- as.Date(model$heading)
        dth_obs <- as.integer(hd_obs - sowing_dates)
    }
    hd_pred <- as.Date(hd_pred)
    dth_pred <- as.integer(hd_pred - sowing_dates)

    object$prediction <- data.frame(HD_Pred = hd_pred,
                                    DTH_Obs = dth_obs,
                                    DTH_Pred = dth_pred)

    if(is.na(dth_obs[1])){
        object$correlation <- NA

    } else {
        object$correlation <- cor(dth_obs, dth_pred, use = "complete.obs")
        object$lm <- lm(formula = dth_pred ~ dth_obs)
    }

    class(object) <- c(class(object), "DviEst")
    return(object)
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

    dvs1 <- 0.145 + 0.005 * g
    dvs2 <- 0.345 + 0.005 * g

    # Get the row indices for the sowing date
    s_date_i <- which(climate$date == sowing_date)
    last_i <- length(climate$date)

    func_t_out <- .funcT(fp = fixed_param,
                         t = climate$temp[s_date_i:last_i],
                         alpha = alpha)
    dvr <- func_t_out / g

    dvs_heading <- cumsum(dvr)

    in_range <- dvs1 <= dvs_heading
    if(sum(in_range) > 1){
        in_range[which(in_range)[1]] <- FALSE
        func_p_out <- .funcP(fp = fixed_param,
                             p = climate$day_len[s_date_i:last_i][in_range],
                             beta = beta)
        dvr[in_range] <- func_t_out[in_range] * func_p_out / g
    }

    dvs_heading <- cumsum(dvr)

    out_range <- dvs2 < dvs_heading
    if(sum(out_range) > 1){
        out_range[which(out_range)[1]] <- FALSE
        dvr[out_range] <- func_t_out[out_range] / g
    }

    dvs_heading <- cumsum(dvr)

    heading_date_i <- min(which(dvs_heading > 1))

    if(length(heading_date_i) == 0){
        return(NA)

    } else {
        return(climate$date[s_date_i:last_i][heading_date_i])
    }
}

################################################################################

#################
# Plot function #
#################


#' Plot scatter plot of correlation
#'
#' Plot a scatter plot to visualize correlation between the observed heading
#' dates and predicated heading dates.
#'
#' @param object a DviEst object created by [estDVIparam()].
#'
#' @return A ggplot2 object.
#'
#' @seealso [predict.DviEst()], [summary.DviEst()] and [plotDviEst]
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' fixed_param <- dviFixedParam()
#' ga_param <- gaParam()
#' date <- seq.Date(from = as.Date("2023-06-01"),
#'                  to = as.Date("2023-08-31"),
#'                  by = "day")
#' day_len <- c(seq(12, 13.5, length.out = 30),
#'              seq(13.5, 12, length.out = length(date) - 30))
#' temp = c(sample(x = 25:30, size = 30, replace = TRUE),
#'          sample(x = 25:35, size = 31, replace = TRUE),
#'          sample(x = 30:35, size = 31, replace = TRUE))
#' dvi_model <- buildDVImodel(date = date,
#'                            day_len = day_len,
#'                            temp = temp,
#'                            sowing = "2023-06-05",
#'                            heading = "2023-08-20",
#'                            critical = 13,
#'                            acc = "test")
#'
#' dvi_est <- estDVIparam(object = dvi_model,
#'                        fixed_param = fixed_param,
#'                        ga_param = ga_param)
#'
#' dvi_est <- predict(object = dvi_est)
#'
#' p <- plotDviEst(object = dvi_est)
#'
#' print(p)
#'

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

    r <- signif(x = object$correlation, digits = 3)
    lm <- summary(object$lm)
    r2 <- signif(x = lm$r.squared, digits = 3)
    acc <- object$model$acc

    lm_pred <- predict(object = object$lm)
    lm_pred_edge <- c(which.min(lm_pred), which.max(lm_pred))

    # Make a ggplot object with minimum aesthetics and output the object
    # to allow users additional modification of the plot.
    p <- ggplot(data.frame(Predicted = object$prediction$DTH_Pred,
                           Observed = object$prediction$DTH_Obs)) +
        geom_point(aes(x = Observed, y = Predicted)) +
        geom_segment(aes(x = object$lm$model$dth_obs[lm_pred_edge[1]],
                         y = lm_pred[lm_pred_edge[1]],
                         xend = object$lm$model$dth_obs[lm_pred_edge[2]],
                         yend = lm_pred[lm_pred_edge[2]]), color = "magenta") +
        labs(title = bquote(.(acc) ~ ", " ~ r ~ "=" ~ .(r) ~ ", " ~
                 " Adjusted " ~ R^2 ~ "=" ~ .(r2)))
    return(p)
}


################################################################################

#####################
# Utility functions #
#####################

#' Show summary of a DviModel object
#'
#' @param x a DviModel object created by [buildDVImodel()].
#'
#' @return Invisibly return the input object itself.
#'
#' @export
#'
summary.DviModel <- function(x){
    cat("Accession: ", x$acc, "\n\n")
    if(length(x$sowing) > 10){
        cat("Sowing date: \n", x$sowing[1:10],
            "... (additional", length(x$sowing) - 10, "data)\n\n")

    } else {
        cat("Sowing date: \n", x$sowing, "\n\n")
    }
    if(length(x$heading) > 10){
        cat("Heading date: \n", x$heading[1:10],
            "... (additional", length(x$heading) - 10, "data)\n\n")

    } else {
        cat("Heading date: \n", x$heading, "\n\n")
    }
    cat("Critical day length: ", x$critical, "\n")
    invisible(x)
}

#' Show summary of a DviEst object
#'
#' @param x a DviEst object created by [estDVIparam()].
#'
#' @return Invisibly return a data.frame of the summary.
#'
#' @export
#'

summary.DviEst <- function(x){
    cat("Accession: ", x$model$acc, "\n")
    cat("Estimated parameters:\n")
    print(x$estimated_param)
    if(!is.null(x$correlation)){
        cat("Correlation: ", x$correlation)
        lm_summary <- summary(x$lm)
        cat("\nR_squared: ", lm_summary$r.squared)
        cat("\nAdjusted R_squared: ", lm_summary$adj.r.squared)
        invisible(data.frame(acc = x$model$acc,
                             x$estimated_param,
                             cor = x$correlation,
                             r_squared = lm_summary$r.squared,
                             adj_r_squared =lm_summary$adj.r.squared))
    } else {
        invisible(data.frame(acc = x$model$acc,
                             x$estimated_param))
    }
}


#' Set parameters of a DVI model
#'
#' Set arbitrary values as estimated parameters of a DVI model.
#'
#' @param object a DviModel object created by [buildDVImodel()].
#' @param fixed_param a DviParam object created by [dviFixedParam()].
#' @param alpha a numeric value of the temperature-sensitivity coefficient.
#' @param beta a numeric value of the photoperiod-sensitivity coefficient.
#' @param g a numeric value of the earliness of flowering under optimal
#' photoperiod and temperature.
#'
#' @return a DviEst object.
#'
#' @export
#'
setDviParam <- function(object, fixed_param, alpha, beta, g){
    stopifnot(inherits(x = object, "DviModel"))
    stopifnot(inherits(x = fixed_param, "DviParam"))

    if(!is.numeric(alpha) | length(alpha) != 1){
        stop("alpha should be a numeric value.")
    }
    if(!is.numeric(beta) | length(beta) != 1){
        stop("beta should be a numeric value.")
    }
    if(!is.numeric(g) | length(g) != 1){
        stop("g should be a numeric value.")
    }
    est_param <- c(alpha, beta, g)
    names(est_param) <- c("Estimated Alpha",
                          "Estimated Beta",
                          "Estimated G")
    if(is.null(object$fixed_param$critical)){
        fixed_param$critical <- NULL

    } else {
        fixed_param$critical <- object$fixed_param$critical
    }

    out <- list(model = object,
                fixed_param = fixed_param,
                estimated_param = est_param)

    class(out) <- c(class(out), "DviEst")
    return(out)
}


#' Get estimated parameters of a DVI model
#'
#' Extracted esetimted parameter values from the DviEst object.
#'
#' @param object a DviEstl object created by [estDVIparam()].
#'
#' @return a numeric vector of estimated alpha, beta, and g.
#'
#' @export
#'
getEstParam <- function(object){
    stopifnot(inherits(x = object, "DviEst"))
    return(object$estimated_param[1:3])
}

#' Get predicted heading dates
#'
#' Extracted predicted heading dates and days-to-heading data
#' from the DviEst object.
#'
#' @param object a DviEstl object created by [estDVIparam()].
#'
#' @return a data.frame extracted from the prediction slot of the DviEst object.
#'
#' @export
#'
getPred <- function(object){
    stopifnot(inherits(x = object, "DviEst"))
    return(object$prediction)
}



#' Get the correlation coefficient
#'
#' Extracted the correlation coefficient from the DviEst object.
#'
#' @param object a DviEstl object created by [estDVIparam()].
#'
#' @return a numeric value of the correlation coefficient.
#'
#' @export
#'
getCor <- function(object){
    stopifnot(inherits(x = object, "DviEst"))
    return(object$correlation)
}


#' Get the linear model
#'
#' Extracted the linear model that represents the relationship between
#' the observed and predicted heading dates.
#'
#' @param object a DviEstl object created by [estDVIparam()].
#'
#' @return a lm class object.
#'
#' @export
#'
getLM <- function(object){
    stopifnot(inherits(x = object, "DviEst"))
    return(object$lm)
}
