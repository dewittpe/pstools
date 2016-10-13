#' Propensity Summary
#'
#' Given a logistic regression model, generate the propensity weights for the
#' subjects (rows) in the data set.
#'
#' The regression model is expected to estimate the probability of an exposure
#' (Z = 1) given a set of predictors, X, i.e., Pr[Z = 1 | X].
#'
#' NOTE: for binary predictors coded as 0/1, such as male, the default action of
#' propensity_summary will return a mean (sd), that is, if
#' propensity_summary(geeglm(y ~ male)) will summarize the variable male as if it
#' was continous predictor. To get the summary for the proportion of males, and
#' females will be reported too, use propensity_summary(geeglm(y ~
#' factor(male))).
#'
#' @details
#' \code{weight_method} controls the weight given to each row of the data set.
#' Let ps = Pr[Z = 1 | X].  The default setting is \code{weight_method = 1}.
#' 
#' \describe{
#'   \item{\code{weight_method = 1}}{\code{w =  min(c(ps, 1-ps)) / (z*ps + (1-z)*(1-ps))}}
#'   \item{\code{weight_method = 2}}{\code{w = 1/(z*ps + (1-z)*(1-ps))}}
#'   \item{\code{weight_method = 3}}{\code{w = z + (1-z)*ps/(1-ps))}}
#' }
#' 
#' @param fit a regression object such as 
#' @param weight_method one of three methods for defining the weight for each
#' row of the data set.  See Details.
#'
#' @return A \code{pstools_propensity} object which is a
#' \code{data.frame} with columns summarizing each varaible used as a predictor in the
#' propensity model. The function will determine if each variable is a
#' categorical or continous variable. if the variable is continuous the mean (sd)
#' is returned.  If the variable is categorical the the percentage of each level
#' is returned.  Standard differences are reported as well. Unadjusted and
#' Matched Weight (see Li and Greene (2012)) values are reported.
#'
#' @references
#' Li, Liang, and Tom Greene. "A weighting analogue to pair matching in
#' propensity score analysis." The international journal of biostatistics 9.2
#' (2013): 215-234.
#'
#' @seealso \code{\link[stats]{glm}} for fitting logistic regression models, or 
#' \code{\link[geepack]{geeglm}} for fitting GEEs.
#'
#' @examples
#'
#' data(pride)
#' glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                      data = pride,
#'                      family = stats::binomial())
#' propensity(glmfit)
#' summary(propensity(glmfit))
#' plot(propensity(glmfit))
#'
#' @export
propensity <- function(fit, weight_method = 1) {
  UseMethod("propensity")
}

#' @export
propensity.glm <- function(fit, weight_method = 1) {
  if (!(weight_method %in% 1:3)) { 
    stop("weight_method needs to be 1, 2, or 3.")
  }

  if (stats::family(fit)$family != "binomial") {
    stop("expected binomial family regression model.")
  }

  # Exposure 1, no-exposure 0
  z <- stats::model.frame(fit)[[1]]
  if (!all(z %in% c(0, 1))) { 
    stop("Expected Outcome vector is a 0/1 integer vector.")
  } 

  # propensity scores
  ps <- qwraps2::invlogit(stats::fitted(fit))

  # get all the variable from the model fit and build the needed model matrix,
  # including the dumby variables for factors and character variables.
  av <- all.vars(stats::formula(fit))

  mm <- lapply(attr(stats::terms(stats::formula(fit)), "term.labels"),
               function(x) { stats::as.formula(paste("~", x, "+0")) }) 
  mm <- lapply(mm, stats::model.matrix, data = dplyr::select_(fit$data, .dots = av))
  mm <- do.call(cbind, mm)

  # build a summary data.frame and add weights
  out <- dplyr::bind_cols(dplyr::data_frame(z, ps), dplyr::as_data_frame(mm))
  out <- dplyr::rowwise(out)

  if (weight_method == 1) { 
    out <- dplyr::mutate_(out, .dots = list("w" = ~ min(c(ps, 1-ps)) / (z*ps + (1-z)*(1-ps))))
  } else if (weight_method == 2) { 
    out <- dplyr::mutate_(out, .dots = list("w" = ~  1/(z*ps + (1-z)*(1-ps))))
  } else {
    out <- dplyr::mutate_(out, .dots = list("w" = ~ z + (1-z)*ps/(1-ps)))
  } 

  out <- dplyr::ungroup(out)

  # Calculate the mean and variance for (unadjusted) values and adjusted values
  # out <- tidyr::gather(out, key, value, -z, -w)
  out <- tidyr::gather_(out, key_col = "key", value_col = "value",
                        gather_cols = c(colnames(mm)))

  out <- dplyr::group_by_(out, .dots = list( ~ key, ~ z))
  out <- dplyr::summarize_(out, 
                           .dots = list(~ mean(value), 
                                        ~ var(value),
                                        `mean(adjvalue)` = ~ sum(w * value) / sum(w),
                                        `var(adjvalue)`  = ~ sum(w * (value - `mean(adjvalue)`)^2) / (sum(w) - 1)))
  out <- dplyr::ungroup(out)

  # Split the data.frame by exposure and column bind
  out <- dplyr::left_join(dplyr::filter_(out, .dots = list( ~ z == 0)),
                          dplyr::filter_(out, .dots = list( ~ z == 1)),
                          by = "key",
                          suffix = c(".z0", ".z1"))
  out <- dplyr::select_(out, .dots = list( ~ -z.z0, ~ -z.z1))

  # calculate the standardized difference (%)
  # add a column to not if the variable is continous or a factor
  out <- dplyr::mutate_(out, 
                        .dots = list(
                       `unadj std diff (%)` = ~ 100 * (`mean(value).z1` - `mean(value).z0`)       / sqrt((`var(value).z1` + `var(value).z0`)/2),
                       `adj std diff (%)`   = ~ 100 * (`mean(adjvalue).z1` - `mean(adjvalue).z0`) / sqrt((`var(adjvalue).z1` + `var(adjvalue).z0`)/2),
                       `continuousvar` =  ~ as.integer(key %in% names(stats::model.frame(fit)))))

  class(out) <- c("pstools_propensity", class(out))
  out
}

#' @export
summary.pstools_propensity <- function(object, ...) {
  dplyr::summarize_(dplyr::rowwise(object),
                    .dots = list(
                                 key = ~ key,
                                 `unadj z1` = ~ if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(value).z1`)), " (", qwraps2::frmt(sqrt(`var(value).z1`)), ")") } else {qwraps2::frmt(mean(`mean(value).z1`)*100, 1)},
                                 `unadj z0` = ~ if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(value).z0`)), " (", qwraps2::frmt(sqrt(`var(value).z0`)), ")") } else {qwraps2::frmt(mean(`mean(value).z0`)*100, 1)},
                                 `unadj std diff (%)` = ~ qwraps2::frmt(`unadj std diff (%)`, 1),
                                 `adj z1` = ~ if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(adjvalue).z1`)), " (", qwraps2::frmt(sqrt(`var(adjvalue).z1`)), ")") } else {qwraps2::frmt(mean(`mean(adjvalue).z1`)*100, 1)},
                                 `adj z0` = ~ if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(adjvalue).z0`)), " (", qwraps2::frmt(sqrt(`var(adjvalue).z0`)), ")") } else {qwraps2::frmt(mean(`mean(adjvalue).z0`)*100, 1)},
                                 `adj std diff (%)` = ~ qwraps2::frmt(`adj std diff (%)`, 1)
                                 )) 
}

#' @export
plot.pstools_propensity <- function(x, y, ...) { 
  plotting_data <-
    dplyr::select_(x, .dots = list("key", "`unadj std diff (%)`", "`adj std diff (%)`"))

  plotting_data <-
    tidyr::gather_(plotting_data,
                   key_col = 'variable',
                   value_col = 'stddiff (%)',
                   gather_cols = list("unadj std diff (%)", "adj std diff (%)"))

  plotting_data$variable <- factor(plotting_data$variable)

  plotting_data$key <- factor(plotting_data$key)
  plotting_data$key <- factor(plotting_data$key, levels = rev(levels(plotting_data$key)))

  ggplot2::ggplot(plotting_data) +
  ggplot2::theme_bw() +
  ggplot2::aes_string(x = "`stddiff (%)`", y = "key") +
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = "variable", shape = "variable")) +
  ggplot2::geom_vline(xintercept = 0) +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())
}

