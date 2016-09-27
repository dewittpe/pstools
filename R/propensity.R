#' Propensity Summary
#'
#' Given a logistic regression model generate a the propensity weights for the
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
#' @param fit a regression object such as 
#'
#' @return A \code{propensity_summary} object which is a
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
#' @export
propensity_summary <- function(fit) {
  UseMethod("propensity_summary")
}

#' @export
propensity_summary.glm <- function(fit) {
  if (stats::family(fit)$family != "binomial") {
    stop("expected binomial family regreesion model.")
  }

  # Exposure 1, no-exposure 0
  z <- stats::model.frame(fit)[[1]]
  if (!all(z %in% c(0, 1))) { 
    stop("Expected Outcome vector is a 0/1 integer vector.")
  } 

  # propensity scores
  ps <- qwraps2::invlogit(fitted(fit))


}



propensity_summary.default <- function(fit) {
  z <- as.integer(stats::model.frame(fit)[[1]])


  w <- qwraps2::invlogit(fitted(fit))

  av <- all.vars(formula(fit))

  mm <-
    lapply(attr(terms(formula(fit)), "term.labels"),
           function(x) {
             as.formula(paste("~", x, "+0"))
           }) %>%
    lapply(., model.matrix, data = dplyr::select_(fit$data, .dots = av)) %>%
    do.call(cbind, .)

  out <-
    dplyr::bind_cols(dplyr::data_frame(z, w), as_data_frame(mm)) %>%

    dplyr::rowwise() %>%
    dplyr::mutate(w = min(c(w, 1-w)) / (z*w + (1-z)*(1-w))) %>%
    # dplyr::mutate(w = 1/(z*w + (1-z)*(1-w))) %>%
    # dplyr::mutate(w = z + (1-z)*w/(1-w)) %>%
    dplyr::ungroup() %>%

    tidyr::gather(key, value, -z, -w) %>%

    dplyr::group_by(key, z) %>%
    dplyr::summarize(mean(value), var(value),
                     `mean(adjvalue)` = sum(w * value) / sum(w),
                     `var(adjvalue)`  = sum(w * (value - `mean(adjvalue)`)^2) / (sum(w) - 1)
                     ) %>%
    dplyr::ungroup() %>%

    {
    dplyr::left_join(dplyr::filter(., z == 0),
                     dplyr::filter(., z == 1),
                     by = "key",
                     suffix = c(".z0", ".z1"))
    } %>%

    dplyr::select(-z.z0, -z.z1) %>%

    dplyr::mutate(`unadj std diff` = 100 * (`mean(value).z1` - `mean(value).z0`)       / sqrt((`var(value).z1` + `var(value).z0`)/2),
                  `adj std diff`   = 100 * (`mean(adjvalue).z1` - `mean(adjvalue).z0`) / sqrt((`var(adjvalue).z1` + `var(adjvalue).z0`)/2),
                  `continuousvar` = as.integer(key %in% names(stats::model.frame(fit))))  %>%

    dplyr::rowwise() %>%

    dplyr::summarize(key,
                     `unadj z1` = if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(value).z1`)), " (", qwraps2::frmt(sqrt(`var(value).z1`)), ")") } else {qwraps2::frmt(mean(`mean(value).z1`)*100, 1)},
                     `unadj z0` = if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(value).z0`)), " (", qwraps2::frmt(sqrt(`var(value).z0`)), ")") } else {qwraps2::frmt(mean(`mean(value).z0`)*100, 1)},
                     `unadj std diff` = qwraps2::frmt(`unadj std diff`, 1),
                     `adj z1` = if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(adjvalue).z1`)), " (", qwraps2::frmt(sqrt(`var(adjvalue).z1`)), ")") } else {qwraps2::frmt(mean(`mean(adjvalue).z1`)*100, 1)},
                     `adj z0` = if (continuousvar) { paste0(qwraps2::frmt(mean(`mean(adjvalue).z0`)), " (", qwraps2::frmt(sqrt(`var(adjvalue).z0`)), ")") } else {qwraps2::frmt(mean(`mean(adjvalue).z0`)*100, 1)},
                     `adj std diff` = qwraps2::frmt(`adj std diff`, 1)
                     )

  class(out) <- c("propensity_summary", class(out))
  out
}

################################################################################
# plot.propensity_summary
#
# An plotting method for propensity_summary objects
plot.propensity_summary <- function(x, y, ...) {

  plotting_data <-
    dplyr::select_(x, .dots = list("key", "`unadj std diff`", "`adj std diff`"))

  plotting_data <-
    tidyr::gather_(plotting_data,
                   key_col = 'variable',
                   value_col = 'stddiff',
                   gather_cols = list("unadj std diff", "adj std diff"))
  plotting_data$stddiff <- as.numeric(plotting_data$stddiff)

  plotting_data$variable <- factor(plotting_data$variable)

  plotting_data$key <- factor(plotting_data$key)
  plotting_data$key <- factor(plotting_data$key, levels = rev(levels(plotting_data$key)))

  ggplot2::ggplot(plotting_data) +
  ggplot2::theme_bw() +
  ggplot2::aes_string(x = "stddiff", y = "key") +
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = "variable", shape = "variable")) +
  ggplot2::geom_vline(xintercept = 0) +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())
}



