#' Descriptive Statistics
#'
#' Given a logistic regression model, generate the propensity weights for the
#' subjects (rows) in the data set and report the descriptive statistics for
#' each predictor in the regression model.
#'
#' @details
#' The regression model is expected to estimate the probability of an exposure
#' (Z = 1) given a set of predictors, X, i.e., Pr[Z = 1 | X].
#'
#' NOTE: for binary predictors coded as 0/1, such as male, the default action of
#' \code{dstats} will return a mean (sd), that is, if the \code{formula} for the
#' regression model is of the form
#' \code{y ~ male}, \code{dstats} will summarize the variable \code{male} as if it
#' was continous predictor. To get the summary of the proportion of males, and
#' females will be reported too, use a regression formula of the form \code{y ~
#' factor(male)}.
#'
#' See \code{\link{psweights}} for details on the propensity score based weights
#' used by \code{dstats}.
#' 
#' @param x a regression object with class \code{glm} or a \code{data.frame}
#' @param exposure_col column in \code{x} indicating the exposure (1) and
#' non-exposure (0).  Ignored if \code{x} is a regression model.
#' @param ps_col column in \code{x} indicating the propensity scores.  Ignored
#' if \code{x} is a regression model.
#' @param ... only used if \code{x} is a \code{data.frame}.  List the variables
#' in \code{x} to summarize. Default: summarize all columns.
#'
#' @return A \code{pstools_dstats} object which is a
#' \code{data.frame}, with descriptive statistics for each variable in and out
#' of the exposure group.
#'
#' @references
#' Li, Liang, and Tom Greene. "A weighting analogue to pair matching in
#' propensity score analysis." The international journal of biostatistics 9.2
#' (2013): 215-234.
#'
#' @seealso \code{\link[stats]{glm}} for fitting logistic regression models, or 
#' \code{\link[geepack]{geeglm}} for fitting GEEs.
#' \code{\link{psweights}} for the weights.
#'
#' @examples
#'
#' data(pride)
#' glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                      data = pride,
#'                      family = stats::binomial())
#' dstats(glmfit)
#'
#' @export
dstats <- function(x, ...) {
  UseMethod("dstats")
}

#' @export
dstats.data.frame <- function(x, exposure_col, ps_col, ...) {
  stop("To be built") 
}

#' @export
dstats.glm <- function(x, ...) {
  if (stats::family(x)$family != "binomial") {
    stop("expected binomial family regression model.")
  }

  # Exposure 1, no-exposure 0
  z <- stats::model.frame(x)[[1]]
  if (!all(z %in% c(0, 1))) { 
    stop("Expected Outcome vector is a 0/1 integer vector.")
  } 

  # propensity scores
  ps <- qwraps2::invlogit(stats::fitted(x))

  # get all the variable from the model fit and build the needed model matrix,
  # including the dumby variables for factors and character variables.
  av <- all.vars(stats::formula(x))

  mm <- lapply(attr(stats::terms(stats::formula(x)), "term.labels"),
               function(x) { stats::as.formula(paste("~", x, "+0")) }) 
  mm <- lapply(mm, stats::model.matrix, data = dplyr::select_(x$data, .dots = av))
  mm <- do.call(cbind, mm)

  # build a summary data.frame and add weights
  out <- dplyr::bind_cols(dplyr::data_frame(z, ps), dplyr::as_data_frame(mm)) 
  out <- psweights_(out, "z", "ps")

  # Calculate the mean and variance for (unadjusted) values and adjusted values
  # out <- tidyr::gather(out, key, value, -z, -w)
  out <- tidyr::gather_(out, key_col = "key", value_col = "value",
                        gather_cols = c(colnames(mm)))

  out <- dplyr::group_by_(out, .dots = list( ~ key, ~ z))
  out <- dplyr::summarize_(out, 
                           .dots = list("unadj_mean" = ~ mean(value), 
                                        "unadj_var"  = ~ var(value),
                                        "adj_mean_IPW" = lazyeval::interp(~ sum(w * value) / sum(w), w = as.name("psw_IPW")),
                                        "adj_var_IPW"  = lazyeval::interp(~ sum(w * (value - adj_mean_IPW)^2) / (sum(w) - 1), w = as.name("psw_IPW")), 
                                        "adj_mean_ACE_Exposed" = lazyeval::interp(~ sum(w * value) / sum(w), w = as.name("psw_ACE_Exposed")),
                                        "adj_var_ACE_Exposed"  = lazyeval::interp(~ sum(w * (value - adj_mean_ACE_Exposed)^2) / (sum(w) - 1), w = as.name("psw_ACE_Exposed")),
                                        "adj_mean_ACE_Unexposed" = lazyeval::interp(~ sum(w * value) / sum(w), w = as.name("psw_ACE_Unexposed")),
                                        "adj_var_ACE_Unexposed"  = lazyeval::interp(~ sum(w * (value - adj_mean_ACE_Unexposed)^2) / (sum(w) - 1), w = as.name("psw_ACE_Unexposed")),
                                        "adj_mean_ACE_MostInfo" = lazyeval::interp(~ sum(w * value) / sum(w), w = as.name("psw_ACE_MostInfo")),
                                        "adj_var_ACE_MostInfo"  = lazyeval::interp(~ sum(w * (value - adj_mean_ACE_MostInfo)^2) / (sum(w) - 1), w = as.name("psw_ACE_MostInfo")),
                                        "adj_mean_ACE_MWP" = lazyeval::interp(~ sum(w * value) / sum(w), w = as.name("psw_ACE_MWP")),
                                        "adj_var_ACE_MWP"  = lazyeval::interp(~ sum(w * (value - adj_mean_ACE_MWP)^2) / (sum(w) - 1), w = as.name("psw_ACE_MWP")) 
                                        )) 
  out <- dplyr::ungroup(out)

  attr(out, "continuous") <-  as.integer(unique(out$key) %in% names(stats::model.frame(x)))
  class(out) <- c("pstools_dstats", class(out))
  out
}

