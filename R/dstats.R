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
dstats.formula <- function(formula, data, ps) { 
  lhs_vars <- all.vars(lazyeval::f_lhs(formula))
  rhs_vars <- all.vars(lazyeval::f_rhs(formula))

  exposure_col <- lhs_vars
  ps_col       <- deparse(substitute(ps))

  mm <- lapply(c(lhs_vars, rhs_vars, ps_col),
               function(x) { stats::as.formula(paste("~", x, "+0")) }) 
  
  mm <- lapply(mm, stats::model.matrix, data = dplyr::select_(data, .dots = c(lhs_vars, rhs_vars, ps_col)))
  mm <- do.call(cbind, mm)

  out <- psweights_(as.data.frame(mm), exposure_col, ps_col)

  out <- tidyr::gather_(out, key_col = "key", value_col = "value",
                        gather_cols = setdiff(colnames(mm), c(exposure_col, ps_col))
                        )

  out <- dplyr::group_by_(out, .dots = list( ~ key, lazyeval::interp( ~ z, z = as.name(exposure_col))))

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

  attr(out, "continuous")   <- as.integer(unique(out$key) %in% names(data))
  attr(out, "exposure_col") <- exposure_col
  attr(out, "ps_col")       <- ps_col
  class(out) <- c("pstools_dstats", class(out))
  out 
}

#' @export
dstats.glm <- function(x, ...) {
  if (stats::family(x)$family != "binomial") {
    stop("expected binomial family regression model.")
  }

  # propensity scores
  ps <- qwraps2::invlogit(stats::fitted(x))

  dstats.formula(formula(x), cbind(x[["data"]], ps = ps), ps) 
} 
