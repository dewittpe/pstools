#' Propensity Score Weights
#'
#' Five different matching weights based on an exposure and propensity score.
#'
#' @details
#' Let p denote the propensity score for a subject.  The five weights are as follows.
#' 
#' \tabular{llcc}{
#'                                                                          \tab                          \tab Exposed        \tab Non-Exposed        \cr
#' Average causal effect in study population                                \tab \code{psw_IPW}           \tab 1/p            \tab 1/(1-p)            \cr
#' Average causal effect in exposed                                         \tab \code{psw_ACE_Exposed}   \tab 1              \tab p/(1-p)            \cr
#' Average causal effect in unexposed                                       \tab \code{psw_ACE_Unexposed} \tab (1-p)/p        \tab 1                  \cr
#' Average causes effect in population for which sample is most informative \tab \code{psw_ACE_MostInfo}  \tab 1-p            \tab p                  \cr
#' Average causal effect in mathcin weight population                       \tab \code{psw_ACE_MWP}       \tab min(p,(1-p))/p \tab min(p,(1-p))/(1-p) \cr
#' }
#'
#' @param .data a \code{data.frame}
#' @param exposure the bare name for the column within \code{.data} indicating
#' the exposure or non-exposure.  Expect this to be an integer, or coercible to
#' an integer with values of zero or one.
#' @param ps the propensity scores.  Expected values between 0 and 1.
#'
#' @return a \code{data.frame} with the exposure and ps vectors returned along
#' with five different weights.  See Details for information on the five
#' weights.
#'
#' @examples
#' glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                      data = pride,
#'                      family = stats::binomial()) 
#'
#' ourpride <- dplyr::mutate(pride, pp = qwraps2::invlogit(fitted(glmfit))) 
#'
#' psw <- psweights(dplyr::tbl_df(ourpride), PCR_RSV, pp)
#' psw
#' 
#' # Mirrored Histograms
#' plot(psw)
#' 
#' @references
#' Li, Liang, and Tom Greene. "A weighting analogue to pair matching in
#' propensity score analysis." The international journal of biostatistics 9.2
#' (2013): 215-234.
#'
#' @export
psweights <- function(.data, exposure, ps) {
  psweights_(.data, deparse(substitute(exposure)), deparse(substitute(ps)))
}

#' @param exposure_col a character string
#' @param ps_col a character string
#' @export
#' @rdname psweights
psweights_ <- function(.data, exposure_col, ps_col) {
  UseMethod("psweights_")
}

#' @export
psweights_.grouped_df <- function(.data, exposure_col, ps_col) {
  warning("grouped_df will be ungrouped and a row-wise operation will take place.", call. = FALSE)
  psweights_(dplyr::ungroup(.data), deparse(substitute(exposure_col)), deparse(substitute(ps)))
}

#' @export
psweights_.rowwise_df <- function(.data, exposure_col, ps_col) {
  warning("rowwise_df will be ungrouped.", call. = FALSE)
  psweights_(dplyr::ungroup(.data), deparse(substitute(exposure_col)), deparse(substitute(ps)))
}


#' @export
psweights_.data.frame <- function(.data, exposure_col, ps_col) {
  
  if (!all(.data[[exposure_col]] %in% c(0, 1))) {
    stop("Exposure vector needs to be 0/1", call. = FALSE)
  }

  if (any(.data[[ps_col]] < 0 | .data[[ps_col]] > 1)) {
    stop("Propensity score < 0 or > 1 detected.", call. = FALSE)
  }

  if (any(.data[[ps_col]] == 0 | .data[[ps_col]] == 1)) {
    warning("Propensity score of 0 or 1 detected.", call. = FALSE)
  }

  out <-dplyr::rowwise(.data)
  out <- dplyr::mutate_(out,
                        .dots = list(
                                     "psw_IPW"           = lazyeval::interp( ~ 1/(z*ps + (1-z)*(1-ps)),                  ps = as.name(ps_col), z = as.name(exposure_col)),
                                     "psw_ACE_Exposed"   = lazyeval::interp( ~ z + (1-z)*ps/(1-ps),                      ps = as.name(ps_col), z = as.name(exposure_col)),
                                     "psw_ACE_Unexposed" = lazyeval::interp( ~ z * (1-ps)/ps + (1-z),                    ps = as.name(ps_col), z = as.name(exposure_col)),
                                     "psw_ACE_MostInfo"  = lazyeval::interp( ~ z * (1-ps) + (1-z) * ps,                  ps = as.name(ps_col), z = as.name(exposure_col)),
                                     "psw_ACE_MWP"       = lazyeval::interp( ~ min(c(ps, 1-ps)) / (z*ps + (1-z)*(1-ps)), ps = as.name(ps_col), z = as.name(exposure_col))
                                     )
                        )
  out <- dplyr::ungroup(out)
  attr(out, "exposure_col") <- exposure_col
  attr(out, "ps_col") <- ps_col
  class(out) <- c("pstools_psweights", class(out))
  out 
}


#' @export
#' @param x a \code{pstools_psweights} object
#' @param ... ignored
#' @describeIn psweights Mirrored histograms
plot.pstools_psweights <- function(x, ...) {

  dat <- tidyr::gather_(x, 
                        key_col = "method",
                        value_col = "value",
                        gather_cols = c("psw_IPW",
                                        "psw_ACE_Exposed",
                                        "psw_ACE_Unexposed",
                                        "psw_ACE_MostInfo",
                                        "psw_ACE_MWP"),
                        factor_key = TRUE) 

  ggplot2::ggplot() +
  ggplot2::aes_string(x = attr(x, "ps_col")) +

  ggplot2::geom_histogram(data = dplyr::filter_(x, .dots = lazyeval::interp( ~ z == 1, z = as.name(attr(x, "exposure_col"))))) +
  ggplot2::geom_histogram(data = dplyr::filter_(x, .dots = lazyeval::interp( ~ z == 0, z = as.name(attr(x, "exposure_col")))),
                          mapping = ggplot2::aes_string(y = "-..count..")) +

  ggplot2::geom_histogram(data = dplyr::filter_(dat, .dots = lazyeval::interp( ~ z == 1, z = as.name(attr(x, "exposure_col")))),
                          mapping = ggplot2::aes_string(weight = "value", fill = "method"),
                          alpha = 0.5
                          ) +
  ggplot2::geom_histogram(data = dplyr::filter_(dat, .dots = lazyeval::interp( ~ z == 0, z = as.name(attr(x, "exposure_col")))),
                          mapping = ggplot2::aes_string(weight = "-value", fill = "method"),
                          alpha = 0.5
                          ) +
  ggplot2::facet_wrap( ~ method, nrow = 1)
                     

}
