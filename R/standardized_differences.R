#' Standardized Differences
#' 
#' Calculate the standardized differences (percentages).
#'
#' @param x a logistic regression object, or a \code{pstools_dstats} object.
#' @param ... ignored at this time.
#'
#' @return A list of \code{data.frame}s each element of the list is for a
#' different weighting method.
#' 
#' @seealso \code{\link{psweights}}, \code{\link{dstats}}
#'
#' @examples
#' data(pride)
#' glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                      data = pride,
#'                      family = stats::binomial())
#' 
#' standardized_differences(glmfit)
#' @export
standardized_differences <- function(x, ...) {
  UseMethod("standardized_differences")
}

#' @export
standardized_differences.glm <- function(x, ...) {
  standardized_differences(dstats(x))
}

#' @export
standardized_differences.pstools_dstats <- function(x, ...) {
  sets <-
    list(dplyr::select_(x, ~ key, ~ z, ~ unadj_mean, ~ unadj_var),
         dplyr::select_(x, ~ key, ~ z, ~ adj_mean_IPW, ~ adj_var_IPW),
         dplyr::select_(x, ~ key, ~ z, ~ adj_mean_ACE_Exposed, ~ adj_var_ACE_Exposed),
         dplyr::select_(x, ~ key, ~ z, ~ adj_mean_ACE_Unexposed, ~ adj_var_ACE_Unexposed),
         dplyr::select_(x, ~ key, ~ z, ~ adj_mean_ACE_MostInfo, ~ adj_var_ACE_MostInfo),
         dplyr::select_(x, ~ key, ~ z, ~ adj_mean_ACE_MWP, ~ adj_var_ACE_MWP)
         )
  
  std_diffs <-
    lapply(sets,
           function(x) { 
             dplyr::summarize_(dplyr::group_by_(x, .dots = ~ key),
                               .dots = list("std_diff" = lazyeval::interp( ~ 100 * diff(w1) / sqrt(sum(w)/2), 
                                                                          w1 = as.name(names(x)[3]), w = as.name(names(x)[4])) 
                                            )) 
           })

  sets <-
    lapply(sets,
           function(x) {
             rtn <-
               dplyr::left_join(dplyr::filter_(x, .dots = ~ z == 0),
                                dplyr::filter_(x, .dots = ~ z == 1),
                                by = "key",
                                suffix = c(".z0", ".z1"))
             names(rtn) <- gsub(".*mean.*(\\.z\\d)", "mean\\1", names(rtn))
             names(rtn) <- gsub(".*var.*(\\.z\\d)", "var\\1", names(rtn))
             dplyr::select_(rtn, .dots = list( ~ - z.z0, ~ - z.z1))
           })

  out <-
    mapply(dplyr::left_join,
           x = sets, 
           y = std_diffs, 
           by = "key", 
           SIMPLIFY = FALSE)

  attr(out, "continuous") <- attr(x, "continuous")
  class(out) <- c("pstools_std_diffs", class(out))
  out
}

