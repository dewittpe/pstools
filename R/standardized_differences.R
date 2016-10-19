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
  names(out) <- c("unadj", "adj_IPW", "adj_ACE_Exposed", 
                  "adj_ACE_Unexposed", "adj_ACE_MostInfo", "adj_ACE_MWP")

  attr(out, "continuous") <- attr(x, "continuous")
  class(out) <- c("pstools_std_diffs", class(out))
  out
}

frmt_pstools_std_diffs <- function(x, ...) {
  mean_sd <- function(x, y) {
    paste0(qwraps2::frmt(x), " (", qwraps2::frmt(y), ")")
  }

  percent <- function(x) {
    paste0(qwraps2::frmt(x * 100, 1), "%")
  }

  x <- lapply(x, dplyr::mutate_, .dots = list("cnt" = ~ attr(x, "continuous")))

  x <- 
    mapply(function(dat, cnt) {
             dat <- dplyr::rowwise(dat)
             dat <- dplyr::transmute_(dat,
                                      .dots = list( 
                                                   key = ~ key,
                                                   unexposed = ~ if (cnt) { mean_sd(mean.z0, sqrt(var.z0)) } else { percent(mean.z0) },
                                                   exposed = ~ if (cnt) { mean_sd(mean.z1, sqrt(var.z1)) } else { percent(mean.z1) },
                                                   std_diff = ~ paste0(qwraps2::frmt(std_diff, 2), "%")
                                                   )
                                      )
             dat <- dplyr::ungroup(dat)
                  },
           dat = x, 
           MoreArgs = list(attr(x, "continuous")), 
           SIMPLIFY = FALSE)
}

#' @export
print.pstools_std_diffs <- function(x, ...) {
  x <- frmt_pstools_std_diffs(x, ...) 
  print(x)
}

#' @export
plot.pstools_std_diffs <- function(x, adj_methods = c("adj_IPW", "adj_ACE_Exposed", "adj_ACE_Unexposed", "adj_ACE_MostInfo", "adj_ACE_MWP"), ...) { 
  dat <- dplyr::bind_rows(x, .id = "adjustment")
  dat <- dplyr::filter_(dat, .dots = ~ adjustment %in% c("unadj", adj_methods))
  dat$key <- factor(dat$key)
  levels(dat$key) <- rev(levels(dat$key))

  ggplot2::ggplot(dat) +
  ggplot2::aes_string(x = "std_diff", y = "key", 
                      color = "adjustment", shape = "adjustment") +
  ggplot2::geom_point() +
  ggplot2::geom_vline(xintercept = 0)
}
