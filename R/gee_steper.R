#' Stepwise GEE Model Selection
#'
#' Starting with a model of y ~ 1, use a forward/backward step selection process
#' to find the preferable fit.
#'
#' @param fit a \code{geeglm} object.  the response, id, data, corstr, ... will
#' be reused in all the the fits called by \code{gee_stepper}
#' @param upper A formula with the rhs containing all the variables that *could*
#' be the the model.
#'
#' @return the "best" model based on QIC from the stepwise selection process.
#'
#' @examples
#' 
#' data(pride)
#' library(geepack)
#' geefit <- geepack::geeglm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                   id = REGION, data = pride, family = stats::binomial())
#' gee_stepper(geefit, formula(geefit))
#' 
#' @export
gee_stepper <- function(fit, upper) {
  UseMethod("gee_stepper")
}

#' @export
gee_stepper.geeglm <- function(fit, upper) {
  preds     <- all.vars(as.list(upper)[[3]])
  preds_out <- preds
  preds_in  <- character(0)

  fit0 <- stats::update(fit, formula = . ~ 1)

  while(TRUE) {
    fits <-
      c(list(fit0),
        lapply(preds_out, function(p) {stats::update(fit0, formula = stats::as.formula(paste(". ~ . +", p))) }),
        lapply(preds_in,  function(p) {stats::update(fit0, formula = stats::as.formula(paste(". ~ . -", p))) })
        )

    minqic <- which.min(sapply(fits, MESS::QIC)[1, ])
    fit0 <- fits[[minqic]]

    status <-
      dplyr::data_frame(preds = sapply(fits, function(f) {as.character(stats::formula(f))[3]}),
                        qic   = sapply(fits, MESS::QIC)[1, ],
                        pick  = seq_along(fits) == minqic)

    preds_in  <- all.vars(as.list(stats::formula(fit0))[[3]])
    preds_out <- dplyr::setdiff(preds, preds_in)

    print(status)

    if (minqic == 1) {
      break
    }
  }
  fit0
}

