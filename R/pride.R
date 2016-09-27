#' Pediatric Respiratory Infection, Deutschland (pride)
#'
#' Data dealing with lower respiratory tract infections (LRTI) in n=3.078
#' infants and children aged less than three years in the observational study
#' Pri.DE (Pediatric Respiratory Infection, Deutschland) in Germany.
#'
#' @source This data set was aquireed from the \code{nonrandom} R package.
#' https://cran.r-project.org/web/packages/nonrandom/
#'
#' @format a data.frame with n = 3078 observations on 15 variables.
#' \describe{
#'   \item{PCR_RSV}{current respiratory syncytial virus (RSV) infection (1) or not (0)}
#'   \item{SEX}{gender: male (1) and female (2)}
#'   \item{ETHNO}{ethnic group: German (1), European Union (2) and others (3)}
#'   \item{FRUEHG}{preterm delivery (1) or not (0)}
#'   \item{RSVINF}{former RSV infection (1) or not (0)}
#'   \item{HERZ}{congenital heart defect (1) or not (0)}
#'   \item{REGION}{German region: South (1), East (2), West (3) and North (4)}
#'   \item{AGE}{age in years}
#'   \item{VOLLSTIL}{current breast feeding or longer than two months (1) or not (0)}
#'   \item{EINZ}{siblings (1) or not (0)}
#'   \item{TOBACCO}{passive tobacco smoke exposure at home (1) or not (0)}
#'   \item{EXT}{external care (1) or not (0)}
#'   \item{ELTATOP}{parental atopy (1) or not (0)}
#'   \item{SEVERE}{severe LRTI (1) or not (0)}
#'   \item{KRANKSUM}{number of diagnosed LRTI}
#' }
#'
#' @examples
#' library(ggplot2)
#' data(pride, package = "icptools")
#' 
#' glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
#'                      data = pride, family = stats::binomial())
#' plotdata <- cbind(pride, propensity = qwraps2::invlogit(fitted(glmfit)))
#' 
#' ggplot2::ggplot(plotdata) + 
#' ggplot2::aes_string(x = "propensity", fill = "factor(PCR_RSV)") + 
#' ggplot2::geom_density(alpha = 0.5)
"pride" 
