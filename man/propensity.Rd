% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propensity.R
\name{propensity}
\alias{propensity}
\title{Propensity Summary}
\usage{
propensity(fit, weight_method = 1)
}
\arguments{
\item{fit}{a regression object such as}

\item{weight_method}{one of three methods for defining the weight for each
row of the data set.  See Details.}
}
\value{
A \code{pstools_propensity} object which is a
\code{data.frame} with columns summarizing each varaible used as a predictor in the
propensity model. The function will determine if each variable is a
categorical or continous variable. if the variable is continuous the mean (sd)
is returned.  If the variable is categorical the the percentage of each level
is returned.  Standard differences are reported as well. Unadjusted and
Matched Weight (see Li and Greene (2012)) values are reported.
}
\description{
Given a logistic regression model, generate the propensity weights for the
subjects (rows) in the data set.
}
\details{
The regression model is expected to estimate the probability of an exposure
(Z = 1) given a set of predictors, X, i.e., Pr[Z = 1 | X].

NOTE: for binary predictors coded as 0/1, such as male, the default action of
propensity_summary will return a mean (sd), that is, if
propensity_summary(geeglm(y ~ male)) will summarize the variable male as if it
was continous predictor. To get the summary for the proportion of males, and
females will be reported too, use propensity_summary(geeglm(y ~
factor(male))).

\code{weight_method} controls the weight given to each row of the data set.
Let ps = Pr[Z = 1 | X].  The default setting is \code{weight_method = 1}.

\describe{
  \item{\code{weight_method = 1}}{\code{w =  min(c(ps, 1-ps)) / (z*ps + (1-z)*(1-ps))}}
  \item{\code{weight_method = 2}}{\code{w = 1/(z*ps + (1-z)*(1-ps))}}
  \item{\code{weight_method = 3}}{\code{w = z + (1-z)*ps/(1-ps))}}
}
}
\examples{
\dontrun{
data(pride)
glmfit <- stats::glm(PCR_RSV ~ SEX + RSVINF + REGION + AGE + ELTATOP + EINZ + EXT, 
                     data = pride,
                     family = stats::binomial())
propensity(glmfit)
summary(propensity(glmfit))
plot(propensity(glmfit))
}
}
\references{
Li, Liang, and Tom Greene. "A weighting analogue to pair matching in
propensity score analysis." The international journal of biostatistics 9.2
(2013): 215-234.
}
\seealso{
\code{\link[stats]{glm}} for fitting logistic regression models, or 
\code{\link[geepack]{geeglm}} for fitting GEEs.
}

