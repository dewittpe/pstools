# pstools

Tools for Propensity Score Analyses

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/dewittpe/pstools.svg?branch=master)](https://travis-ci.org/dewittpe/pstools)
[![Coverage Status](https://img.shields.io/codecov/c/github/dewittpe/pstools/master.svg)](https://codecov.io/github/dewittpe/pstools?branch=master)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/pstools)](https://cran.r-project.org/package=pstools)

[![License](https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.0.2-6666ff.svg)](https://cran.r-project.org/)

## Examples

## Install

### From CRAN
This package is currently *not* on CRAN.  Development is moving forward.

### Developmental
Install the development version of `pstools` directly from github via the 
[`devtools`](https://github.com/hadley/devtools/) package:

    if (!("devtools" %in% rownames(installed.packages()))) { 
      warning("installing devtools from https://cran.rstudio.com")
      install.packages("devtools", repo = "https://cran.rstudio.com")
    }

    devtools::install_github("dewittpe/pstools", build_vignettes = TRUE)

*NOTE:* If you are working on a Windows machine you will need to download and
install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/) before
`devtools` will work for you.

