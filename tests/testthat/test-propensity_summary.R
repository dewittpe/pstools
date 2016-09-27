test_that("model family is detected",
          {
            fit <- glm(mpg > 20 ~ wt + cyl, data = mtcars)
            expect_error(propensity_summary(fit), "binomial") 
          })

test_that("exposure is 0/1",
          { 
            fit <- suppressWarnings(glm(y ~ x, data = list(x = 1:10, y = runif(10)), family = binomial()))
            expect_error(propensity_summary(fit), "0/1 integer vector")
          })

test_that(
