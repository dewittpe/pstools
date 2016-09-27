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

test_that("means", 
          {
            fit <- glm(mpg > 20 ~ wt + cyl, data = mtcars, family = binomial())
            obj <- propensity_summary(fit)
            obj

            expect_equal(obj$key, c("cyl", "wt"))
            expect_equivalent(unclass(obj[2, c(2, 6)]), 
                              aggregate(wt ~ I(mpg > 20), data = mtcars, FUN = mean)[, 2]
                              )
          })
