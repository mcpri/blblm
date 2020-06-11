test_that("blblm_par_test", {
  fit<-blblm_par(mpg~wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE, cl = 3)

  expect_s3_class("blblm_par")
})