test_that("blblm", {
  fit<-blblm(mpg~wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)

  expect_s3_class("blblm")
})
