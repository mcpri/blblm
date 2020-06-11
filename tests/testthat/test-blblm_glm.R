test_that("blblm_glm", {
  fit<-blblm_glm(mpg~wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)

  expect_s3_class("blblm_glm")
})