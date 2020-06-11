#' @import purrr
#' @import stats
#' @import furrr
#' @import parallel
#' @import tidyverse
#' @import future
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Creating the blblm fucntion for Linear regression
#'
#' This function computes the bootstrap sample without parallelization
#' @param formula in the model
#' @param B number for bootstrap
#' @param m a subsample number
#' @param data part of the paramterst
#' @return  model with cluster and blb techniques
#' @export
#' @examples
#'fit<-blblm(mpg~wt * hp, data = mtcars, m = 3, B = 100)
#'coef(fit)
#'confint(fit, c("wt", "hp"))
blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
fit <- blblm(mpg~wt * hp, data = mtcars, m = 3, B = 100)



#' Creating the blblm parallel fucntion for Linear regression
#'
#' This function computes the bootstrap sample with parallelization
#' @param formula in the model
#' @param B number for bootstrap
#' @param m a subsample number
#' @param data part of the dataset
#' @param cl the clusters
#' @return model with cluster and blb techniques
#' @export
#' @examples
#' fit2 <- blblm_par(mpg~wt * hp, data = mtcars, m = 3, B = 100,cl = 3)
#' coef(fit2)
#' confint(fit2,c("wt", "hp"))
blblm_par <- function(formula,data,m,B,cl){
  data_list <- split_data(data,m)
  suppressWarnings(plan(multiprocess, workers = cl))
  estimates <-future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <-"blblm_par"
  invisible(res)
}
fit2 <- blblm_par(mpg~wt * hp, data = mtcars, m = 3, B = 100,cl = 3)


#' Creating the glm for blb fucntion for Linear regression
#'
#' This function computes the bootstrap samplefor glm
#' @param formula in the model
#' @param B number for bootstrap
#' @param m a subsample number
#' @param data part of the dataset
#' @param parallel  Boolean
#' @return model with cluster and blb techniques
#' @export
blblm_glm <- function(formula, data, m = 10, B = 200, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (isTRUE(parallel)){
    plan(multiprocess, workers = 4)
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  } else {
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  m = m
  B = B
  parallel = parallel
  call <- match.call()
  res <- list(estimates = estimates, formula = formula, m=m, B=B, call = call, parallel = parallel)
  class(res) <- "blblm_glm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#' @param m a integer
#' @param data part of the dataframe
#' @return  model with cluster and blb techniques
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @param formula in the model
#' @param B number for bootstrap
#' @param data the data for part of the paramters
#' @param n subsample numeric number
#' @return  model with cluster and blb techniques
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @param formula in the model
#' @param data parnt of the paramters
#' @param n subsample number
#' @return  model with cluster and blb techniques
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula in the model
#' @param data part of the parameters
#' @param freqs part of the paramters
#' @return  model with cluster and blb techniques
#' @export
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit  in the parameters
#' @return  model with cluster
#' @export
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit in the parameters
#' @return  model with cluster and blb techniques
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' Creating the Signa
#'
#' This function computes the boothstrap coefficient
#'
#' @param object the blblm interval
#' @param level for confidence
#' @param confidence false
#' @param ... continuation
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Creating the Coefficient
#'
#' This function computes the boothstrap coefficient
#' @param ... continuation
#' @param object the blblm interval
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Creating Confint Interval
#'
#' This function computes the boothstrap confidence Interval
#' @param object the blblm interval
#' @param level for confidence interval
#' @param parm null and factor-variables
#' @param ... continuation
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Creating Predict Interval
#'
#' This function computes the boothstrap Prediction for a linear model
#'
#' @param object the blblm interval
#' @param level level of confidence
#' @param confidence of the interval
#' @param new_data part of the data
#' @param ... continuation
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}
map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
