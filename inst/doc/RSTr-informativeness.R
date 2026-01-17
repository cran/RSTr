## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
is_cran <- identical(Sys.getenv("NOT_CRAN"), "false")
if (Sys.getenv("NOT_CRAN") == "") is_cran <- TRUE
library(RSTr)
library(ggplot2)

## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# data_u <- lapply(miheart, \(x) x[, "55-64", "1979", drop = FALSE])
# mod_car <- car("my_test_model", data_u, miadj, tempdir(), seed = 1234)

## ----eval = is_cran-----------------------------------------------------------
# For computational reasons, full model fitting is not run during CRAN checks.
# When building on CRAN, this vignette loads a pre-fitted example model included with the package.
# The pkgdown website shows the full model-fitting workflow.
example_dir <- system.file("extdata", package = "RSTr")
mod_car <- load_model("car_example", example_dir)

## -----------------------------------------------------------------------------
estimates <- get_estimates(mod_car)
estimates_supp <- estimates[estimates$rel_prec > 1 & estimates$events < 10, ]
plot(estimates$events, estimates$rel_prec, xlab = "Events", ylab = "Relative Precision")
points(estimates_supp$events, estimates_supp$rel_prec, col = "red")
abline(h = 1, col = "blue")
abline(v = 10, col = "blue")

## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# mod_rcar <- rcar("my_test_model", data_u, miadj, tempdir(), seed = 1234, A = 6)

## ----eval = is_cran-----------------------------------------------------------
# Same as above, but for RCAR model
example_dir <- system.file("extdata", package = "RSTr")
mod_rcar <- load_model("rcar_example", example_dir)

## -----------------------------------------------------------------------------
estimates_rcar <- get_estimates(mod_rcar)
plot(estimates_rcar$events, estimates_rcar$rel_prec, xlab = "Events", ylab = "Relative Precision", col = "purple")
points(estimates$events, estimates$rel_prec)
abline(h = 1, col = "blue")
abline(v = 10, col = "blue")

## ----eval = !is_cran----------------------------------------------------------
# ggplot(mishp) +
#   geom_sf(aes(fill = estimates$medians)) +
#   labs(
#     title = "Spatially Smoothed Estimates, Unrestricted CAR Model",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()
# ggplot(mishp) +
#   geom_sf(aes(fill = estimates_rcar$medians)) +
#   labs(
#     title = "Spatially Smoothed Estimates, RCAR Model",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# data_u <- lapply(miheart, \(x) x[, c("65-74", "75-84", "85+"), "1988", drop = FALSE])
# A <- 6 * colSums(data_u$Y) / sum(data_u$Y)
# mod_rcar <- rcar("test_rcar", data_u, miadj, tempdir(), seed = 1234, A = A)

## ----eval = !is_cran----------------------------------------------------------
# std_pop <- c(68775, 34116, 9888)
# mod_rcar <- age_standardize(mod_rcar, std_pop, "65up")
# est_rcar <- get_estimates(mod_rcar)
# ggplot(mishp) +
#   geom_sf(aes(fill = est_rcar$medians)) +
#   labs(
#     title = "Age-Standardized Spatially Smoothed Estimates, RCAR Model",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

## ----eval = !is_cran----------------------------------------------------------
# mod_rcar <- suppress_estimates(mod_rcar)
# est_rcar <- get_estimates(mod_rcar)
# ggplot(mishp) +
#   geom_sf(aes(fill = est_rcar$medians_suppressed)) +
#   labs(
#     title = "Spatially Smoothed Estimates, RCAR Model",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

