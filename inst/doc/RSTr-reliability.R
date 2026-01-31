## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
is_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")

library(RSTr)

## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# mod_mst <- mstcar(name = "my_test_model", data = miheart, adjacency = miadj, seed = 1234, perc_ci = 0.95)

## ----eval = is_cran, echo = FALSE---------------------------------------------
# For computational reasons, full model fitting is not run during CRAN checks.
# When building on CRAN, this vignette loads a pre-fitted example model included with the package.
# The pkgdown website shows the full model-fitting workflow.
example_dir <- system.file("extdata", package = "RSTr")
mod_mst <- load_model("mstcar_example", example_dir)

## -----------------------------------------------------------------------------
mod_mst <- suppress_estimates(mod_mst, threshold = 1e3)
mod_mst
estimates <- get_estimates(mod_mst)
head(estimates)

## ----eval = !is_cran----------------------------------------------------------
# library(ggplot2)
# est_3544 <- estimates$medians_suppressed[estimates$group == "35-44" & estimates$year == "1988"]
# ggplot(mishp) +
#   geom_sf(aes(fill = est_3544)) +
#   labs(
#     title = "Smoothed Myocardial Infarction Death Rates in MI, Ages 35-44, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

## -----------------------------------------------------------------------------
std_pop <- c(113154, 100640, 95799)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35-64", groups = c("35-44", "45-54", "55-64"))
estimates <- get_estimates(mod_mst)

## ----eval = !is_cran----------------------------------------------------------
# est_3564 <- estimates$medians_suppressed[estimates$group == "35-64" & estimates$year == "1988"]
# ggplot(mishp) +
#   geom_sf(aes(fill = est_3564)) +
#   labs(
#     title = "Smoothed Myocardial Infarction Death Rates in MI, Ages 35-64, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

