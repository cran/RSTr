## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)
library(RSTr)
library(ggplot2)
is_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")


## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# mod_mst <- mstcar(
#   name = "my_test_model",
#   data = miheart,
#   adjacency = miadj,
#   dir = tempdir(),
#   seed = 1234
# )

## ----eval = is_cran, echo = FALSE---------------------------------------------
# For computational reasons, full model fitting is not run during CRAN checks.
# When building on CRAN, this vignette loads a pre-fitted example model included with the package.
# The pkgdown website shows the full model-fitting workflow.
example_dir <- system.file("extdata", package = "RSTr")
mod_mst <- load_model("mstcar_example", example_dir)

## -----------------------------------------------------------------------------
mod_mst

## -----------------------------------------------------------------------------
mst_estimates <- get_estimates(mod_mst, rates_per = 1e5)
head(mst_estimates)

## -----------------------------------------------------------------------------
std_pop <- c(68775, 34116, 9888)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35-64", groups = c("35-44", "45-54", "55-64"))
mod_mst

## -----------------------------------------------------------------------------
mst_estimates_as <- get_estimates(mod_mst)
head(mst_estimates_as)

## -----------------------------------------------------------------------------
mod_mst <- suppress_estimates(mod_mst, threshold = 1e3)
mod_mst
mst_estimates_as <- get_estimates(mod_mst)
head(mst_estimates_as)

## ----eval = !is_cran----------------------------------------------------------
# # Original Myocardial Infarction Death Rates in MI, Ages 35-64, 1988
# estimates_88 <- mst_estimates_as[mst_estimates_as$year == "1988", ]
# estimates_3564 <- estimates_88[estimates_88$group == "35-64", ]
# raw_3564 <- (estimates_3564$events / estimates_3564$population * 1e5)
# ggplot(mishp) +
#   geom_sf(aes(fill = raw_3564)) +
#   labs(
#     title = "Raw Myocardial Infarction Death Rates in MI, Ages 35-64, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()
# # Spatially Smoothed MI Death Rates in MI
# est_3564 <- estimates_3564$medians
# ggplot(mishp) +
#   geom_sf(aes(fill = est_3564)) +
#   labs(
#     title = "Smoothed Myocardial Infarction Death Rates in MI, Ages 35-64, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# Mapping examples are shown on the RSTr pkgdown website.
NULL

