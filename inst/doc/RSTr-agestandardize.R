## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(RSTr)
is_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")


## ----eval = !is_cran, results = "hide", fig.keep = "last"---------------------
# mod_mst <- mstcar(name = "my_test_model", data = miheart, adjacency = miadj)

## ----eval = is_cran, echo = FALSE---------------------------------------------
# For computational reasons, full model fitting is not run during CRAN checks.
# When building on CRAN, this vignette loads a pre-fitted example model included with the package.
# The pkgdown website shows the full model-fitting workflow.
example_dir <- system.file("extdata", package = "RSTr")
mod_mst <- load_model("mstcar_example", example_dir)

## -----------------------------------------------------------------------------
estimates <- get_estimates(mod_mst, rates_per = 1e5)
head(estimates)

## -----------------------------------------------------------------------------
std_pop <- c(113154, 100640, 95799)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35-64", groups = c("35-44", "45-54", "55-64"))
mod_mst

## -----------------------------------------------------------------------------
std_pop <- c(68775, 34116, 9888)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "65up", groups = c("65-74", "75-84", "85+"))
mod_mst

## -----------------------------------------------------------------------------
std_pop <- c(113154, 100640, 95799, 68775, 34116, 9888)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35up")
mod_mst
mst_estimates_as <- get_estimates(mod_mst)
head(mst_estimates_as)

