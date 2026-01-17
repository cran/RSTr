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
# mod_mst <- mstcar(name = "my_test_model", data = miheart, adjacency = miadj, seed = 1234)

## ----eval = is_cran-----------------------------------------------------------
# For computational reasons, full model fitting is not run during CRAN checks.
# When building on CRAN, this vignette loads a pre-fitted example model included with the package.
# The pkgdown website shows the full model-fitting workflow.
example_dir <- system.file("extdata", package = "RSTr")
mod_mst <- load_model("mstcar_example", example_dir)

## -----------------------------------------------------------------------------
samples <- load_samples(mod_mst, param = "lambda", burn = 2000) * 1e5

## -----------------------------------------------------------------------------
dim(samples)

## -----------------------------------------------------------------------------
margin_time <- 3
pop <- mod_mst$data$n
samples_7988 <- aggregate_samples(samples, pop, margin_time)

## -----------------------------------------------------------------------------
samples <- aggregate_samples(samples, pop, margin_time, bind_new = TRUE, new_name = "1979-1988")

## -----------------------------------------------------------------------------
age <- c("35-44", "45-54", "55-64")
std_pop <- c(113154, 100640, 95799)

## -----------------------------------------------------------------------------
dim(samples)

## -----------------------------------------------------------------------------
margin_age <- 2
groups <- c("35-44", "45-54", "55-64")
samples_3564 <- standardize_samples(samples, std_pop, margin_age, groups)

## -----------------------------------------------------------------------------
samples <- standardize_samples(
  samples,
  std_pop,
  margin_age,
  groups,
  bind_new = TRUE,
  new_name = "35-64"
)

## -----------------------------------------------------------------------------
medians <- get_medians(samples)

## -----------------------------------------------------------------------------
ci <- get_credible_interval(sample = samples, perc_ci = 0.95)
rel_prec <- get_relative_precision(medians, ci)
low_rel_prec <- rel_prec < 1

## -----------------------------------------------------------------------------
pop <- aggregate_count(pop, margin_age, groups = 1:3, bind_new = TRUE, new_name = "35-64")
pop <- aggregate_count(pop, margin_time, bind_new = TRUE, new_name = "1988-1988")
low_population <- pop < 1000
medians_supp <- medians
medians_supp[low_rel_prec | low_population] <- NA

## ----eval = !is_cran----------------------------------------------------------
# est_3544 <- medians_supp[, "35-44", "1988"]
# 
# ggplot(mishp) +
#   geom_sf(aes(fill = est_3544)) +
#   labs(
#     title = "Smoothed Myocardial Infarction Death Rates in MI, Ages 35-44, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# Mapping examples are shown on the RSTr pkgdown website.
NULL

## ----eval = !is_cran----------------------------------------------------------
# ci <- get_credible_interval(samples, 0.995)
# rel_prec50 <- get_relative_precision(medians, ci)
# low_rel_prec <- rel_prec50 < 1
# medians_supp <- medians
# medians_supp[low_rel_prec | low_population] <- NA
# 
# est_3544 <- medians_supp[, "35-44", "1988"]
# 
# ggplot(mishp) +
#   geom_sf(aes(fill = est_3544)) +
#   labs(
#     title = "Smoothed Myocardial Infarction Death Rates in MI, 99.5% CI, Ages 35-44, 1988",
#     fill = "Deaths per 100,000"
#   ) +
#   scale_fill_viridis_c() +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# Mapping examples are shown on the RSTr pkgdown website.
NULL

## ----eval = !is_cran----------------------------------------------------------
# crude_3544 <- sum(mod_mst$data$Y[, "35-44", "1988"]) / sum(mod_mst$data$n[, "35-44", "1988"]) * 1e5
# sample_3544 <- samples[, "35-44", "1988", ]
# p_higher <- apply(sample_3544, 1, \(county) mean(county > crude_3544)) * 100
# 
# ggplot(mishp) +
#   geom_sf(aes(fill = p_higher)) +
#   labs(
#     title = "Probability that County Rate > State Rate MI, Ages 35-44, 1988",
#     fill = "Probability"
#   ) +
#   scale_fill_continuous(palette = "RdBu", trans = "reverse") +
#   theme_void()

## ----eval = is_cran, echo = FALSE---------------------------------------------
# Mapping examples are shown on the RSTr pkgdown website.
NULL

