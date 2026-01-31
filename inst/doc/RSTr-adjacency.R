## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)
library(RSTr)
if (!requireNamespace("sf", quietly = TRUE) ||
    !requireNamespace("spdep", quietly = TRUE)) {
  knitr::knit_exit()
}
is_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")


## ----results = "hide"---------------------------------------------------------
head(mamap)

## ----eval = !is_cran----------------------------------------------------------
# ma_shp <- sf::st_as_sf(mamap[order(mamap$GEOID), ])
# ma_adj <- spdep::poly2nb(ma_shp)

## ----eval = is_cran, echo = FALSE---------------------------------------------
# The lines above cause segfault issues when trying to run checks on clang-san. This chunk will load in ma_adj identically to above.
ma_adj <- readRDS(system.file("extdata", "ma_adj.Rds", package = "RSTr"))
ma_shp <- mamap

## -----------------------------------------------------------------------------
ma_adj

## -----------------------------------------------------------------------------
no_neigh <- spdep::card(ma_adj) == 0

## -----------------------------------------------------------------------------
ma_shp[no_neigh, ]

## ----warning = FALSE, eval = !is_cran-----------------------------------------
# ggplot(mamap) +
#   geom_sf(aes(fill = NAME)) +
#   geom_sf_label(aes(label = NAME))

## ----eval = is_cran, echo = FALSE---------------------------------------------
# This map is shown on the pkgdown website, but skipped for CRAN checks
NULL

## -----------------------------------------------------------------------------
county_key <- seq_along(ma_shp$NAME)
names(county_key) <- ma_shp$NAME
county_key

## -----------------------------------------------------------------------------
ma_adj <- add_neighbors(ma_adj, c(1, 4, 10))

